import numpy as np
from obspy import read, UTCDateTime, Stream
from zlib import adler32
import ctypes

from calc import compressSpectrum, smoothSpectrum, getInstrumentResponse, psdWelch
from constants import (DB_REFERENCE, MINIMUM_PERIOD, NUMBER_OF_FREQUENCIES,
                        OCTAVE_PERIOD_STEP, SEGMENT_LENGTH)
from sdsfile import SDSFile

class PSDCollector():
    """
    Calculates spectra for an ObsPy trace
    Settings are internal in calculateSpectrum
    """

    def __init__(self, connect_sql=True):
        """Initialize a PSDCollector.

        Parameters
        ----------
        connect_sql : `bool`
            Whether or not to connect to an SQL database (default `True`).

        """
        self.frequencies = self.setupFrequencies()

        if connect_sql:
            import pymysql.cursors
            self.connection = pymysql.connect(
                user="root",
                password="password",
                host="localhost",
                db="test"
            )

    def setupFrequencies(self):
        """
        def setupFrequencies
        Sets up the frequencies to calculate the PSDs for
        """

        # This is going to give us nice frequencies at 1, 2, 4Hz
        steps = np.arange(NUMBER_OF_FREQUENCIES) * OCTAVE_PERIOD_STEP
        numbers = MINIMUM_PERIOD * 2 ** steps

        return np.reciprocal(numbers)

    def calculateSpectrum(self, trace, frequencies,
                          frequency_responses=None, corresponding_frequencies=None):
        """
        def calculateSpectrum
        Calculates power spectrum of a trace using Welch's method
        Can be either psdWelch, psdWelchScipy, psdWelchMlab
        """

        fs = trace.stats.sampling_rate

        # Number of points in the Welch FFT (13 segments per hour with 75% overlap)
        # ObsPy goes to previous power of 2 (small difference)
        nfft = int(SEGMENT_LENGTH * fs / 4.0)

        # Load the instrument response from cacher if necessary
        resp = frequency_responses
        freqs = corresponding_frequencies
        if resp is None or freqs is None:
            resp, freqs = getInstrumentResponse(trace.stats)

        # Calculate the spectra: import right function from calc.py
        Pxx = psdWelch(
            trace.data,
            fs,
            resp,
            nfft=nfft,
            overlap=0.75,
            reference=DB_REFERENCE
        )

        # Smooth the spectrum over a full octave
        return smoothSpectrum(freqs, Pxx, frequencies)

    def readData(self, SDSFile):
        """
        def readData
        Reads mSEED data from disk from multiple files and creates a single stream
        """

        # Create an empty stream to fill
        ObspyStream = Stream()

        # Read neighbouring files
        for neighbour in SDSFile.neighbours:

            # Read from 0h of this day, most likely in the previous day file,
            # until half an hour in the next day [psd segment end]
            st = read(
                neighbour.filepath,
                starttime=UTCDateTime(SDSFile.start),
                endtime=UTCDateTime(SDSFile.end) + 0.5 * SEGMENT_LENGTH,
                nearest_sample=False
            )

            # Concatenate all the traces
            for tr in st:
                if tr.stats.npts != 0:
                    ObspyStream.extend([tr])

        # No data found
        if not ObspyStream:
            raise ValueError("No data for processing.")

        # Simple clean up to remove overlaps and merge traces
        # This does not fill any gaps
        ObspyStream.merge(-1)

        return ObspyStream

    def storeObjects(self, psdObjects):
        """
        def storeObject
        Stores the prepared object to the database

        # Make table
        CREATE TABLE PSD (
          network CHAR(2) NOT NULL,
          station CHAR(5) NOT NULL,
          location CHAR(2) NOT NULL,
          channel CHAR(3) NOT NULL,
          quality CHAR(1) NOT NULL,
          start DATETIME NOT NULL,
          shift SMALLINT(2) SIGNED NOT NULL,
          offset TINYINT(1) SIGNED NOT NULL,
          spectrum VARBINARY(255) NOT NULL
        ) ENGINE = InnoDB;

        # Add index
        CREATE INDEX idx_network ON PSD (network, station, location, channel, start);

        """

        try:
            with self.connection.cursor() as cursor:

                # Currently store Infrasound / Seismic in same table called PSD
                # TODO Need to store something else? E.g. data hash? Can be in another table pointing to daily 48 segments
                variables = "network, station, location, channel, quality, start, shift, offset, spectrum"
                values = "%s, %s, %s, %s, %s, %s, %s, %s, %s"

                sql = "INSERT INTO PSD (%s) VALUES (%s)" % (variables, values)

                cursor.executemany(sql, psdObjects)

            # TODO check if object already exists
            self.connection.commit()

        finally:
            self.connection.close()

    def getTrace(self, data, segmentStart):
        """
        def getTrace
        Extracts a trace from the data that fits the PSD segment
        """

        # TODO TEST TEST TEST
        # Cut trace to an hour long segment
        stream = data.slice(
            segmentStart,
            segmentStart + SEGMENT_LENGTH - 1E-6,
            nearest_sample=False
        )

        # Only a single trace may remain
        if len(stream) != 1:
            stream._cleanup(misalignment_threshold=0.5)

            # If _cleanup() doesn't fix it, throw an error
            if len(stream) != 1:
                raise ValueError("Number of traces is not one.")

        # Get the first trace
        trace = stream[0]

        # Confirm the number of points is as expected
        # Otherwise there must be a gap and thus we decide to skip the segment
        if trace.stats.npts != int(trace.stats.sampling_rate * SEGMENT_LENGTH):
            raise ValueError(
                "Number of samples in segment length is not expected.")

        return trace

    def getResponseFromFDSN(self, sds_file):
        """Get the instrumental response from FDSN.

        Taken from responseCacher/cache.py.

        Parameters
        ----------
        sds_file : `SDSFile`

        Returns
        -------
        responses : `list` of `dict`
            Each element of the list is a dictionary with the following fields:
                * "start_time" : `obspy.UTCDateTime`
                * "end_time" : `obspy.UTCDateTime`
                * "response" : `numpy.ndarray`
                * "frequencies" : `numpy.ndarray`

        """
        inventory = sds_file.inventory
        if inventory is None:
            raise Exception("Could not find inventory for this file")

        # Verify that only one station is returned
        # (different responses are modeled as separate channels with same code)
        contents = inventory.get_contents()
        if (len(contents["networks"]) != 1 or len(contents["stations"]) != 1):
            raise Exception("Inventory not valid for this file")

        # Get all responses
        response_list = []
        for channel in inventory[0][0]:
            response_dict = {}

            # Store time frame for this response
            response_dict["start_time"] = channel.start_date
            response_dict["end_time"] = channel.end_date  # TODO is end_date None or inexistent?

            # Evaluate the response to VEL for infrasound channels
            # Seismic stations should go to ACC to be consistent with the NLNM, NHNM
            output = "ACC"
            if channel.code.endswith("DF"):
                output = "VEL"

            # Call evalresp to evaluate the response
            # NFFT must be same as in Welch's method!
            # We use 13 segments with 75% overlap
            fs = channel.sample_rate
            resp, freqs = channel.response.get_evalresp_response(
                t_samp=np.reciprocal(fs),
                nfft=int(fs * SEGMENT_LENGTH / 4.0),
                output=output
            )

            response_dict["response"] = resp
            response_dict["frequencies"] = freqs

            # Drop phase imaginary information
            resp = np.abs(resp)

            response_list.append(response_dict)

        return response_list

    def process(self, SDSFile, cache_response=True):
        """Process a given SDSFile to extract PSDs.

        Parameters
        ----------
        SDSFile : `SDSFile`
            SDS file data/metadata
        cache_response : `bool`
            Whether to use the `responseCacher` (default `True`).

        Returns
        -------
        `list` of `tuple`
            Each element of the list corresponds to a segment of the PPSD computation,
            and it contains (network, station, location, channel, quality, start of segment,
            shift, offset, and the binary data).

        """
        # Load data to an ObsPy trace
        data = self.readData(SDSFile)

        # Start and end times
        segmentStart = UTCDateTime(SDSFile.start)
        segmentEnd = UTCDateTime(SDSFile.end) + 0.5 * SEGMENT_LENGTH

        psdObjects = list()

        # Get instrument response if needed
        all_responses = []
        if not cache_response:
            all_responses = self.getResponseFromFDSN(SDSFile)

        # Loop over all the hour long segments until the end is reached
        while segmentStart + SEGMENT_LENGTH <= segmentEnd:
            # Find the appropriate response for the segment time
            resp = None
            freqs = None
            for response in all_responses:
                if (response["start_time"] < segmentStart
                        and (response["end_time"] is None
                             or segmentStart + SEGMENT_LENGTH < response["end_time"])):
                    resp = response["response"]
                    freqs = response["frequencies"]

            # If no response contains segment and is not cached, skip this segment
            if not cache_response and (resp is None or freqs is None):
                segmentStart = segmentStart + 0.5 * SEGMENT_LENGTH
                continue

            # Get the trace in this time range
            # TODO a lot of testing!
            try:
                trace = self.getTrace(data, segmentStart)
            except ValueError:
                segmentStart = segmentStart + 0.5 * SEGMENT_LENGTH
                continue

            # Internal spectrum calculation: the returned spectrum is smoothed
            spectrum = self.calculateSpectrum(trace, self.frequencies, resp, freqs)

            # Compress data to a BLOB for MySQL
            offset, shift, binary = compressSpectrum(spectrum)

            # Compute a short and fast but not very safe checksum from response
            resp_checksum = adler32(resp) & 0xffffffff
            resp_checksum = ctypes.c_int32(resp_checksum).value

            # Save all metadata in a record
            psd_record = {
                "fileId": SDSFile.filename,
                "checksum": SDSFile.checksum,
                "checksumInventory": resp_checksum,
                "net": SDSFile.net,
                "sta": SDSFile.sta,
                "loc": SDSFile.loc,
                "cha": SDSFile.cha,
                "quality": SDSFile.quality,
                "ts": segmentStart.datetime,
                "te": (segmentStart + SEGMENT_LENGTH).datetime,
                "shift": shift,
                "offset": offset,
                "bin": binary
            }

            psdObjects.append(psd_record)

            # Modify the start time of the trace with 50% overlap
            segmentStart = segmentStart + 0.5 * SEGMENT_LENGTH

        # Check if at least one document has been processed
        if not psdObjects:
            raise Exception("Unable to process PSD for any segment")

        # Add checksum_prev to first segment and checksum_next to the last one,
        # no matter if they are at the beggining of the day or not. This way,
        # we can check if the file needs to be processed when previous/next
        # files are added/modified, in all possible cases.
        psdObjects[0]["checksum_prev"] = SDSFile.previous.checksum
        psdObjects[-1]["checksum_next"] = SDSFile.next.checksum

        return psdObjects

    def processAndStore(self, SDSFile):
        """Process a given SDSFile to extract PSDs and saves the data in the
        configured SQL database.

        Parameters
        ----------
        SDSFile : `SDSFile`
            SDS file data/metadata

        """
        psd_objects = self.process(SDSFile)
        psd_objects = [(rec["net"],
                        rec["sta"],
                        rec["loc"],
                        rec["cha"],
                        rec["quality"],
                        rec["ts"].isoformat(),
                        rec["shift"],
                        rec["offset"],
                        rec["bin"]) for rec in psd_objects]
        self.storeObjects(psd_objects)
