"""

  orfeus/sdsfile class for handling SDS type files.

  Author: Mathijs Koymans, 2017
  Copyright: ORFEUS Data Center, 2017
  Modified: 2019

"""

import os
import requests
import subprocess
import base64
import logging
import ctypes

from datetime import datetime, timedelta
from zlib import adler32

from obspy import read_inventory, UTCDateTime

class SDSFile():

    """
    Public Class SDSFile
    Class for handling files in SDS structure.

    Attributes
    ----------
    filename : `str`
        Name of file.
    net : `str`
        Network code.
    sta : `str`
        Station code.
    loc : `str`
        Location code.
    cha : `str`
        Channel code.
    quality : `str`
        Quality parameter.
    year : `str`
        Year in YYYY format.
    day : `str`
        Day of the year, in DDD format (i.e., it goes from "001" to "366").
    """

    # Save some configuration to the class
    fdsnws = "https://www.orfeus-eu.org/fdsnws/station/1/query/"

    def __init__(self, filename, archive_root):
        """
        Create a filestream from a given filename
        """

        try:
            # Extract stream identification
            (self.net,
             self.sta,
             self.loc,
             self.cha,
             self.quality,
             self.year,
             self.day) = filename.split(".")
        except ValueError:
            raise ValueError("Invalid SDS file submitted.")

        self.archive_root = archive_root

        # Initialize logger
        self.logger = logging.getLogger("RuleManager")

        # Initialize costly properties
        self._checksum = None
        self._inventory = None
        self._traces = None
        self._location = None

    # Returns the filename
    @property
    def filename(self):
        return self.custom_quality_filename(self.quality)

    def custom_quality_filename(self, quality):
        """Returns the filename of the file related to this one but with another quality."""
        return ".".join([self.net,
                         self.sta,
                         self.loc,
                         self.cha,
                         quality,
                         self.year,
                         self.day])

    def custom_quality_subdir(self, quality):
        """Returns the subdirectory of the archived file related to this one but
        with another quality."""
        return os.path.join(
            self.year,
            self.net,
            self.sta,
            ".".join([self.cha, quality])
        )

    # Returns custom filepath for a given file
    def custom_path(self, root):
        return os.path.join(self.custom_directory(root), self.filename)

    # Returns filepath for a given file
    @property
    def filepath(self):
        return self.custom_path(self.archive_root)

    # Returns the stream identifier
    @property
    def id(self):
        return ".".join([
            self.net,
            self.sta,
            self.loc,
            self.cha
        ])

    # Returns the subdirectory
    @property
    def sub_directory(self):
        return os.path.join(
            self.year,
            self.net,
            self.sta,
            self.channel_directory
        )

    def custom_directory(self, root):
        return os.path.join(
            root,
            self.sub_directory
        )

    # Returns the file directory based on SDS structure
    @property
    def directory(self):
        return self.custom_directory(self.archive_root)

    # Returns channel directory
    @property
    def channel_directory(self):
        return ".".join([self.cha, self.quality])

    # Returns next file in stream
    @property
    def next(self):
        return self._get_adjacent_file(1)

    # Returns previous file in stream
    @property
    def previous(self):
        return self._get_adjacent_file(-1)

    # Returns start time of file
    @property
    def start(self):
        return datetime.strptime(self.year + " " + self.day, "%Y %j")

    # Returns end time of file
    @property
    def end(self):
        return self.start + timedelta(days=1)

    # Start for dataselect pruning (start is INCLUSIVE)
    @property
    def sample_start(self):
        return self.start.strftime("%Y,%j,00,00,00.000000")

    # End for dataselect pruning (end is INCLUSIVE)
    @property
    def sample_end(self):
        return self.start.strftime("%Y,%j,23,59,59.999999")

    # Returns list of files neighbouring a file
    @property
    def neighbours(self):
        return filter(lambda x: os.path.isfile(x.filepath), [self.previous, self, self.next])

    @property
    def stats(self):
        try:
            return os.stat(self.filepath)
        except FileNotFoundError:
            return None

    def get_stat(self, enum):

        # Check if none and propagate
        if self.stats is None:
            return None

        if enum == "size":
            return self.stats.st_size
        elif enum == "created":
            return datetime.fromtimestamp(self.stats.st_ctime)
        elif enum == "modified":
            return datetime.fromtimestamp(self.stats.st_mtime)

    @property
    def size(self):
        return self.get_stat("size")

    @property
    def created(self):
        return self.get_stat("created")

    @property
    def modified(self):
        return self.get_stat("modified")

    @property
    def checksum(self):
        """
        def SDSFile::checksum
        Calculates the Adler-32 checksum for a given file, converting it to a
        signed int32 to save space in MongoDB documents
        """

        if self._checksum is not None:
            return self._checksum

        if self.stats is None:
            return None

        checksum = None
        with open(self.filepath, "rb") as f:
            checksum = adler32(f.read()) & 0xffffffff
            checksum = ctypes.c_int32(checksum).value
        self._checksum = checksum
        return self._checksum

    @property
    def query_string_txt(self):

        return self.query_string + "&" + "&".join([
            "format=text",
            "level=channel"
        ])

    @property
    def query_string_xml(self):

        return self.query_string + "&" + "&".join([
            "format=fdsnxml",
            "level=response"
        ])

    @property
    def query_string(self):
        """Return the query string for a particular SDS file."""

        return "?" + "&".join([
            "start=%s" % self.start.isoformat(),
            "end=%s" % self.end.isoformat(),
            "network=%s" % self.net,
            "station=%s" % self.sta,
            "location=%s" % self.loc,
            "channel=%s" % self.cha
        ])

    @property
    def samples(self):
        """
        def SDSFile::samples
        Returns number of samples in the SDS file
        """

        return sum(map(lambda x: x["samples"], self.traces))

    @property
    def continuous(self):
        """
        def SDSFile::continuous
        Returns True when a SDS file is considered continuous and the
        record start time & end time come before and after the file ending respectively
        """

        return len(
            self.traces) == 1 and (
            self.traces[0]["start"] <= self.start) and (
            self.traces[0]["end"] >= self.end)

    @property
    def traces(self):
        """
        def SDSFile::traces
        Returns a list of traces
        """

        if self._traces is not None:
            return self._traces

        def parse_msi_output(line):
            """Parse the MSI output."""

            # Format for seed dates e.g. 2005,068,00:00:01.000000
            SEED_DATE_FMT = "%Y,%j,%H:%M:%S.%f"

            (stream, start, end, rate, samples) = map(lambda x: x.decode("ascii"), line.split())

            # Return a simple dict with some information
            return {
                "samples": int(samples),
                "rate": float(rate),
                "start": datetime.strptime(start, SEED_DATE_FMT),
                "end": datetime.strptime(end, SEED_DATE_FMT)
            }

        # Cut to day boundary on sample level
        dataselect = subprocess.Popen([
            "dataselect",
            "-ts", self.sample_start,
            "-te", self.sample_end,
            "-Ps",
            "-szs",
            "-o", "-",
        ] + list(map(lambda x: x.filepath, self.neighbours)), stdout=subprocess.PIPE)

        lines = subprocess.check_output([
            "msi",
            "-ts", self.sample_start,
            "-te", self.sample_end,
            "-T",
            "-"
        ], stdin=dataselect.stdout, stderr=subprocess.DEVNULL).splitlines()

        # Not sure why we need this
        dataselect.stdout.close()

        # Avoid warning when status code for child process is not read (for Python 3.6):
        # Introduced in https://github.com/python/cpython/commit/5a48e21ff163107a6f1788050b2f1ffc8bce2b6d#diff-cc136486b4a8e112e64b57436a0619eb
        dataselect.wait()

        # Skip first header & final line
        self._traces = list(map(parse_msi_output, lines[1:-1]))
        return self._traces

    @property
    def is_pressure_channel(self):
        """Return true when the channel is an infrasound channel."""

        return self.cha.endswith("DF")

    @property
    def inventory(self):
        """
        def SDSFile::inventory
        Returns the FDSNWSXML inventory
        """

        if self._inventory is not None:
            return self._inventory

        # Query our FDSNWS Webservice for the station location
        request = os.path.join(self.fdsnws, self.query_string_xml)

        try:
            self._inventory = read_inventory(request)

        # Re-raise in case this is the Rule Manager timeout going off
        except TimeoutError:
            raise

        # Deal with an Exception from read_inventory
        except Exception:
            return None

        return self._inventory

    @property
    def location(self):
        """
        def SDSFile::location
        Returns the geographical location of the stream
        """

        if self._location is not None:
            return self._location

        # Query our FDSNWS Webservice for the station location
        try:
            request = requests.get(os.path.join(self.fdsnws, self.query_string_txt))
        except requests.exceptions.RequestException:
            return None

        # Any error just ignore
        if request.status_code != 200:
            return None

        lines = request.text.split("\n")

        # Multiple lines means that the location is somehow ambiguous
        if len(lines) != 3:
            return None

        # Some magic parsing: fields 4, 5, 6 on the 2nd line
        (latitude, longitude, elevation) = lines[1].split("|")[4:7]

        self._location = {
            "longitude": float(longitude),
            "latitude": float(latitude),
            "elevation": float(elevation)
        }
        return self._location

    @property
    def psd_bins(self):
        """Return 48 times starting at the start of the SDSFile with 30 minute increments."""

        return map(UTCDateTime, map(lambda x: self.start + timedelta(minutes=(30 * x)), range(48)))

    def prune(self, cut_boundaries=True, remove_overlap=False, repack=False, record_length=4096):
        """Preprocess file using IRIS dataselect and msrepack, and saves the resulting file
        with the quality indicator set to Q.

        Due to `dataselect` running as the first step, always sorts the records. `msrepack`
        only runs when `repack` is set to True.

        Qualities:
        D - The state of quality control of the data is indeterminate
        R - Raw Waveform Data with no Quality Control
        Q - Quality Controlled Data, some processes have been applied to the data.
        M - Data center modified, time-series values have not been changed.

        Parameters
        ----------
        `cut_boundaries` : `bool`
            Whether or not to cut the file at the day boundaries ---
            00:00 of the day and of the following day. (default `True`)
        `remove_overlap` : `bool`
            Whether or not to prune the file and remove overlaps at the
            sample level --- equates to the "-Ps" option of dataselect. (default `False`)
        `repack` : `bool`
            Whether or not to repack records using msrepack. (default `False`)
        `record_length` : `int`
            Size of record to repack if `repack` is `True`. (default 4096)

        """
        # Record length within some bounds
        if record_length < 512 or record_length > 65536:
            raise ValueError("Record length is invalid")

        # Confirm record length is power of two
        if record_length & (record_length - 1) != 0:
            raise ValueError("Record length is not is a power of two")

        # Create a phantom SDSFile with a different quality idenfier
        quality_file = SDSFile(self.filename, self.archive_root)
        quality_file.quality = "Q"

        # Create directories for the pruned file (quality Q)
        if not os.path.exists(quality_file.directory):
            os.makedirs(quality_file.directory)

        # Get neighbours
        neighbours = list(map(lambda x: x.filepath, self.neighbours))

        # Define dataselect arguments
        # -Ps prunes to sample level
        # -Q set quality indicator to Q
        # -ts, -te are start & end time of sample respectively (INCLUSIVE)
        # -szs remove records with 0 samples (may result in empty prune	d files)
        # -o - write to stdout
        dataselect_args = ["-Q", "Q", "-szs"]

        # Cut file at the day boundaries if requested
        if cut_boundaries:
            dataselect_args.extend(["-ts", self.sample_start, "-te", self.sample_end])

        # Check if overlap needs to be removed
        if remove_overlap:
            dataselect_args.append("-Ps")
        else:
            dataselect_args.append("-Pe")

        # If the file is going to be repacked, write to stdout
        if repack:
            dataselect_args.extend(["-o", "-"])
        else:
            dataselect_args.extend(["-o", quality_file.filepath])

        # Create a dataselect process
        dataselect = subprocess.Popen(["dataselect"] + dataselect_args + neighbours,
                                      stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)

        # Open a msrepack process and connect stdin to dataselect stdout
        # -R repack record size to record_length
        # -o output file for pruned data
        # - read from STDIN
        if repack:
            msrepack = subprocess.Popen([
               "msrepack",
               "-R", str(record_length),
               "-o", quality_file.filepath,
               "-"
            ], stdin=dataselect.stdout, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Not sure why we need this
        dataselect.stdout.close()

        # Wait for child processes to terminate
        if dataselect.wait() != 0 or (repack and msrepack.wait() != 0):
            if repack:
                raise Exception(("Unable to prune file "
                                 "(dataselect returned %s, msrepack returned %s)")
                                % (str(dataselect.returncode), str(msrepack.returncode)))
            else:
                raise Exception("Unable to prune file (dataselect returned %s)"
                                % (str(dataselect.returncode)))

        # Check that quality file has been created
        if os.path.exists(quality_file.filepath):
            self.logger.debug("Created pruned file %s" % quality_file.filename)
        else:
            raise Exception("Pruned file %s has not been created!" % quality_file.filename)

    def _get_adjacent_file(self, direction):
        """Private function that returns adjacent SDSFile based on direction."""

        new_date = self.start + timedelta(days=direction)

        # The year and day may change
        new_year = new_date.strftime("%Y")
        new_day = new_date.strftime("%j")

        new_filename = ".".join([
            self.net,
            self.sta,
            self.loc,
            self.cha,
            self.quality,
            new_year,
            new_day
        ])

        return SDSFile(new_filename, self.archive_root)

    def __str__(self):
        return "%s (%s)" % (self.filename, self.modified)
