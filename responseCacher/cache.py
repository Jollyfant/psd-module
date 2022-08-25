"""

Script to cache instrument response to file for calculation in PSDs
The instrument response is always evaluated at the same frequencies and can thus be re-used

Should be run once to download all instruments from the FDSNWS webservices and evaluate the responses
These are written to files

Author: Mathijs Koymans, 2020

"""

import numpy as np
import os
import json
import requests

from obspy import read_inventory
from obspy.core import UTCDateTime
from config import config

if __name__ == "__main__":

    # Read the existing hashmap
    if os.path.exists(config["HASHMAP_FILE"]):
        with open(config["HASHMAP_FILE"], "r") as infile:
            hashMap = json.load(infile)
    else:
        hashMap = {}

    # Create response file folder
    if not os.path.exists(config["RESPFILE_DIR"]):
        os.mkdir(config["RESPFILE_DIR"])

    # Channels for which to cache evaluated instrument response
    r = requests.get(
        "%s?level=channel&format=text&network=NL&station=HGN" % config["API_URL"])

    for line in r.text.split("\n")[1:-1]:

        # Extract parameters from FDSNWS
        (network, station, location, channel,
         latitude, longitude, elevation, depth,
         azimuth, dip, description, scale, freq,
         units, fs, starttime, endtime) = line.split("|")

        # Skip
        if channel not in config["PSD_CHANNELS"]:
            continue

        # Already exists within the hashmap: remove when forced
        if "%s.%s.%s.%s" % (network, station, location, channel) in hashMap:
            continue

        # Sampling rate quickly to float
        fs = float(fs)

        # Set end time far in to the future
        if endtime == "":
            endtime = "2030-01-01T00:00:00"

        # Evaluate the response to VEL for infrasound channels
        # Seismic stations should go to ACC to be consistent with the NLNM, NHNM
        if channel.endswith("DF"):
            output = "VEL"
        else:
            output = "ACC"

        # Segment length in seconds: more for gravity to see longer periods
        # XXX IMPORTANT! These must be the same length as the PSD calculation
        # For infrasound and seismic we use 1h segments (ObsPy, McNamara default)
        if channel == "LGZ":
            segment_length = 24 * 3600
        else:
            segment_length = 3600

        t0 = starttime
        t1 = endtime
        starttime = (UTCDateTime(starttime)+1.0).isoformat()
        endtime = (UTCDateTime(endtime)-1.0).isoformat()

        # Prepare to read inventory using ObsPy
        requestOptions = (
            config["API_URL"],
            network, station,
            location, channel,
            starttime, endtime
        )

        # May return 204
        try:
            inventory = read_inventory(
                "%s?network=%s&station=%s&location=%s&channel=%s&starttime=%s&endtime=%s&level=response" % requestOptions)
        except:
            continue

        # Verify that only one channel is returned
        contents = inventory.get_contents()
        if len(contents["networks"]) != 1 or len(contents["stations"]) != 1 or len(contents["channels"]) != 1:
            continue

        # Get the first response object
        response = inventory[0][0][0].response

        # Call evalresp to evaluate the response
        # NFFT must be same as in Welch's method!
        # We use 13 segments with 75% overlap
        resp, freqs = response.get_evalresp_response(
            t_samp=np.reciprocal(fs),
            nfft=int(fs * segment_length / 4.0),
            output=output
        )

        # Drop phase imaginary information
        resp = np.abs(resp)

        # Save the file as compressed NPZ
        filename = "%s.%s.%s.%s.%s.%s.response" % (
            network, station, location, channel, starttime, endtime)
        filepath = os.path.join(config["RESPFILE_DIR"], filename)

        print("Saving %s." % filename)

        # Write to compressed file
        np.savez_compressed(
            filepath,
            freqs=freqs,
            resp=resp
        )
