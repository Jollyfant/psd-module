"""
Convert the download instrument responses to a hashfiles
and include a hashmap for quick lookup
"""

import os
import json
import shutil
import datetime
import hashlib
from dateutil.parser import parse
from config import config


def shasum256(filename):
    """
    def shasum256
    Calculates the sha256 hash of a file
    """

    sha256_hash = hashlib.sha256()

    with open(filename, "rb") as f:
        for block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(block)
        return sha256_hash.hexdigest()


if __name__ == "__main__":

    # Dict for lookup!
    lookupDict = dict()

    # Wrong order?
    if not os.path.exists(config["RESPFILE_DIR"]):
        raise ValueError("No response folder found. Run cache.py first.")

    # Create necessary directories
    if not os.path.exists(config["HASHFILE_DIR"]):
        os.mkdir(config["HASHFILE_DIR"])

    # Go over all response files
    for filename in os.listdir(config["RESPFILE_DIR"]):

        print("Converting %s." % filename)

        filepath = os.path.join(config["RESPFILE_DIR"], filename)

        # Extract the seed identifier
        network, station, location, channel, starttime, endtime, * \
            _ = filename.split(".")
        identifier = "%s.%s.%s.%s" % (network, station, location, channel)

        # Add the seed identifier to the lookup dict
        if not identifier in lookupDict:
            lookupDict[identifier] = list()

        # Calculate the hash of the file
        shasum = shasum256(filepath)

        # Add an entry a reference pointer to disk
        lookupDict[identifier].append({
            "starttime": parse(starttime),
            "endtime": parse(endtime),
            "file": shasum
        })

        # Lets move this file to a new path if it does not exist
        newpath = os.path.join(config["HASHFILE_DIR"], shasum)

        if not os.path.exists(newpath):
            shutil.copyfile(filepath, newpath)

    # Dump the pointer map
    with open(config["HASHMAP_FILE"], "w") as outfile:
        json.dump(lookupDict, outfile, default=str)
