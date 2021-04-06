# Response Cacher

Queries FDSNWS for inventories and evaluates frequencies, amplitudes for a given segment length. The segment length must be consistent to the NFFT used in the calculation of the PSD.

## Running

1. Run `cache.py` to get StationXML response files from FDSNWS, evalauted to frequencies, amplitudes and saved as numpy compressed files
2. Run `convert.py` to convert the previous files to hashfiles to eliminate duplicates and create a metadata reference map for all channels to files.
3. Use `response-map.json` to discover channels in a certain timewindow and find the respectively hashfile on disk.
