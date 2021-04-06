# Introduction

PSD calculation and storage module extracted from ORFEUS Data Center infrastructure. This module only calculates and compresses PSD segments following `Koymans et al., 2021`.

# PSD Calculation

The input is given as an `SDSFile` which needs to be instantiated with a filename and SDS archive root using the SeisComP3 archive structure. Then call `processAndStore` in psd.py with the SDSFile.

    PSD = PSDCollector()
    sdsfile = SDSFile("NL.HGN.02.BHZ.001.2020", "/data/archive/SDS/")
    PSD.processAndStore(sdsfile)

# SDSFile

SDSFile is an internal library used to handle SDS (mSEED) files in the SeisComP3 archive directory format.

# Storage

The currently configured storage uses a MySQL table. This table will need to be set up and connected to in `PSDCollector.__init__()`.

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

# Plotting

The script `plot.py` has functions to plot spectogram, PPSD, PSD from segments stored in the compressed database.

# ResponseCacher

Module to cache the response from an FDSNWS webservice so it does not need to be evaluated for every PSD segment. See the README.md inside.