import numpy as np
import matplotlib.pyplot as plt
import pymysql.cursors

from constants import (INVALID, MINIMUM_PERIOD, NUMBER_OF_DB_BINS,
                       NUMBER_OF_FREQUENCIES, OCTAVE_PERIOD_STEP)


class PSDSettings():

    """
    Class PSDSettings
    Shared class for all settings for PSD plotting 
    """

    def __init__(self):
        self.setupFrequencies()

    def plotNoiseModels(self, mode):
        """
        def plotNoiseModels
        Plots seismic or infrasound noise models
        """

        if mode == "infrasound":
            periods, (low, high) = self.plotNoiseModelsInfrasound()
        elif mode == "seismic":
            periods, (low, high) = self.plotNoiseModelsSeismic()

        plt.plot(periods, low, c="grey", linestyle="--")
        plt.plot(periods, high, c="grey", linestyle="--")

    def setupFrequencies(self):

        # This is going to give us nice frequencies at 1, 2, 4Hz
        steps = np.arange(NUMBER_OF_FREQUENCIES) * OCTAVE_PERIOD_STEP
        numbers = MINIMUM_PERIOD * 2 ** steps

        self.frequencies = np.reciprocal(numbers)

    def plotNoiseModelsSeismic(self):
        """
        def plotNoiseModelsSeismic
        Plots Peterson, 1993 NHNM, NLNM
        """

        data = np.load("data/noise_models.npz")
        periods = data["model_periods"]
        low = data["low_noise"]
        high = data["high_noise"]

        return periods, (low, high)

    def plotNoiseModelsInfrasound(self):
        """
        def plotNoiseModelsInfrasound
        Infrasound noise models given by Jelle
        """

        freqs = np.array([
            0.01000, 0.02000, 0.03000, 0.04000, 0.05000, 0.06000, 0.07000,
            0.08000, 0.09000, 0.10000, 0.20000, 0.30000, 0.40000, 0.50000,
            0.60000, 0.70000, 0.80000, 0.90000, 1.00000, 2.00000, 3.00000,
            4.00000, 5.00000, 6.00000, 7.00000, 8.00000, 9.00000
        ])

        high = np.array([
            6.30957, 7.95271, 6.89542, 6.32456, 5.61675, 4.35588, 3.73722,
            3.17355, 3.17776, 3.16228, 2.51487, 2.18052, 1.58866, 1.12069,
            0.86911, 0.74567, 0.56435, 0.44887, 0.37584, 0.33536, 0.23097,
            0.21185, 0.22361, 0.23125, 0.23580, 0.17846, 0.16870
        ])

        low = np.array([
            0.00944, 0.00596, 0.00488, 0.00356, 0.00316, 0.00259, 0.00250,
            0.00225, 0.00212, 0.00200, 0.00399, 0.00173, 0.00126, 0.00126,
            0.00123, 0.00112, 0.00089, 0.00075, 0.00063, 0.00017, 0.00010,
            0.00007, 0.00006, 0.00006, 0.00005, 0.00004, 0.00003
        ])

        # TODO MAKE SURE
        low = 20 * np.log10(low)
        high = 20 * np.log10(high)

        return np.reciprocal(freqs), (low, high)


class PSDTimeSeries(PSDSettings):

    """
    Class PSDTimeSeries
    Plots time series of a single frequency (time, ampl)
    """

    def __init__(self):

        PSDSettings.__init__(self)

        # Empty container for the spectrogram
        self.spectrogram = list()

    def add(self, shift, offset, binary):
        thing = -np.frombuffer(binary, dtype=np.uint8) + shift
        prepend = np.empty(offset) * np.nan
        append = NUMBER_OF_FREQUENCIES - len(thing) - len(prepend)
        append = np.empty(append) * np.nan
        self.spectrogram.append(np.concatenate([
            prepend, thing, append
        ]))

    def plot(self, frequency):

        spectrogram = np.array(self.spectrogram)

        # Get the frequency closest to what is requested
        closest = np.abs(self.frequencies - frequency).argmin()

        plt.title("Frequency " + str(self.frequencies[closest]))

        # Extract the particular row from the matrix
        plt.plot(spectrogram[:, closest])
        plt.show()


class PSDSpectrogram(PSDSettings):

    """
    Class PSDSpectrogram
    Plots a PSD Spectrogram (time, freq, ampl)
    """

    def __init__(self):

        PSDSettings.__init__(self)

        # Empty container for the spectrogram
        self.spectrogram = list()

    def add(self, shift, offset, binary):
        thing = -np.frombuffer(binary, dtype=np.uint8) + shift
        prepend = np.empty(offset) * np.nan
        append = NUMBER_OF_FREQUENCIES - len(thing) - len(prepend)
        append = np.empty(append) * np.nan
        self.spectrogram.append(np.concatenate([
            prepend, thing, append
        ]))

    def plot(self):

        self.spectrogram = np.array(self.spectrogram)
        freqsEdges = self.frequencies * (2 ** 0.125) ** 0.5

        # Matplotlib options
        plt.title("Spectrogram")
        plt.xlabel("Segment Number")
        plt.ylabel("Frequency (Hz)")

        plt.semilogy()

        # Show the spectrogram
        plt.pcolormesh(np.arange(len(self.spectrogram)),
                       freqsEdges, self.spectrogram.T)

        # Add the spectrogram colorbar
        colorbar = plt.colorbar()
        colorbar.set_label("Amplitude [(m/s^2)^2/Hz][dB]")

        plt.show()


class PSDHistogram(PSDSettings):

    """
    Class PSDHistogram
    Container for plotting a PSD histogram based on PSD segments
    """

    def __init__(self):

        PSDSettings.__init__(self)

        # Create an empty 2D array
        self.histogram = np.zeros((
            NUMBER_OF_DB_BINS,
            NUMBER_OF_FREQUENCIES
        ))

        # Track number of segments used
        self.nSegments = 0

    def add(self, shift, offset, binary):
        """
        PSDHistogram.add
        Add a segment to the histogram
        """

        if shift != 0:
            raise ValueError(
                "TODO: The histogram cannot accomodate shifted data.")

        self.nSegments += 1

        # Go over the spectrum and fill the histogram
        for i, value in enumerate(binary):

            # Do not show invalid data
            if value != INVALID:
                self.histogram[255 - value][offset + i] += 1

    def plot(self):
        """
        PSDHistogram.plot
        Normalizes and plots the histogram
        """

        # Could be hoisted
        freqsEdges = self.frequencies * (2 ** 0.125) ** 0.5

        # Noise models
        # if trace.stats.channel.endswith("DF"):
        # self.plotNoiseModels("infrasound")
        # else:
        self.plotNoiseModels("seismic")

        # Normalize to percentage and clip between 0% - 30%
        hist = 100 * np.divide(self.histogram, self.nSegments)
        hist = np.clip(hist, 0, 30)

        # Set up the plot
        plt.ylim(-200, -50)
        plt.xlim(1E-2, 200)
        plt.xlabel("Period (s)")
        plt.ylabel("Amplitude [(m/s^2)^2/Hz][dB]")

        # Prepare the axes
        plt.semilogx()
        plt.pcolormesh(np.reciprocal(freqsEdges),
                       np.arange(-255, 0) + 0.5, hist)
        plt.colorbar()
        plt.grid(True, which="both")
        plt.show()


if __name__ == "__main__":

    connection = pymysql.connect(
        user="root",
        password="password",
        host="localhost",
        db="test"
    )

    mode = "histogram"
    #mode = "spectrogram"
    #mode = "timeseries"

    # Create a PSDTimeSeries, PSDSpectrogram, PSDHistogram (PPSD)
    if mode == "spectrogram":
        figure = PSDSpectrogram()
    elif mode == "histogram":
        figure = PSDHistogram()
    elif mode == "timeseries":
        figure = PSDTimeSeries()

    try:

        with connection.cursor() as cursor:

            # May need to be ordered when spectrogram / timeseries
            # Remember to index the table
            if mode == "spectrogram" or mode == "timeseries":
                cursor.execute(
                    "SELECT shift, offset, spectrum FROM PSD ORDER BY start")
            else:
                cursor.execute("SELECT shift, offset, spectrum FROM PSD")

            for (shift, offset, spectrum) in cursor.fetchall():
                figure.add(shift, offset, spectrum)

    finally:
        connection.close()

    # Plot timeseries for a given frequency (4)
    if mode == "timeseries":
        figure.plot(4)
    else:
        figure.plot()
