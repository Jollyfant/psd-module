import numpy as np
import os
import json
from matplotlib import mlab
import matplotlib.pyplot as plt
from dateutil.parser import parse
from scipy.signal import detrend, welch, tukey
import logging


logger = logging.getLogger("RuleManager")


def getInstrumentResponse(stats):
    """
    def getInstrumentResponse
    Loads the cached instrument response
    """

    # Settings
    HASHMAP_FILENAME = "responseCacher/response-map.json"
    HASHMAP_FOLDER = "responseCacher/hashfiles"

    # Create the SEED identifier from the trace.stats
    # The hashmap uses this as a channel identifier
    identifier = "%s.%s.%s.%s" % (
        stats.network,
        stats.station,
        stats.location,
        stats.channel
    )

    # Load the inventory lookup dictionary
    with open(HASHMAP_FILENAME) as json_file:
        lookupDict = json.load(json_file)

    # Go over all ranges for the SEED identified and return the first one that matches
    for ranges in lookupDict[identifier]:
        if (parse(ranges["starttime"]) < stats.starttime) and (stats.endtime < parse(ranges["endtime"])):
            cache = np.load(os.path.join(HASHMAP_FOLDER, ranges["file"]))
            return cache["resp"], cache["freqs"]

    raise ValueError(
        "Could not the find instrument response in cached folder.")


def cutSpectrum(spectrum):
    """
    def cutSpectrum
    Cuts off unnecessary np.nan values from start, end
    """

    frontOffset = np.argmin(np.isnan(spectrum))
    backOffset = len(spectrum) - np.argmin(np.flip(np.isnan(spectrum), 0))

    # Cut the spectrum
    return frontOffset, spectrum[frontOffset:backOffset]


def compressSpectrum(spectrum):
    """
    def compressSpectrum
    Compresses the spectrumt to uint8
    """

    # Invalid values are 0
    INVALID = 0

    # Cut the spectrum
    frontOffset, spectrum = cutSpectrum(spectrum)

    # Get the dB shift to bring the minimum value to 0
    # Most effective way to shift the spectrum
    shift = np.nanmin(spectrum)

    # Shift the full spectrum in to the positive range
    spectrum -= shift

    # Set all np.nan values and values outside of range to invalid
    # Some nan values may actually remain in between other values
    spectrum[np.isnan(spectrum)] = INVALID
    spectrum[(spectrum < 0) | (spectrum > 255)] = INVALID

    # Cast to binary uint8 array for MySQL storage and add front offset
    # and shift for reconstruction
    return int(frontOffset), int(shift), spectrum.astype(np.uint8).tobytes()

def prevPower2(nfft):
    """
    def prevPower2
    Returns the previous power of 2 of the input number
    """

    return int(2 ** (np.floor(np.log2(nfft))))


def smoothSpectrum(freqs, Pxx, frequencies):
    """
    def smooth
    Smoothes spectrum over a fraction of octaves
    """

    OCTAVE_PERIOD_SMOOTHING = 1
    ROUND_POWER = True

    # Smoothing factor
    factor = 2 ** OCTAVE_PERIOD_SMOOTHING

    smoothed = list()

    # Go over the requested frequency bins
    for f in frequencies:

        # Get the lower and upper bounds of the range that should be smoothed over
        lowerEdge = f / factor ** 0.5
        upperEdge = lowerEdge * factor

        # Extract all powers within the bounds and take the mean value
        segment = Pxx[(freqs >= lowerEdge) & (freqs <= upperEdge)]

        # No elements withinin the range: set nan otherwise mean
        if segment.size == 0:
            smoothed.append(np.nan)
        else:
            smoothed.append(np.mean(segment))

    smoothed = np.array(smoothed)

    # Determines whether to round to full dB bins
    if ROUND_POWER:
        return np.round(smoothed)
    else:
        return smoothed


def _welch(signal, fs, resp, nfft=1024, overlap=0.5):
    """
    def _welch
    Calculates power spectral densities without the help of SciPy, Mlab
    """

    # Create a tukey (cosine taper) window of size nfft
    window = tukey(nfft, alpha=0.2)

    # Calculate number of points from the overlap percentage
    overlap = int(nfft * overlap)

    # Normalize the power for the window gain (squared)
    # and frequency resolution (Fs / NFFT)
    scaling = np.reciprocal(fs * np.sum(window * window))

    # Unneeded since evalresp already calculates the frequencies!
    #freqs = np.fft.rfftfreq(nfft, d=np.reciprocal(fs))

    # Container for all the spectra
    Pxx = list()

    # Welch method: slide using a segment of length NFFT
    while(signal.size >= nfft):

        # Cut the segment to the segment length
        segment = signal[:nfft]

        # Detrend the trace by removing a linear relationship
        segment = detrend(segment, type="linear")

        # Apply the window and do the FFT
        # Deconvolve the instrument response by divison in the frequency domain
        with np.errstate(divide="ignore", invalid="ignore"):
            yf = np.fft.rfft(window * segment) / resp

        # Square FFT amplitudes to get the power
        # Half the spectrum is calculating using rfft; so double the power
        yf = 2 * np.abs(yf) ** 2

        # Square the amplitude to get the power
        Pxx.append(yf.real)

        # Cut to the next segment
        signal = signal[nfft - overlap:]

    # Normalize the power for the sampling rate, half power, etc.
    return scaling * np.array(Pxx)


def psdWelchMlab(signal, fs, resp, nfft=1024, overlap=0.5, reference=1.0):
    """
    def psdWelchMlab
    Calculates power spectral densities using matplotlib.mlab
    """

    overlap = int(nfft * overlap)

    # Use mlab
    Pxx, _ = mlab.psd(
        signal,
        NFFT=nfft,
        Fs=fs,
        detrend=mlab.detrend_linear,
        noverlap=overlap,
        sides="onesided",
        scale_by_freq=True
    )

    # Save convention as power
    Pxxn = Pxx / resp ** 2

    # Convert to decibels
    return 10 * np.log10(Pxxn / reference ** 2)


def tukeyWindow(data):

    data *= tukey(len(data), alpha=0.2)
    return data


def psdWelchScipy(signal, fs, resp, nfft=1024, overlap=0.5, reference=1.0):
    """
    def psdWelchScipy
    Calculates power spectral densities using scipy
    """

    overlap = int(nfft * overlap)

    _, Pxx = welch(
        signal,
        fs=fs,
        nperseg=nfft,
        noverlap=overlap,
        nfft=nfft,
        detrend="linear",
        window=tukeyWindow,
        return_onesided=True,
        scaling="density"
    )

    # Instrument deconvolution (as power)
    Pxxn = Pxx / resp ** 2

    # Convert to decibels
    return 10 * np.log10(Pxxn / reference ** 2)


def spectrogram(*args, reference=1.0, **kwargs):
    """
    def spectrogram
    Calculates spectrogram using Welch's method
    """

    Pxx = _welch(*args, **kwargs)

    return 10 * np.log10(Pxx / reference ** 2)


def psdWelch(*args, reference=1.0, **kwargs):
    """
    def psdWelch
    Calculates power spectral densities
    """

    # PSD can be averaged from the spectrogram
    Pxx = _welch(*args, **kwargs)

    # Get the mean from all the windows
    mean = np.mean(Pxx, axis=0)

    return 10 * np.log10(mean / reference ** 2)
