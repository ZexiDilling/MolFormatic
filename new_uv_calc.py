import numpy as np
import pandas as pd
import scipy


def _sample_controller(sample, wavelength, slope_threshold, uv_threshold, data, rt_solvent_peak, blank_run):
    """
    Processes the UV data for a single sample and performs peak integration if necessary.

    :param sample: The name of the sample to process.
    :type sample: str
    :param wavelength: The wavelength to use for processing. If "all", the UV data will be summed across all wavelengths
    :type wavelength: int or str
    :param slope_threshold: The slope threshold used to identify peaks in the UV data.
    :type slope_threshold: float
    :param uv_threshold: The UV intensity threshold used to identify peaks in the UV data.
    :type uv_threshold: float
    :param data: A dictionary containing the UV data for all samples.
    :type data: dict
    :param rt_solvent_peak: The retention time of the solvent peak.
    :type rt_solvent_peak: float
    :param blank_run: A numpy array for the UV data for a blank sample, if there have been a blank sample in the data,
        else it is None
    :type blank_run: numpy.ndarray or None
    :return: A tuple containing the UV tensor, UV retention times, selected wavelength index,
        and number of points used for integration.
    :rtype: tuple
    """
    uv_tensor = np.array(data[sample]["uv"])
    points = 5

    # Define UV retention times in seconds and wavelengths
    uv_rt = data[sample]["uv"].index * 60
    if wavelength == "all":
        uv_tensor = np.sum(uv_tensor, axis=1)
        wave_indx = None

    else:
        wavelength = int(wavelength)
        uv_wavelengths = data[sample]["uv"].columns
        # Find index of selected wavelength
        wave_number = uv_wavelengths[(np.fabs(uv_wavelengths - wavelength)).argmin(axis=0)]
        wave_indx = [i for i, x in enumerate(uv_wavelengths == wave_number) if x][0]

    if blank_run:
        uv_tensor = _subtract_blank(uv_tensor, blank_run)
    slope, start_slope, end_slope = _calc_slopes(uv_tensor, uv_rt, points, wave_indx, uv_threshold)
    rt_start_peak_new, rt_end_peak_new = _find_peaks(slope, slope_threshold, start_slope, end_slope, rt_solvent_peak)
    integrals = _calc_integrals(rt_start_peak_new, rt_end_peak_new, uv_rt, wave_indx, uv_tensor)
    df_integrate_overview = _saving_the_data(rt_start_peak_new, rt_end_peak_new, integrals)

    return df_integrate_overview


def _calc_slopes(uv_tensor, uv_rt, points, wave_indx=None, uv_threshold=0):
    """
    Calculates the slope of the regression line for a set of points in a UV-Spectra tensor.

    :param uv_tensor: numpy array of shape (N, M) representing the UV spectra tensor.
    :type uv_tensor: numpy.ndarray
    :param uv_rt: numpy array of shape (N,) representing the retention times for the UV spectra.
    :type uv_rt: numpy.ndarray
    :param points: The number of points to be used for calculating the slope of the regression line.
    :type points: int
    :param wave_indx: The index of the specific wavelength to be used for the calculations. If None, all wavelengths will be used.
    :type wave_indx: int, Optional
    :param uv_threshold: The threshold below which the UV intensity is not considered.
    :type uv_threshold: float, Optional
    :return: A tuple of three lists representing the calculated slopes, start retention times and end retention times respectively.
    :rtype: Tuple[List[float], List[float], List[float]]
    """
    slope = []
    start_slope = []
    end_slope = []
    for rt in range(len(uv_tensor) - points + 1):
        # Create an array of x values centered at 0 that extend to points - 1
        x = np.arange(points)
        # Select the appropriate y values from the uv_tensor array
        if wave_indx:
            y = uv_tensor[rt:(rt + points), wave_indx]
        else:
            y = uv_tensor[rt:(rt + points)]
        # Skip calculations if the first value of y is below the uv_threshold
        if y[0] > uv_threshold:
            # Calculate the slope of the regression line
            temp_slope, _, _, _, _ = scipy.stats.linregress(x, y)
            slope.append(temp_slope)
            # Calculate the index value for the middle of the window
            index_value = int(round(points/2))
            # Append the start and end retention times for the window to their respective lists
            start_slope.append(uv_rt[rt])
            end_slope.append(uv_rt[rt + points - 1])
    return slope, start_slope, end_slope


def _find_peaks(slope, slope_threshold, start_slope, end_slope, rt_solvent_peak):
    """
    Find the start and end times based on the slope and remove the solvent peak.

    :param slope: The slope values.
    :type slope: numpy array
    :param slope_threshold: The threshold used to determine a slope change.
    :type slope_threshold: float
    :param start_slope: The start time of each slope calculation.
    :type start_slope: numpy array
    :param end_slope: The end time of each slope calculation.
    :type end_slope: numpy array
    :param rt_solvent_peak: The RT for the solvent peak.
    :type rt_solvent_peak: float
    :return: The start and end RTs of the peaks, without the solvent peak.
    :rtype: tuple of two numpy arrays
    """

    # Find the start and end times based on the slope
    rt_start_peak = []
    rt_end_peak = []

    x = -1  # x makes sure the RT's are order correctly
    for s in np.arange(1, len(slope)-1):
        if slope[s] > slope_threshold > slope[s-1] and x < 0:
            rt_start_peak.append(start_slope[s])
            x = 1
        elif slope[s] < -slope_threshold < slope[s+1] and x > 0:
            rt_end_peak.append(end_slope[s])
            x = -1

    # Check if only rt_start_peak was found
    if len(rt_start_peak) > len(rt_end_peak):
        rt_start_peak.pop()

    # Remove the solvent peak
    mask = np.array(rt_start_peak)/60 > rt_solvent_peak
    rt_start_peak_new = np.array(rt_start_peak)[mask]
    rt_end_peak_new = np.array(rt_end_peak)[mask]

    return rt_start_peak_new, rt_end_peak_new


def _calc_integrals(rt_start_peak_new, rt_end_peak_new, uv_rt, wave_indx, uv_tensor):
    """Calculate the integrals of intervals using the trapezoid rule.

    :param rt_start_peak_new: List of start times for each interval.
    :type rt_start_peak_new: list
    :param rt_end_peak_new: List of end times for each interval.
    :type rt_end_peak_new: list
    :param uv_rt: Array of time values.
    :type uv_rt: numpy array
    :param wave_indx: Index of the wave in `uv_tensor` to calculate the integrals for. If `None`, all values are used.
    :type wave_indx: int or None
    :param uv_tensor: Array of UV values.
    :type uv_tensor: numpy array
    :return: List of calculated integrals rounded to the nearest integer.
    :rtype: list
    """
    integrals = []
    for s in np.arange(0, len(rt_start_peak_new)):

        # Find index of start and end times
        indx_start = np.where(uv_rt == rt_start_peak_new[s])[0][0]
        indx_end = np.where(uv_rt == rt_end_peak_new[s])[0][0]

        # Calculate n, a, b and h
        n = len(uv_rt[indx_start:indx_end + 1])
        a = uv_rt[indx_start]
        b = uv_rt[indx_end]
        h = (b - a) / n

        # Define x and y values based on the peak time interval
        x_values = np.linspace(a, b, n)
        if wave_indx:
            y_values = uv_tensor[indx_start:indx_end + 1, wave_indx]
        else:
            y_values = uv_tensor[indx_start:indx_end + 1]

        # Calculate the integrals using the trapezoid integration rule
        integral = h * (np.sum(y_values) - 0.5 * (y_values[0] + y_values[-1]))
        integrals.append(round(integral))

    return integrals


def _saving_the_data(rt_start_peak_new, rt_end_peak_new, integrals):
    """
    This function takes in three arrays as inputs: `rt_start_peak_new`, `rt_end_peak_new`, and `integrals`. It converts
    `rt_start_peak_new` and `rt_end_peak_new` from samples to minutes and calculates the purity of each peak. Finally,
    it saves the information in a pandas DataFrame and returns the DataFrame.

    :param rt_start_peak_new: An array of start times of peaks in samples.
    :type rt_start_peak_new: numpy.ndarray
    :param rt_end_peak_new: An array of end times of peaks in samples.
    :type rt_end_peak_new: numpy.ndarray
    :param integrals: An array of integrals for each peak.
    :type integrals: numpy.ndarray
    :return: A pandas DataFrame containing the peak information.
    :rtype: pandas.DataFrame
    """
    # Create a list of peak names
    peak_list = [f"Peak {rt + 1}" for rt in np.arange(0, len(rt_start_peak_new))]

    # Convert start and end times of peaks from samples to minutes
    rt_start_peak_new = [round(x / 60, 3) for x in rt_start_peak_new]
    rt_end_peak_new = [round(x / 60, 3) for x in rt_end_peak_new]

    # Calculate the purity of each peak
    purity = [round((x / sum(integrals)) * 100, 3) for x in integrals]

    # Create the DataFrame with the peak information
    df_integrate_overview = pd.DataFrame({
        "Peak list": peak_list,
        "Integrals": integrals,
        "Peak start time": rt_start_peak_new,
        "Peak end time": rt_end_peak_new,
        "purity": purity
    })

    # Check if the DataFrame is empty and replace it with a DataFrame containing "No peaks detected" if it is
    if df_integrate_overview.empty:
        df_integrate_overview = pd.DataFrame({
            "Peak list": ["No peaks detected"],
            "Integrals": ["No peaks detected"],
            "Peak start time": ["No peaks detected"],
            "Peak end time": ["No peaks detected"],
            "purity": ["No peaks detected"]
        })

    # Return the DataFrame
    return df_integrate_overview


def _subtract_blank(sample, blank):
    """
    Subtracts the blank HPLC run from a different HPLC run.

    :param sample: The sample HPLC run
    :type sample: numpy.ndarray
    :param blank: The blank HPLC run
    :type blank: numpy.ndarray
    :return: The subtracted HPLC run
    :rtype: numpy.ndarray
    """
    return sample - blank


def uv_integrals_controller(data, slope_threshold, uv_threshold, solvent_peak, sample_data, wavelength_data,
                           sample):
    """
    Calculate integrals for all samples in a plate and combine them to a list of pandas dataframes.
    One dataframe for each sample.
    A dataframe contains the following information, peak list, integrals, peak start time and peak and time.

    :param data: All data for all the samples
    :type data: dict
    :param slope_threshold: The minimum increase in uv intensity before a slope is consideret a peak
    :type slope_threshold: int
    :param uv_threshold: minimum threshold before data is being looked at. everything below is ignored
    :type uv_threshold: int
    :param solvent_peak: retentions time for the solvent peak, to avoid data from that peak.
    :type solvent_peak: float
    :param sample_data: A dict of data for each sample
    :type sample_data: dict or None
    :param wavelength_data: a str or int for all wavelength that are not in sample_data
    :type wavelength_data: str or int
    :param sample: A sample
    :type sample: str or None
    :return: df_combined is a list of pandas dataframes.
    :rtype: dict
    """
    peak_information = {}

    # If there are no Blank runs, then this will be None. #TODO Check if this is working.
    blank_run = None
    # For now only works for all samples. Needs to make a new dict over samples and batches, to run this over,
    # to make it more flexible. It needs to be added to all methods using the module.

    sample_list = []
    if sample:
        sample_list.append(sample)
    else:
        for samples in data:
            sample_list.append(samples)

    if "blank" in sample_list:
        blank_run = np.array(data["blank"]["uv"])

    for samples in sample_list:
        peak_information[samples] = {}
        try:
            wavelength = sample_data[samples]["wavelength"]
        except TypeError or KeyError:
            wavelength = wavelength_data

        if sample == "blank":
            blank_run = np.array(data[sample]["uv"])

        temp_data_frame = _sample_controller(samples, wavelength, slope_threshold, uv_threshold, data, solvent_peak,
                                             blank_run)

        peak_information[samples] = temp_data_frame

    return peak_information
