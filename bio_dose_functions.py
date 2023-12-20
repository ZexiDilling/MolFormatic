import pandas as pd
from scipy.optimize import leastsq, brentq
import numpy as np


def setup_t20_colour_list():
    """ Setup a list of colours for the figures.

    Returns
    -------
    t20 : list
        List of RGB colour tuples, normalised between 0 and 1 as used by python colours.
    """
    """ Setup colours for the figures
    return : t20, a list of colour tuples
    """
    colour_lists_dict = {}
    # define colour lists
    colour_lists = {
        'tableau20': [
            (31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
            (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
            (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
            (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
            (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)
        ],
        'tableau20blind': [
            (0, 107, 164), (255, 128, 14), (171, 171, 171), (89, 89, 89),
            (95, 158, 209), (200, 82, 0), (137, 137, 137), (163, 200, 236),
            (255, 188, 121), (207, 207, 207)
        ]
    }
    # normalise the colours for the colour lists
    for rgb_list in colour_lists:
        colour_array = np.array(colour_lists[rgb_list]) / 255.
        colour_array_tup = tuple(map(tuple, colour_array))
        colour_lists[rgb_list] = colour_array_tup
        # add normalised colours to output dictionary
        colour_lists_dict[rgb_list] = colour_lists[rgb_list]
    # extract the tableau20 lists, join together
    t20 = list(colour_lists_dict["tableau20"] + colour_lists_dict["tableau20blind"])
    # extend the list in case someone need a long list of colours (after 30 colours, they will be redundant!)
    t20 = t20 + t20 + t20 + t20 + t20
    return t20


def residuals(constants, function, x, y):
    """
    Function used to optimise the fit of the curve to the data.
    It calculates the distance between y-value from real data and y-value from the function (sigmoid/sine/etc).
    """
    return y - function(constants, x)


def hill_eq(hill_constants, x):
    """ Four parameter sigmoidal Hill equation.

    y = upper + (lower-upper)/(1+(x/EC50)**-hillslope)

    Parameters
    ----------
    hill_constants : tuple
        Tuple of the four parameters : upper, lower, EC50 and hillslope
    x : float
        x-value for use in the equation

    Constants
    ---------
    upper : float
        Renamed from "bottom" in some sources.
        Will approach the minimum response (perhaps zero) in the dose-response S-curve.
        Will approach the maximum response in an inverse LD50 curve.
        Is not currently used as a parameter to judge whether curve is sigmoidal.
    lower : float
        Renamed from "top" in some sources.
        Will approach the maximum response in the dose-response S-curve.
        Will approach zero in an inverse LD50 curve.
        Is not currently used as a parameter to judge whether curve is sigmoidal.
    EC50 : float
        Does not always accurately reflect the EC50 from the data.
        Is not currently used as a parameter to judge whether curve is sigmoidal.
        Is not currently used as the calculated EC50 value.
    hillslope : float
        A high hillslope generally indicates a steep curve, an accurate EC50 calculation, and a strong sigmoidal shape.
        Note that the hillslope can be strongly negative.
        Hillslope values approaching zero (e.g. -1 > hill_slope > +1) generally indicate an exponential, rather than
        a sigmoidal curve. For this reason, hillslope is currently used as a parameter to judge whether the curve is
        sigmoidal, using the parameters chosen in the settings_excel_file and implemented in the "judge_fit" program.

    Notes
    -----
    The hill equation is used for the EC50 calculation, but there are several variants of the formula.

    Dose-Response Equations
    variant 1: http://www.graphpad.com/guides/prism/6/curve-fitting/index.htm?reg_classic_dr_variable.htm
    variant 2: http://en.wikipedia.org/wiki/EC50
    variant 3: http://pydoc.net/Python/cgptoolbox/0.1.2/cgp.sigmoidmodels.doseresponse/

    Sigmoid equations:
    variant 4: y = c * 1.0 / (1.0 + ((k/x)**g))
    variant 5: y = c / (1 + np.exp(-k*(x-x0))) + y0
        (http://stackoverflow.com/questions/7588371/scipy-leastsq-goodness-of-fit-estimator)

    This equation uses variant 2. See Wikipedia for a more detailed explanation.

    In theory, any formula that accurately fits a sigmoidal curve is appropriate for the purpose. The Hill equation
    is preferred, because under some circumstances the EC50 value is one of the fitted parameters in the hill_constants
    tuple. However this value is extremely unreliable, and regularly gives EC50 values that are out of the range of the
    datapoints. A more accurate EC50 calculation is achieved through root finding using the brent equation.

    """
    upper, lower, EC50, hillslope = hill_constants
    y = upper + (lower - upper) / (1 + (x / EC50) ** -hillslope)
    return y


def hill_eq_brentq(xvalues_for_curve, hill_constants, y_value_curve_center):
    """ Residual function for the four parameter sigmoidal Hill equation.

    For further detail on the Hill equation, see the relevant docstring for hill_eq.

    y = hill_eq(x) - y_value_curve_center

    Parameters
    ----------
    xvalues_for_curve : array
        Numpy array of >250 x-values between the lowest and highest dose.
    hill_constants : tuple
        Tuple of the four parameters : upper, lower, EC50 and hillslope
    y_value_curve_center : array
        Represents the y-value of the curve. The difference between this and the real y-value from experimental
        data (e.g. y_orig_norm) is used to optimise the fit of the curve to the data.

    Returns
    -------
    y - y_value_curve_center : array
        Array that approaches zero, when an optimum fit is found.
    """
    upper, lower, EC50, hillslope = hill_constants

    y = upper + (lower - upper) / (1 + (xvalues_for_curve / EC50) ** -hillslope)

    return y - y_value_curve_center


def normalise_0_1(data_list):
    """ Normalise an array to values between 0 and 1.

    The following linear formula is used.
    norm_array = (orig_array - array_min)/(array_max - array_min)

    The use of this simple linear formula allows the normalised data to be "denormalised" later, so long as
    the min and max values of the original array are known.

    Parameters
    ----------
    arraylike : array
        Numpy array (or other arraylike) dataset of floats or ints to be normalised.

    Returns
    -------
    normalised : array
        Array of floats, containing the normalised datapoints.
    array_min : float
        Minimum value of the original data. Necessary in order to "denormalise" the data later, back to the effective
        original values.
    array_max : float
        Maximum value of the original data. Necessary in order to "denormalise" the data later, back to the effective
        original values.

    Usage
    -----
    normalised_array, min_, max_ = normalise_0_1(original_array)
    # or, if denormalisation is not necessary
    normalised_array = normalise_0_1(original_array)[0]
    # for further usage examples, see the docstring for denormalise_0_1
    """
    arraylike = np.array(data_list)
    array_min = np.min(arraylike)
    array_max = np.max(arraylike)
    normalised = (arraylike - array_min) / (array_max - array_min)
    # convert to float
    normalised = np.array(normalised).astype(float)

    return normalised, array_min, array_max


def denormalise_0_1(value_or_array, array_min, array_max):
    """ Denormalise a value or array to orig values.

    For use after normalisation between 0 and 1 with the normalise_0_1 function.

    The normalisation formula (normalise_0_1):
        norm_array = (orig_array - array_min)/(array_max - array_min)

    The denormalisation formula (denormalise_0_1):
        denormalised_array = norm_array*(array_max - array_min) + array_min

    Parameters
    ----------
    value_or_array : int, float or arraylike
        Int or float to be denormalised.
        Numpy array (or other arraylike) of data (float, int, etc) to be denormalised.

    Returns
    -------
    normalised : float, or numpy array
        Array of floats, containing the normalised datapoints.
    array_min : float
        Minimum value of the original data. Necessary in order to "denormalise" the data later, back to the effective
        original values.
    array_max : float
        Maximum value of the original data. Necessary in order to "denormalise" the data later, back to the effective
        original values.

    Usage
    -----
    from eccpy.tools import normalise_0_1, denormalise_0_1
    import numpy as np
    original_array = np.linspace(10,130,10)
    original_array[2], original_array[4] = 3, 140
    print(original_array)
    # normalise original array
    normalised_array, min_, max_ = normalise_0_1(original_array)
    print(normalised_array)
    # do stuff to normalised array (e.g., multiply by 0.5)
    normalised_array_halved = normalised_array * 0.5
    # denormalise values to match the equivalents in the original array.
    # Note that the min value (3) was normalised to zero, and was therefore not affected by multiplication.
    normalised_array_halved_denorm = denormalise_0_1(normalised_array_halved, min_, max_)
    print(normalised_array_halved_denorm)
    # now calculate average values, and check that they match
    norm_array_mean = np.mean(normalised_array)
    norm_array_mean_denormalised = denormalise_0_1(norm_array_mean, min_, max_)
    orig_array_mean = np.mean(original_array)
    # print the two mean values. They should be equal.
    print(norm_array_mean_denormalised)
    print(orig_array_mean)
    """


    if isinstance(value_or_array, list):
        raise ValueError('this function accepts arraylike data, not a list. '
                         'Please check data or convert list to numpy array')
    elif isinstance(value_or_array, float):
        # print("found a float")
        denormalised = value_or_array * (array_max - array_min) + array_min
    elif isinstance(value_or_array, np.ndarray):
        # print("found an array")
        denormalised = value_or_array * (array_max - array_min) + array_min
    elif isinstance(value_or_array, pd.Series):
        # print("found a series")
        denormalised = value_or_array * (array_max - array_min) + array_min
    else:
        print("Unknown datatype. denormalise_0_1 has been given an input that does not appear to be "
              "an int, float, np.ndarray or pandas Series\n"
              "Attempting to process as if it is arraylike.....")

    return denormalised


def _calc_rsquared(f_vec, reading_normalized):
    """
    obtain the rsquared value for the fit of the curve to the data
    code is from http://stackoverflow.com/questions/7588371/scipy-leastsq-goodness-of-fit-estimator
    :param f_vec:
    :param reading_normalized:
    :return:
    """
    ss_err = np.sum(np.array(f_vec) ** 2)
    ss_tot = np.sum((reading_normalized - reading_normalized.mean()) ** 2)
    rsquared = 1 - (ss_err / ss_tot)

    return rsquared


def _calc_residuals_mean(hill_constants, temp_data):
    dose_normalized = temp_data["dose"]["normalized"]
    reading_normalized = temp_data["reading"]["normalized"]
    reading_min = temp_data["reading"]["min"]
    reading_max = temp_data["reading"]["max"]
    reading_fitted_dose_normalized = hill_eq(hill_constants, dose_normalized)

    residuals_normalized = abs(reading_normalized - reading_fitted_dose_normalized)

    residuals_normalized_mean = residuals_normalized.mean()
    residuals_mean = denormalise_0_1(residuals_normalized_mean, reading_min, reading_max)

    return residuals_mean


def _calc_EC50_brent_eq(sample_name, hill_constants, ec50_normalized):
    """
    The brent equation(in python the brentq function) is used here to find the position on the x-axis at a
    particular y-value. This tends to be more robust than the EC50 parameter derived from the Hill equation.

    In some cases, the brentq function will result in ValueError: f(a) and f(b) must have different signs.
    This is due to a value outside of the initial range (between 0 and 1), and usually indicates a very poor fit
    http://stackoverflow.com/questions/22277982/how-to-find-50-point-after-curve-fitting-using-numpy/22288826#22288826
    :param sample_name:
    :param hill_constants:
    :param ec50_normalized:
    :return:
    """

    # try to determine the EC50 between 0.0 and 1.0, which works for most data

    try:
        EC50_norm_bq = brentq(hill_eq_brentq, 0.0, 1.0, args=(hill_constants, ec50_normalized))
        EC50_calculable = True
    except ValueError:
        # widen the scope of the EC50 to outside of the actual range of datapoints.
        # This generally indicates poor data quality, however there may be exceptions. The need for a wider range
        # is therefore printed in the console, however this wider scope is not a feature in the judge_fit algorithm,
        # and therefore is not a factor used to determine EC50 data quality.
        try:
            print("ValueError encountered in brentq. Attempting wider scope for EC50. "
                  "This is unusual, and may indicate poor data quality.")
            EC50_norm_bq = brentq(hill_eq_brentq, -1.0, 2.0, args=(hill_constants, ec50_normalized))
            EC50_calculable = True
        except ValueError:
            EC50_calculable = False
            print("ValueError encountered in brentq, even with wider scope! EC50 is not calculable. "
                  "Sample : %s_%s" % (sample_name))

    return EC50_norm_bq, EC50_calculable


def _cal_ec50_normalized(temp_data, dose_response_curveshape, method_calc_reading_50, sample):
    # percentage_response = settings["percentage_response"]
    percentage_response = 50
    # calculate fraction response (i.e. convert EC50 to 0.5)
    fract_response = percentage_response / 100
    temp_data["curve_min"] = {}
    temp_data["curve_max"] = {}
    temp_data["y50_curve_min"] = {}
    temp_data["y50_curve_max"] = {}
    if dose_response_curveshape == "S":
        temp_data["curve_min"]["value"] = temp_data["dose"]["normalized"][0]
        temp_data["curve_max"]["value"] = temp_data["dose"]["normalized"][-1]
    elif dose_response_curveshape == "Z":
        temp_data["curve_min"]["value"] = temp_data["dose"]["normalized"][-1]
        temp_data["curve_max"]["value"] = temp_data["dose"]["normalized"][0]

    if method_calc_reading_50 == "ec50 = (curve_max - curve_min)*0.5 + curve_min":
        # use the hill equation to find the y-value of the curve at these positions
        temp_data["y50_curve_min"]["value"] = hill_eq(temp_data["hill_constants"], temp_data["curve_min"]["value"])
        temp_data["y50_curve_max"]["value"] = hill_eq(temp_data["hill_constants"], temp_data["curve_max"]["value"])

        # define reading_50 (yvalue at top of curve - yvalue at bottom of the curve) * 0.5 [or 0.9 for EC90, etc.]
        y50_normalized = (temp_data["y50_curve_max"]["value"] - temp_data["y50_curve_min"]["value"]) * \
                         fract_response + temp_data["y50_curve_min"]["value"]

    # detect extended curve formula (e.g. "reading_50 = (extendedcurve|0.2|_max - extendedcurve|0.2|_min)*0.5")
    elif "extendedcurve" in method_calc_reading_50 and "|" in method_calc_reading_50:

        x_range = temp_data["curve_max"]["value"] - temp_data["curve_min"]["value"]
        # extract the extension (e.g. 0.2, 20% from the text string in the settings file)
        extension_curve_max = float(method_calc_reading_50.split("|")[1])
        extension_curve_min = float(method_calc_reading_50.split("|")[3])
        if dose_response_curveshape == "S":
            # define the x-value for curve_min as the first x-value, minus the xrange * extension
            temp_data["curve_min"]["value"] = temp_data["dose_normalized"][0] - x_range * extension_curve_min
            # define the x-value for curve_max as the last x-value
            temp_data["curve_max"]["value"] = temp_data["dose_normalized"][-1] + x_range * extension_curve_max
        if dose_response_curveshape == "Z":
            # define the x-value for curve_min as the last x-value
            temp_data["curve_min"]["value"] = temp_data["dose_normalized"][-1] + x_range * extension_curve_max
            # define the x-value for curve_max as the first x-value
            temp_data["curve_max"]["value"] = temp_data["dose_normalized"][0] - x_range * extension_curve_min

        # use the hill equation to find the y-value of the curve at these positions
        temp_data["y50_curve_min"]["value"] = hill_eq(temp_data["hill_constants"], temp_data["curve_min"]["value"])
        temp_data["y50_curve_max"]["value"] = hill_eq(temp_data["hill_constants"], temp_data["curve_max"]["value"])
        # define reading_50 (yvalue at top of curve - yvalue at bottom of the curve) * 0.5 [or 0.9 for EC90, etc.]

        y50_normalized = (temp_data["y50_curve_max"]["value"] - temp_data["y50_curve_min"]["value"]) * fract_response

    elif method_calc_reading_50 == "reading_50 = (resp_max - resp_min)*0.5":
        # define ec50 as the centre between min and max datapoints
        # currently the data is normalised between zero and one, so actually ec50 = fract_response
        y50_normalized = (temp_data["reading_max"] - temp_data["reading_min"]) * fract_response

    elif method_calc_reading_50 == "ec50 = (resp_end - resp_start)*0.5":
        y50_normalized = (temp_data["dose_normalized"][-1] - temp_data["dose_normalized"][0]) * fract_response

    elif method_calc_reading_50 == "reading_50 = (resp_start - resp_end)*0.5":
        y50_normalized = (temp_data["dose_normalized"][0] - temp_data["dose_normalized"][-1]) * fract_response

    temp_data["y50_normalized"] = {}
    temp_data["y50_normalized"]["value"] = y50_normalized

    # the y-value of 50% cell density is calculated as the middle position in the curve
    # if the curve is perfectly symmetrical, the EC50 should equal the constant 'k' from the hill_constants
    temp_data["curve_max_norm"] = {}
    temp_data["curve_min_norm"] = {}
    temp_data["curve_max_norm"]["value"] = max(temp_data["reading"]["fitted_normalized"])
    temp_data["curve_min_norm"]["value"] = min(temp_data["reading"]["fitted_normalized"])

    # dfe.loc["EC50_norm_bq{}".format(d),"%s_okay" % sLet]
    brentq_out_tuple = _calc_EC50_brent_eq(sample, temp_data["hill_constants"], temp_data["y50_normalized"]["value"])
    temp_data["EC50_norm_bq"] = {}
    temp_data["EC50_norm_bq"]["value"], temp_data["EC50_calculable"] = brentq_out_tuple
    # add if the EC50 was calculable to the summary dataframe "okay" column
    if temp_data["EC50_calculable"]:
        temp_data["EC50_norm_bq"]["check"] = True
    else:
        temp_data["EC50_norm_bq"]["check"] = False

    return
