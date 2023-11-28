import configparser
from pathlib import Path

import numpy as np
from scipy.optimize import leastsq

from bio_dose_excle_handler import dose_excel_controller
from bio_dose_functions import denormalise_0_1, _calc_EC50_brent_eq, _cal_ec50_normalized, _calc_residuals_mean, \
    _calc_rsquared, normalise_0_1, residuals, hill_eq
from bio_dose_plotting import PlottingDose


def _set_up(temp_data, hill_constants_guess):
    # make an array of >250 datapoints representing the x-axis of the curve
    temp_min = 0
    # max = settings["fitted_curve_xaxis_max"]
    # n_datapoints = settings["fitted_curve_n_datapoints"]
    temp_max = 1.2
    temp_n_datapoints = 1000
    temp_data["dose"]["fitted_normalized"] = np.linspace(temp_min, temp_max, temp_n_datapoints)
    temp_data["n_doseconc_tested"] = len(temp_data["dose"]["raw"])
    hill_constants, cov, infodict, mesg, ier = leastsq(residuals, hill_constants_guess,
                                                       args=(hill_eq, temp_data["dose"]["normalized"],
                                                             temp_data["reading"]["normalized"]),
                                                       full_output=1)
    temp_data["hill_constants"] = list(hill_constants)
    temp_data["hillslope"] = {}
    temp_data["hillslope"]["value"] = temp_data["hill_constants"][3]
    temp_data["EC50_hill_eq_norm"] = temp_data["hill_constants"][2]
    temp_data["rsquared"] = {}
    temp_data["rsquared"]["value"] = _calc_rsquared(infodict["fvec"], temp_data["reading"]["normalized"])
    temp_data["residuals_normalized_mean"] = {}
    temp_data["residuals_normalized_mean"]["value"] = _calc_residuals_mean(hill_constants, temp_data)

    # denormalise
    xmin = temp_data["dose"]["min"]
    xmax = temp_data["dose"]["max"]
    temp_data["EC50_hill_eq"] = denormalise_0_1(hill_constants[2], xmin, xmax)

    # calculate the value for y for the 1500 points
    x_fitted_norm = np.array(temp_data["dose"]["fitted_normalized"])
    temp_data["reading"]["fitted_normalized"] = hill_eq(hill_constants, x_fitted_norm)


def dose_response_controller(config, dose_response_curveshape, all_data, method_calc_reading_50):
    pd = PlottingDose(config, all_data)

    if dose_response_curveshape == "S":
        hill_constants_guess = (0.0, 1.0, 0.5, 10.0)
    elif dose_response_curveshape == "Z":
        hill_constants_guess = (1.0, 0.0, 0.5, 10.0)

    for samples in all_data:
        if samples != "state_data":
            all_data[samples]["reading"]["normalized"], all_data[samples]["reading"]["min"], \
            all_data[samples]["reading"]["max"] = normalise_0_1(all_data[samples]["reading"]["raw"])
            all_data[samples]["dose"]["normalized"], all_data[samples]["dose"]["min"], all_data[samples]["dose"]["max"] = \
                normalise_0_1(all_data[samples]["dose"]["raw"])

            _set_up(all_data[samples], hill_constants_guess)

            _cal_ec50_normalized(all_data[samples], dose_response_curveshape,
                                                   method_calc_reading_50, samples)

            brentq_out_tuple = _calc_EC50_brent_eq(samples, all_data[samples]["hill_constants"],
                                                   all_data[samples]["y50_normalized"]["value"])


            all_data[samples]["EC50"] = {}
            all_data[samples]["EC50"]["value"] = float(
                denormalise_0_1(brentq_out_tuple[0], all_data[samples]["dose"]["min"], all_data[samples]["dose"]["max"]))
            # dfe.loc["y_fitted{}".format(d), sLet] = temp_y_
            # denormalise the y50, the y-value used to calculated the EC50
            all_data[samples]["y50"] = {}
            all_data[samples]["y50"]["value"] = denormalise_0_1(all_data[samples]["y50_normalized"]["value"], all_data[samples]["reading"]["min"],
                                                                all_data[samples]["reading"]["max"])
    print("calc done !!")
    all_data = pd.controller(save_location="save_locations")

    return all_data


if __name__ == "__main__":

    config = configparser.ConfigParser()
    config.read("config.ini")

    dose_response_curveshape = "S"
    all_data = {"control": {
        "reading": {"raw": [107.9699121, 103.5144308, 104.0479892, 104.3916889, 97.27653273, 90.69159613, 56.65732349,
                    37.87317223, 19.78570689, 8.299916445, 11.36102686, 4.144031277, 5.013870229, 3.04536459,
                    4.692352962, 7.622099516, 5.336804101, 9.853183031, 9.120045299, 4.869276415, 6.108990371,
                    3.597853303, 5.994511486, 5.717414715]},
        "dose": {"raw": [0, 0.043478261, 0.086956522, 0.130434783, 0.173913043, 0.217391304, 0.260869565, 0.304347826,
                 0.347826087, 0.391304348, 0.434782609, 0.47826087, 0.52173913, 0.565217391, 0.608695652, 0.652173913,
                 0.695652174, 0.739130435, 0.782608696, 0.826086957, 0.869565217, 0.913043478, 0.956521739, 1]}

    },
    "sample":{
        "reading": {"raw": [109.3214119, 108.6611434, 102.3685731, 103.9247113, 104.9325289, 102.4206455, 98.47691972,
                            102.5875971, 93.42962519, 76.20087959, 60.9377449, 45.92929776, 27.13356751, 15.82598671,
                            12.67835439, 13.02441613, 10.35963326, 5.007706181, 6.146401995, 3.79602644, 4.418804534,
                            9.61913024, 3.43543862, 5.886764018]},
        "dose": {"raw": [0, 0.043478261, 0.086956522, 0.130434783, 0.173913043, 0.217391304, 0.260869565, 0.304347826,
                         0.347826087, 0.391304348, 0.434782609, 0.47826087, 0.52173913, 0.565217391, 0.608695652,
                         0.652173913, 0.695652174, 0.739130435, 0.782608696, 0.826086957, 0.869565217, 0.913043478,
                         0.956521739, 1]}
    }}
    method_calc_reading_50 = "ec50 = (curve_max - curve_min)*0.5 + curve_min"
    all_dose_data = {}

    all_dose_data["plate"] = dose_response_controller(config, dose_response_curveshape, all_data, method_calc_reading_50)
    save_location = Path(r"C:\Users\phch\Desktop\test\dose_response")
    plate_group_to_compound_id = None
    include_id = False
    dose_excel_controller(all_dose_data, plate_group_to_compound_id, include_id, save_location)
