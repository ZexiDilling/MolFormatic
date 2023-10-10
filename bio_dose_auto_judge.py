import configparser

import numpy as np

from bio_dose_functions import hill_eq


class JudgeFit:
    def __init__(self, config):

        # (e.g. max_std_resp_highdose_dp).
        self.n_neighb = config["dose"].getint("n_neighb")
        # maximum standard deviation of the response datapoints at high dose concentration
        self.max_std_resp_highdose_dp = config["dose"].getfloat("max_std_resp_highdose_dp")
        # maximum standard deviation of the response datapoints at low dose concentration
        self.max_std_resp_lowdose_dp = config["dose"].getfloat("max_std_resp_lowdose_dp")
        self.min_flat_lowdose_dp = config["dose"].getfloat("min_flat_lowdose_dp")
        self.min_flat_highdose_dp = config["dose"].getfloat("min_flat_highdose_dp")

        # minimum rsquared of the fit from sigmoidal curve to the data
        self.min_rsquared = config["dose"].getfloat("min_rsquared")
        # minimum acceptable dase concentration stepsizes. Smaller stepsizes give more accurate EC50 values!
        self.max_acceptable_doseconc_stepsize_at_EC50 = config["dose"].getfloat("max_acceptable_doseconc_stepsize_at_EC50")
        self.max_recommended_doseconc_stepsize_at_EC50 = config["dose"].getfloat("max_recommended_doseconc_stepsize_at_EC50")
        # minimum hillslope of the fit from sigmoidal curve to the data (below 1, tends not to be sigmoidal)
        self.weak_hillslope_range = [-1, 1]
        # ToDo get from cofig
        # minimum value for the end of the curve, on the y-axis (below -1, tends not to be sigmoidal)
        self.min_curve_lowresp = config["dose"].getfloat("min_curve_lowresp")

        self.width_lowdose_slope = config["dose"].getfloat("width_lowdose_slope")
        self.width_highdose_slope = config["dose"].getfloat("width_highdose_slope")
        self.max_lowdose_slope = config["dose"].getfloat("max_lowdose_slope")
        self.max_highdose_slope = config["dose"].getfloat("max_highdose_slope")

    def _set_up(self):

        self.temp_data["indices_highdose_datapoints_excl_nearest_EC50"] = {}
        self.temp_data["response_highdose_datapoints"] = {}
        self.temp_data["n_highdose_datapoints"] = {}
        self.temp_data["response_lowdose_datapoints"] = {}
        self.temp_data["n_lowdose_datapoints"] = {}
        self.temp_data["std_resp_highdose_datapoints"] = {}
        self.temp_data["std_resp_lowdose_datapoints"] = {}
        self.temp_data["doseconc_stepsize_at_EC50"] = {}

        self.x = np.array(self.temp_data["dose"]["raw"])
        self.y = np.array(self.temp_data["reading"]["raw"])
        # identify the datapoints at high dose concentrations

        self.temp_data["indices_highdose_datapoints"] = np.where(self.x > self.temp_data["EC50"]["value"])[0]

        # remove the datapoint closest to EC50
        self.temp_data["indices_highdose_datapoints_excl_nearest_EC50"]["value"] = \
            self.temp_data["indices_highdose_datapoints"][self.n_neighb:]

        # slice using the indices to yield the OD600 values for the highdose datapoints
        self.temp_data["response_highdose_datapoints"]["value"] =  \
            self.y[self.temp_data["indices_highdose_datapoints_excl_nearest_EC50"]["value"]]

        # count the number of highdose datapoint
        self.temp_data["n_highdose_datapoints"]["value"] = len(self.temp_data["response_highdose_datapoints"]["value"])

        # identify the lowdose datapoints, count and measure standard deviation
        # identify the lowdose datapoints (x < EC50)
        self.temp_data["indices_lowdose_datapoints"] = np.where(self.x < self.temp_data["EC50"]["value"])[0]

        # exclude datapoint closest to the EC50
        self.temp_data["indices_lowdose_datapoints_excl_nearest_EC50"] = \
            self.temp_data["indices_lowdose_datapoints"][self.n_neighb:]

        # use index to select the y-axis (response) data
        self.temp_data["response_lowdose_datapoints"]["value"] = \
            self.y[self.temp_data["indices_lowdose_datapoints_excl_nearest_EC50"]]

        # count the datapoints
        self.temp_data["n_lowdose_datapoints"]["value"] = len(self.temp_data["response_lowdose_datapoints"]["value"])

    def _check_lower_bound(self):
        # judge whether the data contains enough high and lowdose datapoints
        if self.temp_data["n_highdose_datapoints"]["value"] >= self.min_flat_highdose_dp:
            self.temp_data["n_highdose_datapoints"]["check"] = True
            self.temp_data["n_highdose_datapoints"]["colour"] = "k"
        else:
            self.temp_data["n_highdose_datapoints"]["check"] = False
            self.temp_data["n_highdose_datapoints"]["colour"] = "r"

        # evaluate as "okay" if number of highdose or lowdose datapoints is more than two
        if self.temp_data["n_lowdose_datapoints"]["value"] >= self.min_flat_lowdose_dp:
            self.temp_data["n_lowdose_datapoints"]["check"] = True
            self.temp_data["n_lowdose_datapoints"]["colour"] = "k"
        else:
            self.temp_data["n_lowdose_datapoints"]["check"] = False
            self.temp_data["n_lowdose_datapoints"]["colour"] = "r"

    def _check_response_bound(self):




        # judge whether the standard deviation of the high and lowdose datapoints is acceptable
        if self.temp_data["n_highdose_datapoints"]["value"] > 1:
            # calculate std of highdose datapoints
            self.temp_data["std_resp_highdose_datapoints"]["value"] = \
                np.std(self.temp_data["response_highdose_datapoints"]["value"])
            # evaluate as "okay" if std of highdose datapoints is less than a cutoff value (max_std_resp_highdose_dp)
            if self.temp_data["std_resp_highdose_datapoints"]["value"] < self.max_std_resp_highdose_dp:
                self.temp_data["std_resp_highdose_datapoints"]["check"] = True
                self.temp_data["std_resp_highdose_datapoints"]["colour"] = 'k'
            else:
                self.temp_data["std_resp_highdose_datapoints"]["check"] = False
                self.temp_data["std_resp_highdose_datapoints"]["colour"] = 'r'
        else:
            # there is either insufficient lowresponse or highresponse datapoints(insuff_lowresp_dp,
            # or insuff_highresp_dp).
            # Replace std with 0, and colour black on the figure.
            self.temp_data["std_resp_highdose_datapoints"]["value"] = 0
            self.temp_data["std_resp_highdose_datapoints"]["colour"] = 'k'

        if self.temp_data["n_lowdose_datapoints"]["value"] > 1:
            # calculate std of lowdose datapoints
            self.temp_data["std_resp_lowdose_datapoints"]["value"] = \
                np.std(self.temp_data["response_lowdose_datapoints"]["value"])
            # evaluate as "okay" if std of lowdose datapoints is less than a cutoff value
            if self.temp_data["std_resp_lowdose_datapoints"]["value"] < self.max_std_resp_lowdose_dp:
                self.temp_data["std_resp_lowdose_datapoints"]["check"] = True
                self.temp_data["std_resp_lowdose_datapoints"]["colour"] = 'k'
            else:
                self.temp_data["std_resp_lowdose_datapoints"]["check"] = False
                self.temp_data["std_resp_lowdose_datapoints"]["colour"] = 'r'
        else:
            # there is either insufficient lowresponse or highresponse datapoints(insuff_lowresp_dp, or insuff_highresp_dp).
            # Replace std with 0, and colour black.
            self.temp_data["std_resp_lowdose_datapoints"]["value"] = 0
            self.temp_data["std_resp_lowdose_datapoints"]["colour"] = 'k'

    def _check_step_size(self):

        # identify the tested dose concentration below the EC50
        indices_lowdose_datapoints = np.where(self.x < self.temp_data["EC50"]["value"])[0]
        if indices_lowdose_datapoints.size != 0:
            doseconc_before_EC50 = self.x[indices_lowdose_datapoints[-1]]
            # identify the tested dose concentration after the EC50

            doseconc_after_EC50 = self.x[self.temp_data["indices_highdose_datapoints"][0]]
            # add values to output dataframe, so that the plot can be annotated
            self.temp_data["doseconc_steps_at_EC50"] = (doseconc_before_EC50, doseconc_after_EC50)
            # calculate the stepsize at the EC50. Smaller is better!
            doseconc_stepsize_at_EC50 = doseconc_after_EC50 - doseconc_before_EC50
            self.temp_data["doseconc_stepsize_at_EC50"]["value"] = doseconc_stepsize_at_EC50
            # evaluate as "okay" if the stepsize at the EC50 is smaller than the min acceptable value
            if doseconc_stepsize_at_EC50 <= self.max_acceptable_doseconc_stepsize_at_EC50:
                self.temp_data["doseconc_stepsize_at_EC50"]["check"] = True
                # if the stepsize is small, colour to dark red as a warning that the doseconc should be optimised
                if doseconc_stepsize_at_EC50 <= self.max_recommended_doseconc_stepsize_at_EC50:
                    self.temp_data["doseconc_stepsize_at_EC50"]["colour"] = 'k'
                else:
                    # colour dark red
                    self.temp_data["doseconc_stepsize_at_EC50"]["colour"] = '#990033'
            else:
                # the stepsize is extremely high, and the data therefore has little value. doseconc needs to be optimised.
                self.temp_data["doseconc_stepsize_at_EC50"]["check"] = False
                self.temp_data["doseconc_stepsize_at_EC50"]["colour"] = 'r'
        else:
            # there is either insufficient lowresponse or highresponse datapoints(insuff_lowresp_dp, or insuff_highresp_dp).
            # Stepsize can't be calculated. Replace with 0, and colour grey.
            self.temp_data["doseconc_stepsize_at_EC50"]["check"] = False
            self.temp_data["doseconc_steps_at_EC50"] = (0, 0)
            self.temp_data["doseconc_stepsize_at_EC50"]["value"] = 0
            self.temp_data["doseconc_stepsize_at_EC50"]["colour"] = '0.5'

    def _check_r_squard(self):
        if self.temp_data["rsquared"]["value"] > self.min_rsquared:
            self.temp_data["rsquared"]["check"] = True
            self.temp_data["rsquared"]["colour"] = 'k'
        else:
            self.temp_data["rsquared"]["check"] = False
            self.temp_data["rsquared"]["colour"] = 'r'

    def _check_hillslope(self):
        """
         Does the hillslope parameter from the fitted curve approach 0?
        :return:
        """

        if self.weak_hillslope_range[0] < self.temp_data["hillslope"]["value"] < self.weak_hillslope_range[1]:
            # if it's outside the range, label as "data_needs_checking"
            self.temp_data["hillslope"]["check"] = False
            self.temp_data["hillslope"]["colour"] = 'r'
        else:
            # if it's outside the range, label as okay
            self.temp_data["hillslope"]["check"] = True
            self.temp_data["hillslope"]["colour"] = 'k'

    def _check_curve_start_at_negative(self):
        if self.temp_data["curve_min_norm"]["value"] > self.min_curve_lowresp:
            self.temp_data["curve_min_norm"]["check"] = True
            self.temp_data["curve_min_norm"]["colour"] = 'k'
        else:
            self.temp_data["curve_min_norm"]["check"] = False
            self.temp_data["curve_min_norm"]["colour"] = 'r'

    def _check_saxe(self):
        """
        Does the Slope At X-axis Extremes (SAXE) exceed a user-defined threshold?
        1) calculate average stepsize
        2) define points left and right of the xaxis extremes
        3) calculate slope between points
        4) determine if slope is shallow enough to indicate an S- or Z-shaped curve
        :return:
        """
        hill_constants = self.temp_data["hill_constants"]
        xnorm = self.temp_data["reading"]["normalized"]
        # find the average stepsize by creating two truncated arrays from xnorm, and subtracting
        xnorm_left = np.array(list(xnorm)[:-1])
        xnorm_right = np.array(list(xnorm)[1:])
        xnorm_stepsizes = xnorm_right - xnorm_left
        xnorm_stepsize_mean = xnorm_stepsizes.mean()
        # define the width surrounding the datapoint for the slope measurement
        # calculated as the mean stepsize multiplied by a user value (0.001 to 1.0)
        width_lowdose_slope = xnorm_stepsize_mean / 2 * self.width_lowdose_slope
        width_highdose_slope = xnorm_stepsize_mean / 2 * self.width_highdose_slope
        # define the min and max of normalised datapoints (will simply be 0 and 1 for normalised data)
        xnorm_min, xnorm_max = xnorm.min(), xnorm.max()
        # define SAXE lowdose/highdose x-axis datapoints (to the left and right of the original datapoints)
        saxe_lowdose_x_dp_left = xnorm_min - width_lowdose_slope
        # if it is negative (results in nan in the sigmoidal function), replace with 0
        saxe_lowdose_x_dp_left = saxe_lowdose_x_dp_left if saxe_lowdose_x_dp_left > 0 else 0
        saxe_lowdose_x_dp_right = xnorm_min + width_lowdose_slope
        saxe_highdose_x_dp_left = xnorm_max - width_highdose_slope
        saxe_highdose_x_dp_right = xnorm_max + width_highdose_slope
        # calculate the y-values on the curve, for the x-values surrounding the min and max datapoints
        saxe_lowdose_y_dp_left = hill_eq(hill_constants, saxe_lowdose_x_dp_left)
        saxe_lowdose_y_dp_right = hill_eq(hill_constants, saxe_lowdose_x_dp_right)
        saxe_highdose_y_dp_left = hill_eq(hill_constants, saxe_highdose_x_dp_left)
        saxe_highdose_y_dp_right = hill_eq(hill_constants, saxe_highdose_x_dp_right)
        # calculate the linear slope (y2 - y1)/(x2 - x1) between the chosen datapoints from the curve
        saxe_lowdose_1 = (saxe_lowdose_y_dp_right - saxe_lowdose_y_dp_left) / (
                    saxe_lowdose_x_dp_right - saxe_lowdose_x_dp_left)
        saxe_highdose_1 = (saxe_highdose_y_dp_right - saxe_highdose_y_dp_left) / (
                    saxe_highdose_x_dp_right - saxe_highdose_x_dp_left)
        # convert slopes to positive numbers
        saxe_lowdose = abs(saxe_lowdose_1)
        saxe_highdose = abs(saxe_highdose_1)
        # add to output dataframe, dfe
        self.temp_data["saxe_lowdose"] = {}
        self.temp_data["saxe_highdose"] = {}

        self.temp_data["saxe_lowdose"]["value"] = saxe_lowdose
        self.temp_data["saxe_highdose"]["value"] = saxe_highdose
        self.temp_data["saxe_lowdose"]["values"] = [[saxe_lowdose_x_dp_left, saxe_lowdose_x_dp_right],
                                                            [saxe_lowdose_y_dp_left, saxe_lowdose_y_dp_right]]
        self.temp_data["saxe_highdose"]["values"] = [[saxe_highdose_x_dp_left, saxe_highdose_x_dp_right],
                                                             [saxe_highdose_y_dp_left, saxe_highdose_y_dp_right]]

        # check that the calculated slopes of the curve do not exceed the saxe_max_slope
        if saxe_lowdose < self.max_lowdose_slope:
            self.temp_data["saxe_lowdose"]["check"] = True
            self.temp_data["saxe_lowdose"]["colour"] = 'k'
        else:
            self.temp_data["saxe_lowdose"]["check"] = False
            self.temp_data["saxe_lowdose"]["colour"] = 'r'

        if saxe_highdose < self.max_highdose_slope:
            self.temp_data["saxe_highdose"]["check"] = True
            self.temp_data["saxe_highdose"]["colour"] = 'k'
        else:
            self.temp_data["saxe_highdose"]["check"] = False
            self.temp_data["saxe_highdose"]["colour"] = 'r'

    def judge_controller(self, temp_data):
        self.temp_data = temp_data
        self._set_up()
        self._check_lower_bound()
        self._check_response_bound()
        self._check_step_size()
        self._check_r_squard()
        self._check_hillslope()
        self._check_curve_start_at_negative()
        self._check_saxe()
        temp_data["approved"] = True        #ToDo make this check work


if __name__ == "__main__":
    config = configparser.ConfigParser()
    config.read("config.ini")

    jf = JudgeFit(config)

