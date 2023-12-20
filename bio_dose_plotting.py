import numpy as np
import matplotlib.pyplot as plt

from bio_dose_auto_judge import JudgeFit
from bio_dose_functions import setup_t20_colour_list, denormalise_0_1


class PlottingDose:
    def __init__(self, config, dose_label="chloramphenicol", dose_units="mg/ml", method="EC50"):
        self.config = config
        self.jf = JudgeFit(config)
        self.dose_label = dose_label
        self.dose_units = dose_units
        self.method = method

        self.annotation_font_size = 10
        self.fig_font_size = 6
        # set the default font for the figures
        plt.subplots()
        plt.rc('font', family='sans-serif')
        plt.rc('font', serif='Helvetica Neue')
        plt.rc('text', usetex='false')
        plt.rcParams.update({'font.size': self.fig_font_size})
        plt.close("all")
        self.fig, self.axarr = plt.subplots(nrows=2, ncols=1)
        # create a colour list for the various datasets, with selected colours from the tableau20
        t20 = setup_t20_colour_list()
        self.colour_list = []
        for colour_counter, _ in enumerate(t20):
            self.colour_list.append(t20[colour_counter])
        # create dictionary to hold the size of the plotted datapoints for each dataset (gradually smaller)
        self.size_datapoints = {}
        self.lines_tyles = {}
        self.all_line_styles = ["dashed", "dotted", "-."]
        self.line_style_counter = 0
        self.alpha_transparency_level = 0.8
        self.xy_coordinates = "axes fraction"

        self.y_limit_min_raw = None
        self.x_limit_min_raw = None
        self.y_limit_max_raw = None
        self.x_limit_max_plot1 = None

    def plot_raw_fig1(self, counter, sample, temp_data):

        self.axarr[0].scatter(temp_data["dose"]["raw"], temp_data["reading"]["raw"],
                              color=self.colour_list[counter], s=self.size_datapoints[sample], label=sample,
                              alpha=self.alpha_transparency_level)

        # set xlabel, ylabel, title, etc
        self.axarr[0].set_xlabel(f"{self.dose_label} ({self.dose_units})", fontsize=self.fig_font_size)
        # axarr[0,0].set_ylabel(settings["y-axis (response) label"],rotation='vertical', fontsize = fig_font_size)
        self.axarr[0].set_ylabel("cell density (A600)", rotation='vertical', fontsize=self.fig_font_size)
        self.axarr[0].set_title(f"Raw Data", fontsize=self.fig_font_size)
        self.axarr[0].grid(True, color='0.75')

        ymin, ymax = self.axarr[0].get_ylim()

        # define the y-limit for the non-normalised (raw) data as the minimum value minus a percentage of max value
        self.y_limit_min_raw = ymin - ymax * 0.05
        # define the x-limit for the non-normalised (raw) data as the minimum value minus a percentage of max value

        self.x_limit_min_raw = min(temp_data["dose"]["raw"]) - max(temp_data["dose"]["raw"]) * 0.1
        self.y_limit_max_raw = max(temp_data["reading"]["raw"]) + 0.1 if max(temp_data["reading"]["raw"]) > 1.0 else 1.0
        self.axarr[0].set_ylim(self.y_limit_min_raw, self.y_limit_max_raw)
        # define 110% of the limit of the x-axis as the maximum dose conc.
        self.x_limit_max_plot1 = max(temp_data["dose"]["raw"]) * 1.1
        self.axarr[0].set_xlim(self.x_limit_min_raw, self.x_limit_max_plot1)

        temp_data["dose"]["fitted"] = denormalise_0_1(temp_data["reading"]["fitted_normalized"],
                                                      temp_data["reading"]["min"], temp_data["reading"]["max"])
        temp_data["reading"]["fitted"] = denormalise_0_1(temp_data["dose"]["fitted_normalized"],
                                                         temp_data["dose"]["min"], temp_data["dose"]["max"])

    def plot_fitted_fig1(self, counter, sample, temp_data):

        self.axarr[0].scatter(temp_data["dose"]["raw"], temp_data["reading"]["raw"], color=self.colour_list[counter],
                              s=self.size_datapoints[sample], label=sample, alpha=self.alpha_transparency_level)
        # plot fitted curve on the subplot with the original y values
        self.axarr[0].plot(temp_data["dose"]["fitted"], temp_data["reading"]["fitted"], '-',
                           color=self.colour_list[counter], alpha=self.alpha_transparency_level)
        # extract y50 and EC50 from dataframe
        y50 = temp_data["y50"]["value"]
        EC50 = temp_data["EC50"]["value"]
        # draw horizontal line from y50 to EC50
        self.axarr[0].hlines(y=y50, xmin=self.x_limit_min_raw, xmax=EC50, colors=self.colour_list[counter],
                             linestyles=self.lines_tyles[sample], label='', alpha=self.alpha_transparency_level)
        # draw vertical line at EC50 from y50
        self.axarr[0].vlines(x=EC50, ymin=self.y_limit_max_raw, ymax=y50, colors=self.colour_list[counter],
                             linestyles=self.lines_tyles[sample])

    def plot_normalised_fig2(self, counter, sample, temp_data):
        # set y-axis intercept
        ymin_norm = -0.2
        xmin_norm = -0.05

        # add normalised datapoints as a scattergram
        self.axarr[1].scatter(temp_data["dose"]["normalized"], temp_data["reading"]["normalized"],
                              color=self.colour_list[counter], s=self.size_datapoints[sample], label=sample)
        # add the fitted curve as a line plot
        self.axarr[1].plot(temp_data["dose"]["fitted"], temp_data["reading"]["fitted"],
                           '-', color=self.colour_list[counter], alpha=self.alpha_transparency_level)
        # add horizontal line at y50
        self.axarr[1].hlines(y=temp_data["y50_normalized"]["value"], xmin=xmin_norm, colors=self.colour_list[counter],
                             xmax=temp_data["EC50_norm_bq"]["value"], linestyles=self.lines_tyles[sample])

        # add vertical line at EC50
        self.axarr[1].vlines(x=temp_data["EC50_norm_bq"]["value"], ymin=ymin_norm, colors=self.colour_list[counter],
                             ymax=temp_data["y50_normalized"]["value"], linestyles=self.lines_tyles[sample])

        # set xlabel, ylabel, title, grid, etc
        self.axarr[1].set_xlabel("dose concentration (normalised)",
                                 fontsize=self.fig_font_size)
        self.axarr[1].set_ylabel("response concentration (normalised)", rotation='vertical',
                                 fontsize=self.fig_font_size)
        self.axarr[1].text(0.6, 1.1, "normalised data", horizontalalignment='center',
                           fontsize=self.fig_font_size)
        self.axarr[1].grid(True, color='0.75')
        self.axarr[1].set_ylim(ymin_norm, 1.2)
        self.axarr[1].set_xlim(xmin_norm, 1.2)

    def plot_automatic_judgement_fig4(self, temp_data):

        if temp_data["EC50_calculable"]:
            # add the stepsize near the EC50, which determines whether more dose concentrations are necessary
            stepcolour = temp_data["doseconc_stepsize_at_EC50"]["colour"]

            if stepcolour != "k":

                self.axarr[1].plot(temp_data["doseconc_steps_at_EC50"], (0, 0), color=stepcolour, linestyle="-", lw=2)

            test_max_low_dose_slop = self.config["dose"].getfloat("max_lowdose_slope")
            test_max_high_dose_slop = self.config["dose"].getfloat("max_highdose_slope")

            # if the slope at the lowdose is above the chosen cutoff, draw a line on the normalised plot, axarr[1,0]
            if temp_data["saxe_lowdose"]["value"] > test_max_low_dose_slop:

                saxe_lowdose_values = temp_data["saxe_lowdose_values"]["value"]
                # draw red vertical line showing the slope at the lowdose datapoint
                self.axarr[1, 0].plot(saxe_lowdose_values[0], saxe_lowdose_values[1], 'r-', lw=2)
            # if the slope at the highdose is higher than the chosen cutoff, draw a line on tho normalised plot
            if temp_data["saxe_highdose"]["value"] > test_max_high_dose_slop:
                saxe_highdose_values = temp_data["saxe_highdose_values"]["value"]
                # draw red vertical line showing the slope at the lowdose datapoint
                self.axarr[1, 0].plot(saxe_highdose_values[0], saxe_highdose_values[1], 'r-', lw=2)

    def _setup(self, sample_data):
        self.line_style_counter = 0
        for counter, sample in enumerate(sample_data):

            size_start = 15
            size_increment = -7
            self.size_datapoints[sample] = size_start + size_increment * counter
            self.lines_tyles[sample] = self.all_line_styles[self.line_style_counter]
            self.line_style_counter += 1
            if self.line_style_counter >= len(self.all_line_styles):
                self.line_style_counter = 0

    def controller(self, sample_data):
        self._setup(sample_data)

        yaxis_pos = np.linspace(0.9, 0.1, 10)
        xaxis_left = 0.05

        for counter, sample in enumerate(sample_data):
            if sample != "state_data":
                self.plot_raw_fig1(counter, sample, sample_data[sample])
                self.plot_fitted_fig1(counter, sample, sample_data[sample])
                self.plot_normalised_fig2(counter, sample, sample_data[sample])
                # analyse the curve fit and data to judge whether the EC50 value is accurate
                self.jf.judge_controller(sample_data[sample])
                self.plot_automatic_judgement_fig4(sample_data[sample])

        # plt.legend()
        # plt.show()
        return self.fig, self.axarr


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    plt.subplots()
