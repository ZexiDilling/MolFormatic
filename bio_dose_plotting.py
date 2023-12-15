import numpy as np

from bio_dose_auto_judge import JudgeFit
from bio_dose_functions import setup_t20_colour_list, denormalise_0_1


class PlottingDose:
    def __init__(self, config, all_data, doselabel="chloramphenicol", doseunits="mg/ml", method="EC50"):
        self.all_data = all_data
        self.jf = JudgeFit(config)
        self.doselabel = doselabel
        self.doseunits = doseunits
        self.method = method

        # set an annotation fontsize
        self.annotation_fontsize = 10
        # set the fontsize for the figure
        self.fig_fontsize = 6
        # set the default font for the figures
        plt.rc('font', family='sans-serif')
        plt.rc('font', serif='Helvetica Neue')
        plt.rc('text', usetex='false')
        plt.rcParams.update({'font.size': self.fig_fontsize})
        plt.close("all")
        self.fig, self.axarr = plt.subplots(nrows=2, ncols=2)
        # create a colour list for the various datasets, with selected colours from the tableau20
        t20 = setup_t20_colour_list()
        self.colour_list = []
        for colour_counter, _ in enumerate(t20):
            self.colour_list.append(t20[colour_counter])
        # create dictionary to hold the size of the plotted datapoints for each dataset (gradually smaller)
        self.size_datapoints = {}
        # create dictionary to hold the linestyles for the EC50 locators
        self.linestyles = {}
        all_line_styles = ["dashed", "dotted", "-."]
        line_style_counter = 0

        # set alpha (transparency) level for plotted lines and data in general
        self.alpha_transparency_level = 0.8
        # define xycoordinates for later annotations
        self.xy_coordinates = "axes fraction"

        for counter, sample in enumerate(self.all_data):
            size_start = 15
            size_increment = -7
            self.size_datapoints[sample] = size_start + size_increment * counter
            self.linestyles[sample] = all_line_styles[line_style_counter]
            line_style_counter += 1
            if line_style_counter >= len(all_line_styles):
                line_style_counter = 0

        self.y_limit_min_raw = None
        self.x_limit_min_raw = None
        self.y_limit_max_raw = None
        self.x_limit_max_plot1 = None

    def plot_raw_fig1(self, counter, sample, temp_data):

        self.axarr[0, 0].scatter(temp_data["dose"]["raw"], temp_data["reading"]["raw"],
                                 color=self.colour_list[counter], s=self.size_datapoints[sample], label=sample,
                                 alpha=self.alpha_transparency_level)

        # set xlabel, ylabel, title, etc
        self.axarr[0, 0].set_xlabel(f"{self.doselabel} ({self.doseunits})", fontsize=self.fig_fontsize)
        # axarr[0,0].set_ylabel(settings["y-axis (response) label"],rotation='vertical', fontsize = fig_fontsize)
        self.axarr[0, 0].set_ylabel("cell density (A600)", rotation='vertical', fontsize=self.fig_fontsize)
        self.axarr[0, 0].set_title(f"SET_TITLE!", fontsize=self.fig_fontsize)
        self.axarr[0, 0].grid(True, color='0.75')

        ymin, ymax = self.axarr[0, 0].get_ylim()

        # define the y-limit for the non-normalised (raw) data as the minimum value minus a percentage of max value
        self.y_limit_min_raw = ymin - ymax * 0.05
        # define the x-limit for the non-normalised (raw) data as the minimum value minus a percentage of max value

        self.x_limit_min_raw = min(temp_data["dose"]["raw"]) - max(temp_data["dose"]["raw"]) * 0.1
        self.y_limit_max_raw = max(temp_data["reading"]["raw"]) + 0.1 if max(temp_data["reading"]["raw"]) > 1.0 else 1.0
        self.axarr[0, 0].set_ylim(self.y_limit_min_raw, self.y_limit_max_raw)
        # define 110% of the limit of the x-axis as the maximum dose conc.
        self.x_limit_max_plot1 = max(temp_data["dose"]["raw"]) * 1.1
        self.axarr[0, 0].set_xlim(self.x_limit_min_raw, self.x_limit_max_plot1)

        temp_data["dose"]["fitted"] = denormalise_0_1(temp_data["reading"]["fitted_normalized"], temp_data["reading"]["min"], temp_data["reading"]["max"])
        temp_data["reading"]["fitted"] = denormalise_0_1(temp_data["dose"]["fitted_normalized"], temp_data["dose"]["min"], temp_data["dose"]["max"])

    def plot_fitted_fig1(self, counter, sample, temp_data):

        self.axarr[0, 0].scatter(temp_data["dose"]["raw"], temp_data["reading"]["raw"], color=self.colour_list[counter],
                                 s=self.size_datapoints[sample], label=sample, alpha=self.alpha_transparency_level)
        # plot fitted curve on the subplot with the original y values
        self.axarr[0, 0].plot(temp_data["dose"]["fitted"], temp_data["reading"]["fitted"], '-', color=self.colour_list[counter],
                              alpha=self.alpha_transparency_level)
        # extract y50 and EC50 from dataframe
        y50 = temp_data["y50"]["value"]
        EC50 = temp_data["EC50"]["value"]
        # draw horizontal line from y50 to EC50
        self.axarr[0, 0].hlines(y=y50, xmin=self.x_limit_min_raw, xmax=EC50, colors=self.colour_list[counter],
                                linestyles=self.linestyles[sample], label='', alpha=self.alpha_transparency_level)
        # draw vertical line at EC50 from y50
        self.axarr[0, 0].vlines(x=EC50, ymin=self.y_limit_max_raw, ymax=y50, colors=self.colour_list[counter],
                                linestyles=self.linestyles[sample])


    def plot_normalised_fig2(self, counter, sample, temp_data):
        # set y-axis intercept
        ymin_norm = -0.2
        xmin_norm = -0.05

        # add normalised datapoints as a scattergram
        self.axarr[1, 0].scatter(temp_data["dose"]["normalized"], temp_data["reading"]["normalized"],
                                 color=self.colour_list[counter], s=self.size_datapoints[sample], label=sample)
        # add the fitted curve as a line plot
        self.axarr[1, 0].plot(temp_data["dose"]["fitted"], temp_data["reading"]["fitted"],
                              '-', color=self.colour_list[counter], alpha=self.alpha_transparency_level)
        # add horizontal line at y50
        self.axarr[1, 0].hlines(y=temp_data["y50_normalized"]["value"], xmin=xmin_norm, colors=self.colour_list[counter],
                                xmax=temp_data["EC50_norm_bq"]["value"], linestyles=self.linestyles[sample])

        # add vertical line at EC50
        self.axarr[1, 0].vlines(x=temp_data["EC50_norm_bq"]["value"], ymin=ymin_norm, colors=self.colour_list[counter],
                                ymax=temp_data["y50_normalized"]["value"], linestyles=self.linestyles[sample])

        # set xlabel, ylabel, title, grid, etc
        self.axarr[1, 0].set_xlabel("dose concentration (normalised)",
                                    fontsize=self.fig_fontsize)
        self.axarr[1, 0].set_ylabel("response concentration (normalised)", rotation='vertical',
                                    fontsize=self.fig_fontsize)
        self.axarr[1, 0].text(0.6, 1.1, "normalised data", horizontalalignment='center',
                              fontsize=self.fig_fontsize)
        self.axarr[1, 0].grid(True, color='0.75')
        self.axarr[1, 0].set_ylim(ymin_norm, 1.2)
        self.axarr[1, 0].set_xlim(xmin_norm, 1.2)

    def plot_automatic_judgement_fig3(self, counter, sample, temp_data, annotation_location, yaxis_pos, xaxis_left):

        title_summ = sample
        self.axarr[0, 1].set_title(title_summ, fontsize=self.annotation_fontsize, color="k", alpha=0.75)

        self.axarr[0, 1].annotate(text=f"{self.method} ({self.doseunits})", xy=(xaxis_left, yaxis_pos[3]),
                                  fontsize=self.annotation_fontsize, xycoords=self.xy_coordinates, alpha=0.75)
        self.axarr[0, 1].annotate(text="rsquared", xy=(xaxis_left, yaxis_pos[4]),
                                  fontsize=self.annotation_fontsize, xycoords=self.xy_coordinates, alpha=0.75)
        self.axarr[0, 1].annotate(text="hillslope", xy=(xaxis_left, yaxis_pos[5]),
                                  fontsize=self.annotation_fontsize, xycoords=self.xy_coordinates, alpha=0.75)
        self.axarr[0, 1].annotate(text="n_lowdose_datapoints", xy=(xaxis_left, yaxis_pos[6]),
                                  fontsize=self.annotation_fontsize, xycoords=self.xy_coordinates, alpha=0.75)
        self.axarr[0, 1].annotate(text="std_lowdose_datapoints", xy=(xaxis_left, yaxis_pos[7]),
                                  fontsize=self.annotation_fontsize, xycoords=self.xy_coordinates, alpha=0.75)
        self.axarr[0, 1].annotate(text="n_highdose_datapoints", xy=(xaxis_left, yaxis_pos[8]),
                                  fontsize=self.annotation_fontsize, xycoords=self.xy_coordinates, alpha=0.75)
        self.axarr[0, 1].annotate(text="std_highdose_datapoints", xy=(xaxis_left, yaxis_pos[9]),
                                  fontsize=self.annotation_fontsize, xycoords=self.xy_coordinates, alpha=0.75)

        # add headers to table showing the rsquared and other aspects of the fit and dataset
        self.axarr[0, 1].annotate(text=sample, xy=(annotation_location, yaxis_pos[2]),
                                  fontsize=self.annotation_fontsize, xycoords=self.xy_coordinates, alpha=0.75)
        if temp_data["EC50_calculable"]:
            EC50colour = "k" if temp_data["approved"] == True else "r"
            self.axarr[0, 1].annotate(text=round(temp_data["EC50"]["value"], 2),
                                      xy=(annotation_location, yaxis_pos[3]), fontsize=self.annotation_fontsize,
                                      xycoords=self.xy_coordinates, alpha=0.75, color=EC50colour)
            # rsquared of the fit to the data
            self.axarr[0, 1].annotate(text=round(temp_data["rsquared"]["value"], 2),
                                      xy=(annotation_location, yaxis_pos[4]), fontsize=self.annotation_fontsize,
                                      xycoords=self.xy_coordinates, alpha=0.75, color=temp_data["rsquared"]["colour"])
            # hillslope of the fit to the data
            self.axarr[0, 1].annotate(text=round(temp_data["hillslope"]["value"], 2),
                                      xy=(annotation_location, yaxis_pos[5]), fontsize=self.annotation_fontsize,
                                      xycoords=self.xy_coordinates, alpha=0.75, color=temp_data["hillslope"]["colour"])
            # number of lowdose datapoints
            self.axarr[0, 1].annotate(text=round(temp_data["n_lowdose_datapoints"]["value"], 2),
                                      xy=(annotation_location, yaxis_pos[6]), fontsize=self.annotation_fontsize,
                                      xycoords=self.xy_coordinates, alpha=0.75,
                                      color=temp_data["n_lowdose_datapoints"]["colour"])
            # std of lowdose datapoints
            self.axarr[0, 1].annotate(text=round(temp_data["std_resp_lowdose_datapoints"]["value"], 2),
                                      xy=(annotation_location, yaxis_pos[7]), fontsize=self.annotation_fontsize,
                                      xycoords=self.xy_coordinates, alpha=0.75,
                                      color=temp_data["std_resp_lowdose_datapoints"]["colour"])
            # number of highdose datapoints
            self.axarr[0, 1].annotate(text=round(temp_data["n_highdose_datapoints"]["value"], 2),
                                      xy=(annotation_location, yaxis_pos[8]), fontsize=self.annotation_fontsize,
                                      xycoords=self.xy_coordinates, alpha=0.75,
                                      color=temp_data["n_highdose_datapoints"]["colour"])
            # std of highdose datapoints
            self.axarr[0, 1].annotate(text=round(temp_data["std_resp_highdose_datapoints"]["value"], 2),
                                      xy=(annotation_location, yaxis_pos[9]), fontsize=self.annotation_fontsize,
                                      xycoords=self.xy_coordinates, alpha=0.75,
                                      color=temp_data["std_resp_highdose_datapoints"]["colour"])

        else:
            if isinstance(temp_data["EC50"], str):
                temp_data["EC50_to_insert"] = "N/A"
            elif isinstance(temp_data["EC50"], float):
                temp_data["EC50_to_insert"] = temp_data["EC50"]
            # insert error string or EC50, coloured red to indicate likely poor data
            self.axarr[0, 1].annotate(text=temp_data["EC50_to_insert"], xy=(annotation_location, yaxis_pos[3]),
                                      fontsize=self.annotation_fontsize, xycoords=self.xy_coordinates, alpha=0.75,
                                      color="r")

    def plot_automatic_judgement_fig4(self, counter, sample, temp_data, annotation_location, yaxis_pos, xaxis_left):

        self.axarr[1, 1].annotate(text="dose conc. stepsize", xy=(xaxis_left, yaxis_pos[4]),
                                  fontsize=self.annotation_fontsize, xycoords=self.xy_coordinates, alpha=0.75)
        self.axarr[1, 1].annotate(text="slope at lowdose", xy=(xaxis_left, yaxis_pos[5]),
                                  fontsize=self.annotation_fontsize, xycoords=self.xy_coordinates, alpha=0.75)
        self.axarr[1, 1].annotate(text="slope at highdose", xy=(xaxis_left, yaxis_pos[6]),
                                  fontsize=self.annotation_fontsize, xycoords=self.xy_coordinates, alpha=0.75)
        self.axarr[1, 1].annotate(text=sample, xy=(annotation_location, yaxis_pos[2]),
                                  fontsize=self.annotation_fontsize, xycoords=self.xy_coordinates, alpha=0.75)

        if temp_data["EC50_calculable"]:
            # add the stepsize near the EC50, which determines whether more dose concentrations are necessary
            stepsize = temp_data["doseconc_stepsize_at_EC50"]["value"]
            stepcolour = temp_data["doseconc_stepsize_at_EC50"]["colour"]

            # if dfe.loc["doseconc_stepsize_at_EC50{}".format(d),"%s_okay" % sLet] == True:
            self.axarr[1, 1].annotate(text=round(stepsize, 2), xy=(annotation_location, yaxis_pos[4]),
                                      fontsize=self.annotation_fontsize, xycoords=self.xy_coordinates, alpha=0.75,
                                      color=stepcolour)

            self.axarr[1, 1].annotate(text=round(temp_data["saxe_lowdose"]["value"], 2),
                                      xy=(annotation_location, yaxis_pos[5]), fontsize=self.annotation_fontsize,
                                      xycoords=self.xy_coordinates, alpha=0.75,
                                      color=temp_data["saxe_lowdose"]["colour"])
            self.axarr[1, 1].annotate(text=round(temp_data["saxe_highdose"]["value"], 2),
                                      xy=(annotation_location, yaxis_pos[6]), fontsize=self.annotation_fontsize,
                                      xycoords=self.xy_coordinates, alpha=0.75,
                                      color=temp_data["saxe_highdose"]["colour"])
            if stepcolour != "k":

                self.axarr[1, 0].plot(temp_data["doseconc_steps_at_EC50"], (0, 0), color=stepcolour, linestyle="-",
                                      lw=2)

            # test_max_low_dose_slop = settings["max_lowdose_slope"]
            # test_max_high_dose_slop = settings["max_highdose_slope"]

            test_max_low_dose_slop = 1
            test_max_high_dose_slop = 1


            # # if the slope at the lowdose is above the chosen cutoff, draw a line on the normalised plot, axarr[1,0]
            # if temp_data["saxe_lowdose"]["value"] > test_max_low_dose_slop:
            #
            #     saxe_lowdose_values = temp_data["saxe_lowdose_values"]["value"]
            #     # draw red vertical line showing the slope at the lowdose datapoint
            #     self.axarr[1, 0].plot(saxe_lowdose_values[0], saxe_lowdose_values[1], 'r-', lw=2)
            # # if the slope at the highdose is higher than the chosen cutoff, draw a line on tho normalised plot
            # if temp_data["saxe_highdose"]["value"] > test_max_high_dose_slop:
            #     saxe_highdose_values = temp_data["saxe_highdose_values"]["value"]
            #     # draw red vertical line showing the slope at the lowdose datapoint
            #     self.axarr[1, 0].plot(saxe_highdose_values[0], saxe_highdose_values[1], 'r-', lw=2)

        else:
            # optional: print the dataframe showing which parameters are not acceptable
            if isinstance(temp_data["EC50"]["value"], str):
                # define x_position of annotation. Place the orig at around 0.6, & 2nd dataset(_ful) at around 0.3
                xd_wide = annotation_location-0.3-0.3 * counter
                self.axarr[1, 1].annotate(text=sample, xy=(xd_wide, yaxis_pos[5]),
                                          fontsize=self.annotation_fontsize, xycoords=self.xy_coordinates, color="r")
                self.axarr[1, 1].annotate(text=temp_data["EC50"]["value"], xy=(xd_wide, yaxis_pos[6]),
                                          fontsize=self.annotation_fontsize, xycoords=self.xy_coordinates, color="r")
    @staticmethod
    def _anno_location(counter, start=0.85, dist=0.2):
        annotation_location = start - dist * counter
        return annotation_location

    def controller(self, save_location):

        yaxis_pos = np.linspace(0.9, 0.1, 10)
        xaxis_left = 0.05

        for counter, sample in enumerate(self.all_data):
            if sample != "state_data":
                self.plot_raw_fig1(counter, sample, self.all_data[sample])
                self.plot_fitted_fig1(counter, sample, self.all_data[sample])
                self.plot_normalised_fig2(counter, sample, self.all_data[sample])
                # analyse the curve fit and data to judge whether the EC50 value is accurate
                self.jf.judge_controller(self.all_data[sample])
                annotation_location = self._anno_location(counter)
                self.plot_automatic_judgement_fig3(counter, sample, self.all_data[sample], annotation_location,
                                                   yaxis_pos, xaxis_left)
                self.plot_automatic_judgement_fig4(counter, sample, self.all_data[sample], annotation_location,
                                                   yaxis_pos, xaxis_left)

                # fig_save = save_location/sample.png
                # # save figure with the fitted curve and calculated EC50 value
                # self.fig.tight_layout()
                # self.fig.savefig(fig_save, format='png', dpi=140)
                # if settings["save_as_pdf"] in (True,"TRUE"):
                #     fig.savefig(fig0_single_sample_pdf, format='pdf')
                # plt.close('all')
        # plt.legend()
        # plt.show()
        return self.all_data


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    plt.subplots()
