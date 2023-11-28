from PySimpleGUI import LISTBOX_SELECT_MODE_SINGLE, LISTBOX_SELECT_MODE_MULTIPLE, PopupError, PopupGetText
from matplotlib import pyplot as plt

from lcms_functions import get_peak_information, purity_plotting, add_start_end_time, purity_ops, grab_sample_data
from upstarts_values import window_1_lcms, all_table_data


def sample_selection_mode_update(window, values):
    if values["-PURITY_INFO_SAMPLE_SELECTION-"]:
        window["-PURITY_INFO_SAMPLE_BOX-"].update(select_mode=LISTBOX_SELECT_MODE_MULTIPLE)
    else:
        window["-PURITY_INFO_SAMPLE_BOX-"].update(select_mode=LISTBOX_SELECT_MODE_SINGLE)


def lcms_calculation(config, window, values):
    if values["-PURITY_INFO_SAMPLE_BOX-"]:
        samples = values["-PURITY_INFO_SAMPLE_BOX-"]
        slope_threshold = int(values["-PURITY_INFO_SLOPE_THRESHOLD-"])
        uv_threshold = int(values["-PURITY_INFO_UV_THRESHOLD-"])
        rt_solvent_peak = float(values["-PURITY_INFO_RT_SOLVENT-"])
        wavelength_data = values["-PURITY_INFO_UV_WAVE-"]
        peak_information = {}
        peak_table_data = {}
        sample_peak_dict = {}

        for sample in samples:
            temp_peak_information, temp_peak_table_data, temp_sample_peak_dict = get_peak_information(
                window_1_lcms["purity_data"], slope_threshold, uv_threshold, rt_solvent_peak, None, wavelength_data, sample)

            peak_information[sample] = temp_peak_information[sample]
            peak_table_data[sample] = temp_peak_table_data[sample]
            sample_peak_dict[sample] = temp_sample_peak_dict[sample]

        if values["-PURITY_INFO_USE_MS_DATA-"]:
            if not window_1_lcms["sample_data_file"]:
                mass = PopupGetText("What is the mass of the compound?")
            else:
                mass = None

            if values["-PURITY_INFO_MS_MODE_POS-"]:
                ms_mode = "ms_pos"
            elif values["-PURITY_INFO_MS_MODE_NEG-"]:
                ms_mode = "ms_neg"
            else:
                ms_mode = "both"    # ToDo Make this work ?

            delta_mass = float(values["-PURITY_INFO_MS_DELTA-"])
            mz_threshold = int(values["-PURITY_INFO_MS_THRESHOLD-"])
            peak_amounts = int(values["-PURITY_INFO_MS_PEAKS-"])

            sample_data, _ = grab_sample_data(window_1_lcms["sample_data_file"], window_1_lcms["purity_data"], config)

            all_table_data["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"], \
            purity_peak_list_table_data = purity_ops(sample_data, window_1_lcms["purity_data"], peak_information, ms_mode,
                                                     delta_mass, mz_threshold, peak_amounts, mass)

            add_start_end_time(purity_peak_list_table_data, sample_peak_dict)
            window["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"]. \
                update(values=all_table_data["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"])
            return purity_peak_list_table_data, peak_table_data
    else:
        PopupError("Missing sample Information - Please select data")


def lcms_drawing(check, window, values, event, peak_table_data, lc_graph_showing, purity_peak_list_table_data):
    if check:
        if event == "-PURITY_INFO_SAMPLE_BOX-":
            window_1_lcms["purity_info_samples"] = values["-PURITY_INFO_SAMPLE_BOX-"]
            lc_method = values["-PURITY_INFO_GRAPH_SHOWING-"]

            if lc_method == lc_graph_showing[3]:
                start = values["-PURITY_INFO_RT_START-"]
                end = values["-PURITY_INFO_RT_END-"]
                peak = "None"
                window_1_lcms["canvas_lines"]["peak_lines"][peak] = {"start": start, "end": end}

        elif event == "-PURITY_INFO_GRAPH_SHOWING-":
            lc_method = values["-PURITY_INFO_GRAPH_SHOWING-"]

        elif event == "-PURITY_INFO_PURITY_OVERVIEW_TABLE-" and values["-PURITY_INFO_GRAPH_SHOWING-"]:
            lc_method = lc_graph_showing[0]
            window["-PURITY_INFO_GRAPH_SHOWING-"].update(value=lc_method)
            print(f'All table data - purity info: {all_table_data["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"]}')
            window["-PURITY_INFO_MZ-"].update(value=all_table_data["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"][values[
                "-PURITY_INFO_PURITY_OVERVIEW_TABLE-"][0]][1])
            window_1_lcms["purity_info_samples"] = [all_table_data["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"][values[
                "-PURITY_INFO_PURITY_OVERVIEW_TABLE-"][0]][0]]
            all_table_data["-PURITY_INFO_PURITY_PEAK_LIST_TABLE-"] = \
                purity_peak_list_table_data[window_1_lcms["purity_info_samples"][0]]
            window["-PURITY_INFO_PURITY_PEAK_LIST_TABLE-"].update(
                values=all_table_data["-PURITY_INFO_PURITY_PEAK_LIST_TABLE-"])
            window["-PURITY_INFO_PEAK_LIST_SAMPLE_TEXT-"].update(value=window_1_lcms["purity_info_samples"][0])
            window["-PURITY_INFO_PEAK_LIST_SAMPLE-"].update(value=window_1_lcms["purity_info_samples"])

        elif event == "-PURITY_INFO_PURITY_PEAK_LIST_TABLE-" and values["-PURITY_INFO_PURITY_PEAK_LIST_TABLE-"]:
            lc_method = lc_graph_showing[3]
            window["-PURITY_INFO_GRAPH_SHOWING-"].update(value=lc_method)
            window_1_lcms["purity_info_samples"] = values["-PURITY_INFO_PEAK_LIST_SAMPLE-"]
            window_1_lcms["purity_info_rt_start"] = all_table_data["-PURITY_INFO_PURITY_PEAK_LIST_TABLE-"][
                values["-PURITY_INFO_PURITY_PEAK_LIST_TABLE-"][0]][4]
            window["-PURITY_INFO_RT_START-"].update(value=window_1_lcms["purity_info_rt_start"])

            window_1_lcms["purity_info_rt_end"] = all_table_data["-PURITY_INFO_PURITY_PEAK_LIST_TABLE-"][
                values["-PURITY_INFO_PURITY_PEAK_LIST_TABLE-"][0]][5]
            window["-PURITY_INFO_RT_END-"].update(value=window_1_lcms["purity_info_rt_end"])

            window_1_lcms["purity_info_samples"] = window_1_lcms["purity_info_samples"].strip("',)")
            window_1_lcms["purity_info_samples"] = window_1_lcms["purity_info_samples"].split("'")[1]
            window_1_lcms["purity_info_samples"] = [window_1_lcms["purity_info_samples"]]
            window_1_lcms["canvas_lines"]["peak_lines"][peak] = {"start": window_1_lcms["purity_info_rt_start"],
                                                                 "end": window_1_lcms["purity_info_rt_end"]}
        # Set size of the canvas figure

        elif event == "-PURITY_INFO_DRAW_PEAKS-":
            for data in all_table_data["-PURITY_INFO_PEAK_TABLE-"]:
                peak = data[1]
                start = data[3]
                end = data[4]
                window_1_lcms["canvas_lines"]["peak_lines"][peak] = {"start": start, "end": end}
                window_1_lcms["update_purity_info_peak_table"] = False

        elif event == "-PURITY_INFO_PEAK_TABLE-" and values["-PURITY_INFO_PEAK_TABLE-"]:
            purity_info_samples = [
                all_table_data["-PURITY_INFO_PEAK_TABLE-"][values["-PURITY_INFO_PEAK_TABLE-"][0]][0]]
            peak = all_table_data["-PURITY_INFO_PEAK_TABLE-"][values["-PURITY_INFO_PEAK_TABLE-"][0]][1]
            start = all_table_data["-PURITY_INFO_PEAK_TABLE-"][values["-PURITY_INFO_PEAK_TABLE-"][0]][3]
            end = all_table_data["-PURITY_INFO_PEAK_TABLE-"][values["-PURITY_INFO_PEAK_TABLE-"][0]][4]
            window_1_lcms["canvas_lines"]["peak_lines"][peak] = {"start": start, "end": end}

            if values["-PURITY_INFO_RADIO_PEAKS_UV-"]:
                lc_method = lc_graph_showing[0]
            elif values["-PURITY_INFO_RADIO_PEAKS_MS_SPECTRA-"]:
                lc_method = lc_graph_showing[3]

            window["-PURITY_INFO_GRAPH_SHOWING-"].update(value=lc_method)
            window["-PURITY_INFO_RT_START-"].update(value=start)
            window["-PURITY_INFO_RT_END-"].update(value=end)
            window_1_lcms["update_purity_info_peak_table"] = False

        if values["-PURITY_INFO_DRAW_THRESHOLD-"]:
            uv_line = float(values["-PURITY_INFO_UV_THRESHOLD-"])
            window_1_lcms["canvas_lines"]["uv"] = uv_line

        fig_size = (7, 3)
        print("MISSING GUARD TO PREVENT DIFFERENT LENGTH DATA CRASHING THE PROGRAM!!!")

        purity_info_canvas = window["-PURITY_INFO_CANVAS-"]

        if values["-PURITY_INFO_MS_MODE_POS-"]:
            ms_mode = "ms_pos"
        elif values["-PURITY_INFO_MS_MODE_NEG-"]:
            ms_mode = "ms_neg"
        else:
            ms_mode = "both"        # ToDo make this work

        if not window_1_lcms["purity_info_rt_start"]:
            try:
                window_1_lcms["purity_info_rt_start"] = float(values["-PURITY_INFO_RT_START-"])
            except ValueError:
                window_1_lcms["purity_info_rt_start"] = None
        try:
            wavelength = float(values["-PURITY_INFO_WAVELENGTH-"])
        except ValueError:
            wavelength = None
        try:
            bin_numbers = int(values["-PURITY_INFO_BIN-"])
        except ValueError:
            bin_numbers = None
        try:
            mz_value = float(values["-PURITY_INFO_MZ-"])
        except ValueError:
            mz_value = None

        if window_1_lcms["purity_info_samples"]:
            plot_style = purity_plotting(lc_method, window_1_lcms["purity_data"], purity_info_canvas,
                                         window_1_lcms["purity_info_samples"], fig_size, ms_mode,
                                         window_1_lcms["purity_info_rt_start"], window_1_lcms["purity_info_rt_end"],
                                         wavelength, bin_numbers, mz_value, window_1_lcms["canvas_lines"])

        else:
            fig = plt.gcf()
            plt.close(fig)
            try:
                window_1_lcms["temp_purity_info_canvas"].get_tk_widget().forget()
            except AttributeError:
                print("Attribute Error for info canvas")

        if type(plot_style) == str:
            PopupError(plot_style)
            plot_style = None
        else:

            if window_1_lcms["temp_purity_info_canvas"] is not None:
                if not window_1_lcms["temp_purity_info_canvas"] == plot_style:
                    fig = plt.gcf()
                    plt.close(fig)
                    window_1_lcms["temp_purity_info_canvas"].get_tk_widget().forget()
                else:
                    plot_style.get_tk_widget().forget()
            if purity_info_samples:
                plot_style.draw()

                if window_1_lcms["toolbar"]:
                    window_1_lcms["toolbar"].destroy()

                toolbar = window_1_lcms["Toolbar"](plot_style, window["-PURITY_INFO_CANVAS_TOOLBAR-"].TKCanvas)
                toolbar.update()
                plot_style.get_tk_widget().pack()
            try:
                window.refresh()
            except AttributeError:
                print("Canvas - AttributeError on window.refresh")
            window_1_lcms["temp_purity_info_canvas"] = plot_style
            temp_peak_table_data = []
            if len(purity_info_samples) > 1:
                for sample in purity_info_samples:
                    for rows in peak_table_data[sample]:
                        temp_peak_table_data.append(rows)
            else:
                try:
                    temp_peak_table_data = peak_table_data[purity_info_samples[0]]
                except IndexError:
                    temp_peak_table_data = ""

            if window_1_lcms["update_purity_info_peak_table"]:
                all_table_data["-PURITY_INFO_PEAK_TABLE-"] = temp_peak_table_data
                window["-PURITY_INFO_PEAK_TABLE-"].update(values=all_table_data["-PURITY_INFO_PEAK_TABLE-"])
            all_table_data["-PURITY_INFO_RAW_DATA_TABLE-"] = \
                window_1_lcms["purity_data"][purity_info_samples[0]]["peak_table_raw"][1]
            window["-PURITY_INFO_RAW_DATA_TABLE-"].update(values=all_table_data["-PURITY_INFO_RAW_DATA_TABLE-"])

            window_1_lcms["update_purity_info_peak_table"] = True
            window_1_lcms["purity_info_rt_start"] = None
            window_1_lcms["purity_info_rt_end"] = None
            window_1_lcms["canvas_lines"] = {
                "uv": None,
                "peak_lines": {}
            }
    else:
        PopupError(f"Missing {check}")
