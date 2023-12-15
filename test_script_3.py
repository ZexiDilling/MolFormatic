import configparser

from database_controller import FetchData

from lcms_functions import _compound_list, table_update_tree

from lcms_functions import get_peak_information, lcms_plotting, add_start_end_time, lcms_ops, grab_sample_data
from start_up_values import window_1_lcms

def testing(check, window, values, event, peak_table_data, lc_graph_showing, purity_peak_list_table_data):

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


if __name__ == "__main__":
    config = configparser.ConfigParser()
    config.read("config.ini")

    mp_amount = None
    min_mp = None
    samples_per_plate = None
    ignore_active = True
    sub_search = True
    smiles = "Cc1=c-c=c(N2CCC(C)=N2)-c=c-1"
    sub_search_methode = "finger"
    threshold = 40
    source_table = config["Tables"]["compound_main"]
    fd = FetchData(config)
    search_limiter = {
                config["Tables"]["compound_source"]: {"academic_commercial": {"value": None,
                                                                              "operator": "IN",
                                                                              "target_column": "ac",
                                                                              "use": False},
                                                      "vendor_center": {"value": None,
                                                                        "operator": "IN",
                                                                        "target_column": "origin",
                                                                        "use": False}},
                config["Tables"]["compound_main"]: {"origin_id": {"value": "",
                                                                  "operator": "IN",
                                                                  "target_column": "ac_id",
                                                                  "use": False},
                                                    "volume": {"value": None,
                                                               "operator": "<",
                                                               "target_column": "volume",
                                                               "use": False}},
                "join_tables": {config["Tables"]["compound_main"]: {},
                                config["Tables"]["compound_mp_table"]: {
                                    "compound_id": {"value": "",
                                                    "operator": "IN",
                                                    "target_column": "compound_id",
                                                    "use": False}},
                                "shared_data": "compound_id"}
            }

    treedata, all_data, rows, counter = table_update_tree(mp_amount, min_mp, samples_per_plate, ignore_active,
                                                                   sub_search, smiles, sub_search_methode,
                                                                   threshold, source_table, search_limiter, config)


