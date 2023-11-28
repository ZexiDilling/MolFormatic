from PySimpleGUI import Popup, PopupError, PopupGetText, TreeData

from lcms_functions import table_update_tree, compound_export, dp_creator, compound_info_table_data
from database_functions import grab_table_data
from info import vol_converter
from start_up_values import window_1_search, all_table_data


def tree_database_update(config, window, values, compound_data):
    try:
        temp_id = window.Element("-TREE_DB-").SelectedRows[0]
    except IndexError:
        pass
    # temp_info = window.Element("-TREE_DB-").TreeData.tree_dict[temp_id].values'
    tree_sample = compound_data[temp_id]["compound_id"]
    window["-COMPOUND_INFO_ID-"].update(value=compound_data[temp_id]["compound_id"])
    window["-COMPOUND_INFO_SMILES-"].update(value=compound_data[temp_id]["smiles"])
    window["-COMPOUND_INFO_MP_VOLUME-"].update(value=compound_data[temp_id]["volume"])
    window["-COMPOUND_INFO_PIC-"].update(data=compound_data[temp_id]["png"])
    window["-COMPOUND_INFO_ORIGIN_ID-"].update(value=compound_data[temp_id]["origin_id"])
    window["-COMPOUND_INFO_CONCENTRATION-"].update(value=compound_data[temp_id]["concentration"])

    search_limiter_origin = {"academic_commercial": {"value": values["-SEARCH_AC-"],
                                                     "operator": "IN",
                                                     "target_column": "ac",
                                                     "use": window_1_search["ac_use"]},
                             "vendor_center": {"value": values["-SEARCH_ORIGIN-"],
                                               "operator": "IN",
                                               "target_column": "origin",
                                               "use": window_1_search["origin_use"]}}
    all_data_origin, _ = grab_table_data(config, "origin", search_limiter_origin)
    window["-COMPOUND_INFO_AC-"].update(value=all_data_origin[0][1])
    window["-COMPOUND_INFO_ORIGIN-"].update(value=all_data_origin[0][2])

    compound_id = compound_data[temp_id]["compound_id"]

    # Table updates:
    all_table_data["-COMPOUND_INFO_ALL_PLATE_INFO_TABLE-"], all_table_data["-COMPOUND_INFO_MP_PLATE_INFO_TABLE-"], \
        all_table_data["-COMPOUND_INFO_DP_PLATE_INFO_TABLE-"], all_table_data["-COMPOUND_INFO_BIO_INFO_TABLE-"], \
        all_table_data["-COMPOUND_INFO_PURITY_INFO_TABLE-"] = compound_info_table_data(config, tree_sample)

    window["-COMPOUND_INFO_ALL_PLATE_INFO_TABLE-"].update(
        values=all_table_data["-COMPOUND_INFO_ALL_PLATE_INFO_TABLE-"])
    window["-COMPOUND_INFO_MP_PLATE_INFO_TABLE-"].update(
        values=all_table_data["-COMPOUND_INFO_MP_PLATE_INFO_TABLE-"])
    window["-COMPOUND_INFO_DP_PLATE_INFO_TABLE-"].update(
        values=all_table_data["-COMPOUND_INFO_DP_PLATE_INFO_TABLE-"])
    window["-COMPOUND_INFO_BIO_INFO_TABLE-"].update(
        values=all_table_data["-COMPOUND_INFO_BIO_INFO_TABLE-"])
    window["-COMPOUND_INFO_PURITY_INFO_TABLE-"].update(
        values=all_table_data["-COMPOUND_INFO_PURITY_INFO_TABLE-"])

    return compound_id


def compound_table_refreshed(config, window, values):
    print(window_1_search)
    if not window_1_search["compound_table_clear"]:
        if values["-SEARCH_ALL_COMPOUNDS-"]:
            values["-IGNORE_ACTIVE-"] = True
            values["-SUB_SEARCH-"] = False
            values["-SEARCH_PLATE_AMOUNT_MAX-"] = True
            values["-SEARCH_IGNORE_VOLUME-"] = True
            values["-SEARCH_AC-"] = None
            values["-SEARCH_ORIGIN-"] = None

        if values["-SEARCH_PLATE_PRODUCTION-"] == "Daughter Plates":
            table = "join_main_mp"

        elif values["-SEARCH_PLATE_PRODUCTION-"] == "Mother Plates":
            table = config["Tables"]["compound_main"]
        window_1_search["current_table_data"] = values["-SEARCH_PLATE_PRODUCTION-"]

        if values["-SEARCH_PLATE_AMOUNT-"] == "" and not values["-SEARCH_PLATE_AMOUNT_MAX-"]:
            PopupError("Please fill out plate amount")
            return None, None, None, None
        elif not values["-SEARCH_TRANS_VOL-"] and not values["-SEARCH_IGNORE_VOLUME-"]:
            PopupError("Please specify transferee amount")
            return None, None, None, None
        else:
            if not values["-SEARCH_PLATE_AMOUNT_MAX-"]:
                mp_amount = int(values["-SEARCH_PLATE_AMOUNT-"])
            else:
                mp_amount = None
            if not values["-SEARCH_IGNORE_VOLUME-"]:

                transferee_volume = float(values["-SEARCH_TRANS_VOL-"]) / vol_converter[
                    values["-SEARCH_VOL_PARAMETERS-"]]
            else:
                transferee_volume = None

            ignore_active = values["-SEARCH_IGNORE_PLATED_COMPOUNDS-"]
            sub_search = values["-SUB_SEARCH-"]
            smiles = values["-SUB_SEARCH_SMILES-"]
            sub_search_methode = values["-SUB_SEARCH_METHOD-"]
            threshold = float(values["-SUB_SEARCH_THRESHOLD-"])
            source_table = table

            if values["-SEARCH_AC-"]:
                ac_use = True
            else:
                ac_use = False

            if values["-SEARCH_ORIGIN-"]:
                origin_use = True
            else:
                origin_use = False

            samples_per_plate = int(values["-SEARCH_PLATE_LAYOUT_SAMPLE_AMOUNT-"])
            search_limiter = {
                config["Tables"]["compound_source"]: {"academic_commercial": {"value": values["-SEARCH_AC-"],
                                                                              "operator": "IN",
                                                                              "target_column": "ac",
                                                                              "use": ac_use},
                                                      "vendor_center": {"value": values["-SEARCH_ORIGIN-"],
                                                                        "operator": "IN",
                                                                        "target_column": "origin",
                                                                        "use": origin_use}},
                config["Tables"]["compound_main"]: {"origin_id": {"value": "",
                                                                  "operator": "IN",
                                                                  "target_column": "ac_id",
                                                                  "use": ac_use},
                                                    "volume": {"value": transferee_volume,
                                                               "operator": "<",
                                                               "target_column": "volume",
                                                               "use": not values["-SEARCH_IGNORE_VOLUME-"]}},
                "join_tables": {config["Tables"]["compound_main"]: {},
                                config["Tables"]["compound_mp_table"]: {
                                    "compound_id": {"value": "",
                                                    "operator": "IN",
                                                    "target_column": "compound_id",
                                                    "use": False}},
                                "shared_data": "compound_id"}
            }
            min_mp = values["-SEARCH_MP_MINIMIZED-"]
            table_data = table_update_tree(mp_amount, min_mp, samples_per_plate, ignore_active, sub_search,
                                           smiles,
                                           sub_search_methode, threshold, source_table, search_limiter, config)
            if table_data:
                treedata, all_data, compound_data, counter = table_data
                window['-TREE_DB-'].image_dict.clear()
                window["-TREE_DB-"].update(treedata)
                window["-C_TABLE_COUNT-"].update(f"Compounds: {counter}")
                window["-C_TABLE_REFRESH-"].update(text="Clear Table")
                window_1_search["compound_table_clear"] = True
                return treedata, all_data, compound_data, counter
            else:
                if values["-SEARCH_ALL_COMPOUNDS-"]:
                    PopupError("All compounds are in MotherPlates")
                    return None, None, None, None
                else:
                    PopupError("Database is empty")
                    return None, None, None, None
            #
            # except ValueError:
            #     sg.Popup("Fill in missing data")

    elif window_1_search["compound_table_clear"]:
        window['-TREE_DB-'].image_dict.clear()
        treedata = TreeData()
        window['-TREE_DB-'].update(treedata)
        window["-C_TABLE_REFRESH-"].update(text="Refresh")
        window["-C_TABLE_COUNT-"].update(f"Compounds: 0")
        window['-TREE_DB-'].image_dict.clear()
        window_1_search["compound_table_clear"] = False
        return None, None, None, None

    else:
        return None, None, None, None


def compound_table_export(config, window, values, all_data, archive_plates_dict):
    if not all_data:
        PopupError("Missing table data")
    elif not values["-SEARCH_OUTPUT_FOLDER-"]:
        PopupError("missing folder")
    else:
        if window_1_search["current_table_data"] == "Mother Plates":
            output_folder = values["-SEARCH_OUTPUT_FOLDER-"]
            all_compound_data = all_data["compound_list"]

            compound_export(output_folder, all_compound_data)
            file_location = f"{values['-SEARCH_OUTPUT_FOLDER-']}/comPOUND"

        elif window_1_search["current_table_data"] == "Daughter Plates":

            if not values["-SEARCH_PLATE_AMOUNT-"]:
                PopupError("Please fill out plate Amount")
            elif not values["-SEARCH_TRANS_VOL-"]:
                PopupError("Please fill out Transferee volume")
            elif not values["-SEARCH_PLATE_LAYOUT-"]:
                PopupError("Please select a plate Layout for the DP production")
            dp_name = PopupGetText("Dp name? ")
            if dp_name:
                plate_layout = archive_plates_dict[values["-SEARCH_PLATE_LAYOUT-"]]
                temp_sample_amount = int(values["-SEARCH_PLATE_LAYOUT_SAMPLE_AMOUNT-"])
                vol_converter = {"mL": 1000000, "uL": 10000, "nL": 1}

                transferee_volume = values["-SEARCH_TRANS_VOL-"] * vol_converter[
                    values["-SEARCH_VOL_PARAMETERS-"]]
                output_folder = values["-SEARCH_OUTPUT_FOLDER-"]
                mp_data = all_data["mp_data"]

                dp_creator(config, plate_layout, temp_sample_amount, mp_data, transferee_volume, dp_name,
                           output_folder)

                file_location = f"{values['-SEARCH_OUTPUT_FOLDER-']}/dp_output/"

        Popup(f"Done - files are located {file_location}")