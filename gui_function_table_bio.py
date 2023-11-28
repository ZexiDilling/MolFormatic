from PySimpleGUI import PopupError

from database_handler import DataBaseFunctions
from gui_function_info_compound import update_overview_compound
from database_functions import grab_table_data

from start_up_values import window_tables


def table_tab_group_pressed_update(config, window, values):
    if values["-TABLE_TAB_GRP-"] == "Bio Experiment table":

        if window_tables["all_assays"] is None:
            all_assays_data, _ = grab_table_data(config, "assay")
            all_assays = []
            for rows in all_assays_data:
                all_assays.append(rows[1])

            window["-BIO_EXP_TABLE_ASSAY_LIST_BOX-"].update(values=all_assays)


def experiment_table_assay_list_update(config, window, values):

    if values["-BIO_EXP_TABLE_ASSAY_LIST_BOX-"]:
        # clearing all tables
        window["-BIO_EXP_PLATE_TABLE-"].update(values=[[]])
        window["-BIO_EXP_COMPOUND_TABLE-"].update(values=[[]])
        window["-BIO_EXP_PLATE_NOTE-"].update(value="")

        # updating run overview table,
        table_name = "assay_runs"
        selected_assays = values["-BIO_EXP_TABLE_ASSAY_LIST_BOX-"]
        selected_headlines = ["run_name", "batch", "date", "note"]
        search_list_clm = "assay_name"
        bio_exp_assay_runs, _ = grab_table_data(config, table_name, selected_assays,
                                                specific_rows=selected_headlines, search_list_clm=search_list_clm)

        window["-BIO_EXP_ASSAY_RUN_TABLE-"].update(values=bio_exp_assay_runs)
        return bio_exp_assay_runs


def experiment_table_assay_run_update(config, window, values, bio_exp_assay_runs):
    if values["-BIO_EXP_ASSAY_RUN_TABLE-"]:

        # clearing experiment and plate table:
        window["-BIO_EXP_COMPOUND_TABLE-"].update(values=[[]])
        window["-BIO_EXP_PLATE_NOTE-"].update(value="")

        approval_check = values["-BIO_EXP_APPROVED_PLATES_ONLY-"]

        # updating bio plate table:
        bio_exp_plate_data = _update_bio_exp_plate_table(config, window, values, bio_exp_assay_runs,
                                                         approval_check=approval_check)


def experiment_table_plate_update(config, window, values, bio_exp_plate_data):
    if values["-BIO_EXP_PLATE_TABLE-"]:
        approval_check = values["-BIO_EXP_APPROVED_COMPOUNDS_ONLY-"]
        try:
            float(values["-BIO_EXP_SET_THRESHOLD-"])
        except ValueError:
            threshold = False
        else:
            threshold = float(values["-BIO_EXP_SET_THRESHOLD-"])

        try:
            float(values["-BIO_EXP_SET_COMPOUND_AMOUNT-"])
        except ValueError:
            compound_amount = False
        else:
            compound_amount = float(values["-BIO_EXP_SET_COMPOUND_AMOUNT-"])

        bio_exp_compound_data = _update_bio_exp_compound_table(config,  window, values, bio_exp_plate_data,
                                                               threshold=threshold, compound_amount=compound_amount,
                                                               approval_check=approval_check)
        return bio_exp_compound_data


def assay_table_double_click(window, values):
    print("we are double clicking - runs")


def plate_table_double_click(window, values):
    print("we are double clicking - plates")


def compound_table_double_click(config, window, values, bio_exp_compound_data):
    try:
        table_row = values["-BIO_EXP_COMPOUND_TABLE-"][0]
    except IndexError:
        pass
    else:
        print(bio_exp_compound_data[table_row])
        temp_compound_id = bio_exp_compound_data[table_row][0]
        print(temp_compound_id)
        window["-COMPOUND_INFO_ID-"].update(value=temp_compound_id)
        update_overview_compound(config, window, values, temp_compound_id)


def _update_bio_exp_plate_table(config, window, values, bio_exp_assay_runs, approval_check=False):
    table_name = "biological_plate_data"
    bio_exp_selected_runs = []
    for temp_values in values["-BIO_EXP_ASSAY_RUN_TABLE-"]:
        bio_exp_selected_runs.append(bio_exp_assay_runs[temp_values][0])
    search_list_clm = "assay_run"
    selected_headlines = ["plate_name", "z_prime", "approval", "note", "analysed_method", "assay_run"]
    bio_exp_plate_data, _ = grab_table_data(config, table_name, bio_exp_selected_runs,
                                            specific_rows=selected_headlines, search_list_clm=search_list_clm)

    table_data = []
    if bio_exp_plate_data:
        for rows, row_data in enumerate(bio_exp_plate_data):
            if approval_check and row_data[2] == "True" or not approval_check:

                if type(row_data[1]) == float:
                    row_data[1] = round(bio_exp_plate_data[rows][1], 2)

                if row_data[3] == "Dead Plate - See Run note":
                    row_data[3] = "Dead"

                table_data.append(row_data)

        window["-BIO_EXP_PLATE_TABLE-"].update(values=table_data)
        # Update the plate counter:
        window["-BIO_EXP_PLATE_COUNTER-"].update(value=len(table_data))
        temp_plate = bio_exp_selected_runs[0]
        dbf = DataBaseFunctions(config)
        temp_data = dbf.find_data_single_lookup("assay_runs", temp_plate, "run_name")
        note = temp_data[0][7]
        window["-BIO_EXP_RUN_NOTE-"].update(value=note)
    else:
        table_data = None
        window["-BIO_EXP_PLATE_TABLE-"].update(values=[[]])
        window["-BIO_EXP_PLATE_COUNTER-"].update(value="0")
        bio_exp_plate_data.update(value="No Data")

    return table_data


def _update_bio_exp_compound_table(config, window, values, bio_exp_plate_data,
                                   threshold=False, compound_amount=False, approval_check=False, hit=False):

    bio_exp_selected_plates = []
    table_name = "biological_compound_data"
    for temp_values in values["-BIO_EXP_PLATE_TABLE-"]:
        bio_exp_selected_plates.append(bio_exp_plate_data[temp_values][0])

    search_list_clm = "assay_plate"
    selected_headlines = ["compound_id", "score", "hit", "concentration", "approved", "assay_well", "note"]

    bio_exp_compound_data, _ = grab_table_data(config, table_name, bio_exp_selected_plates,
                                               specific_rows=selected_headlines, search_list_clm=search_list_clm)
    # Sort the list based on the score
    if bio_exp_compound_data:
        bio_exp_compound_data.sort(key=lambda x: x[1])
        table_data = []
        if bio_exp_compound_data:
            for rows, row_data in enumerate(bio_exp_compound_data):
                if threshold:
                    try:
                        row_data[1] > threshold
                    except TypeError:
                        break
                    else:
                        if row_data[1] > threshold:
                            break

                if approval_check and row_data[4] == "1" or not approval_check:
                    if row_data[6] == "No transfer data":
                        row_data[6] = "Nah"
                    row_data.remove(row_data[4])

                    if type(row_data[1]) == float:
                        row_data[1] = round(row_data[1], 2)

                    table_data.append(row_data)
                if compound_amount and rows + 1 >= compound_amount:
                    break

        else:
            print(bio_exp_compound_data)

        # Update compound table data and counter:
        window["-BIO_EXP_COMPOUND_TABLE-"].update(values=table_data)
        window["-BIO_EXP_COMPOUND_COUNTER-"].update(value=len(table_data))
    else:
        table_data = None

    # update the note field for the plate table:
    # try:
    bio_exp_plate_data[values["-BIO_EXP_PLATE_TABLE-"][0]]
    print("Testing error !!!! if it fails")
    # except:
    #
    #     print(values["-BIO_EXP_PLATE_TABLE-"])
    #     note = ""
    # else:
    temp_plate = bio_exp_plate_data[values["-BIO_EXP_PLATE_TABLE-"][0]][0]
    dbf = DataBaseFunctions(config)
    temp_data = dbf.find_data_single_lookup("biological_plate_data", temp_plate, "plate_name")
    note = temp_data[0][9]
    window["-BIO_EXP_PLATE_NOTE-"].update(value=note)

    return table_data

def bio_exp_compound_list(config, event, values, bio_exp_compound_data):

    # Sort the list based on the score
    bio_exp_compound_data.sort(key=lambda x: x[1])

    try:
        float(values["-BIO_EXP_SET_THRESHOLD-"])
    except ValueError:
        threshold = False
    else:
        threshold = float(values["-BIO_EXP_SET_THRESHOLD-"])

    try:
        float(values["-BIO_EXP_SET_COMPOUND_AMOUNT-"])
    except ValueError:
        compound_amount = False
    else:
        compound_amount = float(values["-BIO_EXP_SET_COMPOUND_AMOUNT-"])

    bio_list = []
    for rows, row_data in enumerate(bio_exp_compound_data):
        if threshold and row_data[1] > threshold:
            break

        bio_list.append(row_data[0])
        if compound_amount and rows + 1 >= compound_amount:
            break
    print("testing!!! ")

#ToDo Write information for the whole dose-response module!!!

