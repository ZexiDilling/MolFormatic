from database_handler import DataBaseFunctions
from gui_function_info_bio import bio_info_plate_update, bio_info_window_update, \
    bio_info_plate_list_update
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


def experiment_table_assay_list_update(dbf, config, window, values):

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
        bio_exp_assay_runs, headlines = grab_table_data(config, table_name, selected_assays,
                                                specific_rows=selected_headlines, search_list_clm=search_list_clm)

        if values["-BIO_EXP_ANALYSE_STYLE-"] != "All":
            kill_counter = []
            run_name_counter = headlines.index("run_name")
            for row_counter, rows in enumerate(bio_exp_assay_runs):

                for counter, data in enumerate(rows):
                    if counter == run_name_counter:
                        row_data = dbf.find_data_single_lookup("biological_plate_data", data, "assay_run")[0]
                        if row_data[-1].casefold() != values["-BIO_EXP_ANALYSE_STYLE-"].casefold():
                            kill_counter.append(row_counter)

            if kill_counter:
                kill_counter.reverse()
                for killer_count in kill_counter:
                    del bio_exp_assay_runs[killer_count]

        window["-BIO_EXP_ASSAY_RUN_TABLE-"].update(values=bio_exp_assay_runs)


def experiment_table_assay_run_update(config, window, values):
    if values["-BIO_EXP_ASSAY_RUN_TABLE-"]:

        # clearing experiment and plate table:
        window["-BIO_EXP_COMPOUND_TABLE-"].update(values=[[]])
        window["-BIO_EXP_PLATE_NOTE-"].update(value="")

        approval_check = values["-BIO_EXP_APPROVED_PLATES_ONLY-"]

        # updating bio plate table:
        _update_bio_exp_plate_table(config, window, values, approval_check=approval_check)


def experiment_table_plate_update(config, window, values):
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

        _update_bio_exp_compound_table(config,  window, values, threshold=threshold, compound_amount=compound_amount,
                                       approval_check=approval_check)


def bio_tables_double_clicked(dbf, config, window, values, event, well_dict_bio_info):

    temp_plate = None
    temp_run = None
    temp_assay = None
    if event == "-BIO_EXP_PLATE_TABLE-+-double click-":
        try:
            table_row = values["-BIO_EXP_PLATE_TABLE-"][0]
        except IndexError:
            temp_indicator = False
        else:
            temp_assay = values["-BIO_EXP_TABLE_ASSAY_LIST_BOX-"]
            temp_table_data = window["-BIO_EXP_PLATE_TABLE-"].get()
            temp_run = temp_table_data[table_row][5]
            temp_plate = temp_table_data[table_row][0]
            temp_indicator = True

    elif event == "-BIO_EXP_ASSAY_RUN_TABLE-+-double click-":
        try:
            table_row = values["-BIO_EXP_ASSAY_RUN_TABLE-"][0]
        except IndexError:
            temp_indicator = False
        else:
            temp_table_data = window["-BIO_EXP_ASSAY_RUN_TABLE-"].get()
            temp_run = temp_table_data[table_row][0]
            temp_assay = values["-BIO_EXP_TABLE_ASSAY_LIST_BOX-"]
            temp_indicator = True

    elif event == "-BIO_EXP_TABLE_ASSAY_LIST_BOX-+-double click-":
        temp_assay = values["-BIO_EXP_TABLE_ASSAY_LIST_BOX-"]
        if temp_assay:
            temp_indicator = True
        else:
            temp_indicator = False

    else:
        temp_indicator = False

    if temp_indicator:
        try:
            temp_assay[0]
        except IndexError:
            pass
        else:
            temp_assay = temp_assay[0]

        if temp_plate:
            if not temp_run:
                row_data = dbf.find_data_single_lookup("biological_plate_data", temp_plate, "plate_name")[0]
                temp_run = row_data[2]
            if not temp_assay:
                row_data = dbf.find_data_single_lookup("assay_runs", temp_assay, "run_name")[0]
                temp_assay = row_data[2]
            window["-BIO_INFO_ASSAY_DROPDOWN-"].update(value=temp_assay)
            bio_info_window_update(dbf, window, values)
            window["-BIO_INFO_RUN_DROPDOWN-"].update(value=temp_run)
            bio_info_plate_list_update(dbf, window, values, None)
            window["-BIO_INFO_PLATES_DROPDOWN-"].update(value=temp_plate)
            well_dict_bio_info = bio_info_plate_update(dbf, config, window, values, event, well_dict_bio_info,
                                                       plate=temp_plate)
        elif temp_run:
            if not temp_assay:
                row_data = dbf.find_data_single_lookup("assay_runs", temp_assay, "run_name")[0]
                temp_assay = row_data[2]
            window["-BIO_INFO_ASSAY_DROPDOWN-"].update(value=temp_assay)
            bio_info_window_update(dbf, window, values, assay=temp_assay)
            window["-BIO_INFO_RUN_DROPDOWN-"].update(value=temp_run)
            bio_info_plate_list_update(dbf, window, values, ["All", temp_run])

        else:
            window["-BIO_INFO_ASSAY_DROPDOWN-"].update(value=temp_assay)
            bio_info_window_update(dbf, window, values, assay=temp_assay)
            # window["-BIO_INFO_PLATES_DROPDOWN-"].update(value="All")

    window["-TAB_GROUP_TWO-"].Widget.select(1)
    return well_dict_bio_info


def compound_table_double_click(dbf, config, window, values, event):

    if event == "-BIO_EXP_COMPOUND_TABLE-+-double click-":
        try:
            table_row = values["-BIO_EXP_COMPOUND_TABLE-"][0]
        except IndexError:
            temp_compound_id = None
        else:
            temp_table_data = window["-BIO_EXP_COMPOUND_TABLE-"].get()
            temp_compound_id = temp_table_data[table_row][0]
    elif event == "-PLATE_TABLE_TABLE-+-double click-":
        try:
            table_row = values["-PLATE_TABLE_TABLE-"][0]
        except IndexError:
            temp_compound_id = None
        else:
            temp_table_data = window["-PLATE_TABLE_TABLE-"].get()
            temp_compound_id = temp_table_data[table_row][2]
    elif event == "-SUB_SEARCH_TABLE-+-double click-":
        try:
            table_row = values["-SUB_SEARCH_TABLE-"][0]
        except IndexError:
            temp_compound_id = None
        else:
            temp_table_data = window["-SUB_SEARCH_TABLE-"].get()
            temp_compound_id = temp_table_data[table_row][0]
    elif event == "-MAIN_COMPOUND_TABLE-+-double click-":
        try:
            window.Element("-MAIN_COMPOUND_TABLE-").SelectedRows[0]
        except IndexError:
            temp_compound_id = None
        else:
            temp_compound_id = window.Element("-MAIN_COMPOUND_TABLE-").SelectedRows[0]
    else:
        temp_compound_id = None

    if temp_compound_id:
        window["-TAB_GROUP_TWO-"].Widget.select(0)
        window["-COMPOUND_INFO_ID-"].update(value=temp_compound_id)
        print(temp_compound_id)
        update_overview_compound(dbf, config, window, temp_compound_id)
    else:
        pass


def _update_bio_exp_plate_table(config, window, values, approval_check=False):
    table_name = "biological_plate_data"
    bio_exp_selected_runs = []
    temp_table_data = window["-BIO_EXP_ASSAY_RUN_TABLE-"].get()
    for temp_values in values["-BIO_EXP_ASSAY_RUN_TABLE-"]:
        bio_exp_selected_runs.append(temp_table_data[temp_values][0])
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

    return table_data


def _update_bio_exp_compound_table(config, window, values,
                                   threshold=False, compound_amount=False, approval_check=False, hit=False):

    temp_plate_table_data = window["-BIO_EXP_PLATE_TABLE-"].get()
    bio_exp_selected_plates = []
    for temp_values in values["-BIO_EXP_PLATE_TABLE-"]:
        bio_exp_selected_plates.append(temp_plate_table_data[temp_values][0])

    table_name = "biological_compound_data"
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
            print("No bio_exp_compound_data - for update bio exp table")
            print(bio_exp_compound_data)

        # Update compound table data and counter:
        window["-BIO_EXP_COMPOUND_TABLE-"].update(values=table_data)
        window["-BIO_EXP_COMPOUND_COUNTER-"].update(value=len(table_data))
    else:
        table_data = None

    # update the note field for the plate table:
    temp_table_data = window["-BIO_EXP_PLATE_TABLE-"].get()
    try:
        temp_table_data[values["-BIO_EXP_PLATE_TABLE-"][0]]
    except KeyError as e:
        print(e)
        print(values["-BIO_EXP_PLATE_TABLE-"])
        note = ""
    else:
        temp_plate = temp_table_data[values["-BIO_EXP_PLATE_TABLE-"][0]][0]
        dbf = DataBaseFunctions(config)
        temp_data = dbf.find_data_single_lookup("biological_plate_data", temp_plate, "plate_name")
        note = temp_data[0][9]

    window["-BIO_EXP_PLATE_NOTE-"].update(value=note)


def bio_exp_compound_list(config, window, event, values):

    # Sort the list based on the score
    temp_table_data = window["-BIO_EXP_COMPOUND_TABLE-"].get()
    temp_table_data.sort(key=lambda x: x[1])

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
    for rows, row_data in enumerate(temp_table_data):
        if threshold and row_data[1] > threshold:
            break

        bio_list.append(row_data[0])
        if compound_amount and rows + 1 >= compound_amount:
            break
    print("testing!!! ")

# ToDo Write information for the whole dose-response module!!!


if __name__ == "__main__":
    pass
