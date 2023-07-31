from gui_functions import grab_table_data
from database_controller import DataBaseFunctions


def update_bio_exp_assay_table(config, sg, window, event, values):

    table_name = "assay_runs"
    selected_assays = values["-BIO_EXP_TABLE_ASSAY_LIST_BOX-"]
    selected_headlines = ["run_name", "batch", "date", "note"]
    search_list_clm = "assay_name"
    bio_exp_assay_runs, _ = grab_table_data(config, table_name, selected_assays,
                                            specific_rows=selected_headlines, search_list_clm=search_list_clm)

    window["-BIO_EXP_ASSAY_RUN_TABLE-"].update(values=bio_exp_assay_runs)

    return bio_exp_assay_runs


def update_bio_exp_plate_table(config, sg, window, event, values, bio_exp_assay_runs, approval_check=False):
    table_name = "biological_plate_data"
    bio_exp_selected_runs = []
    for temp_values in values["-BIO_EXP_ASSAY_RUN_TABLE-"]:
        bio_exp_selected_runs.append(bio_exp_assay_runs[temp_values][0])
    search_list_clm = "assay_run"

    selected_headlines = ["plate_name", "z_prime", "approval", "note", "analysed_method", "assay_run"]

    bio_exp_plate_data, _ = grab_table_data(config, table_name, bio_exp_selected_runs,
                                            specific_rows=selected_headlines, search_list_clm=search_list_clm)

    table_data = []
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

    return table_data


def update_bio_exp_compound_table(config, sg, window, event, values, bio_exp_plate_data,
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

    # update the note field for the plate table:
    try:
        bio_exp_plate_data[values["-BIO_EXP_PLATE_TABLE-"][0]]
    except:
        print(values["-BIO_EXP_PLATE_TABLE-"])
        note = ""
    else:
        temp_plate = bio_exp_plate_data[values["-BIO_EXP_PLATE_TABLE-"][0]][0]
        dbf = DataBaseFunctions(config)
        temp_data = dbf.find_data_single_lookup("biological_plate_data", temp_plate, "plate_name")
        note = temp_data[0][9]
    window["-BIO_EXP_PLATE_NOTE-"].update(value=note)


    return table_data


