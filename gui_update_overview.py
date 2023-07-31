from gui_functions import grab_table_data


def update_overview_compound(config, sg, window, event, values, compound_id):

    sample_row = grab_table_data(config, "compound_main", single_row=True, data_value=compound_id,
                                 headline="compound_id")
    if sample_row:
        # Get Academic/commercial information:
        ac_id = sample_row[0][6]
        all_data_ac = grab_table_data(config, "origin", single_row=True, data_value=ac_id, headline="ac_id")

        # Update info frame:
        window["-COMPOUND_INFO_AC-"].update(value=all_data_ac[0][3])
        window["-COMPOUND_INFO_ORIGIN-"].update(value=all_data_ac[0][2])
        window["-COMPOUND_INFO_ORIGIN_ID-"].update(value=sample_row[0][7])
        window["-COMPOUND_INFO_CONCENTRATION-"].update(value=sample_row[0][5])
        window["-COMPOUND_INFO_TUBE_VOLUME-"].update(value=sample_row[0][4])
        # ToDo find info from DB if that is there, else:
        window["-COMPOUND_INFO_TUBE_VOLUME-"].update(value="Missing info")

        # Update Picture frame:
        window["-COMPOUND_INFO_SMILES-"].update(value=sample_row[0][2])
        window["-COMPOUND_INFO_PIC-"].update(data=sample_row[0][3])

        mp_table_data = grab_table_data(config, "compound_mp", single_row=True, data_value=compound_id,
                                        headline="compound_id")

        dp_table_data = grab_table_data(config, "compound_dp", single_row=True, data_value=compound_id,
                                        headline="compound_id")
        assay_compound_table_data = grab_table_data(config, "biological_compound_data", single_row=True,
                                                    data_value=compound_id, headline="compound_id")

        assay_plate = []
        for assays in assay_compound_table_data:
            assay_plate.append(assays[3])
        assay_plate_table_data, _ = grab_table_data(config, "biological_plate_data", assay_plate,
                                                    specific_rows=None,
                                                    search_list_clm="plate_name")

        assay_runs = []
        for plates in assay_plate_table_data:
            assay_runs.append(plates[1])

        assay_run_table_data, _ = grab_table_data(config, "assay_runs", assay_runs,
                                                  specific_rows=None, search_list_clm="run_name")

        assays = []
        for runs in assay_run_table_data:
            assays.append(runs[1])
        assay_table_data, _ = grab_table_data(config, "assay", assays, specific_rows=None,
                                              search_list_clm="assay_name")

        # Get data:
        amount_mp = len(mp_table_data)
        amount_dp = len(dp_table_data)
        amount_assays = len(assay_compound_table_data)
        amount_assay_hits = 0
        amount_transfers = amount_dp + amount_assays
        amount_purity = "Nah"

        updated_mp_table_data = []
        updated_dp_table_data = []
        updated_assay_table_data = []
        updated_hit_table_data = []
        # updated_transfer_table_data = []
        updated_purify_table_data = []

        if len(mp_table_data) != 0:
            for rows, row_data in enumerate(mp_table_data):
                updated_mp_table_data.append([
                    row_data[2],
                    row_data[4],
                    row_data[5]
                ])
        else:
            updated_mp_table_data = [["No", "Data", "Found"]]

        if len(dp_table_data) != 0:
            for rows, row_data in enumerate(dp_table_data):
                updated_dp_table_data.append([
                    row_data[2],
                    row_data[4],
                    row_data[5]
                ])
        else:
            updated_dp_table_data = [["No", "Data", "Found"]]

        if len(assay_plate_table_data) != 0:
            for rows, row_data in enumerate(assay_plate_table_data):
                updated_assay_table_data.append([
                    row_data[9],
                    row_data[2],
                    row_data[10]
                ])
        else:
            updated_assay_table_data = [["No", "Data", "Found"]]

        if len(assay_compound_table_data) != 0:
            for rows, row_data in enumerate(assay_compound_table_data):
                if row_data[8] == 1:
                    amount_assay_hits += 1
                    temp_number = row_data[3].split("_")[-1]

                    assay = row_data[3].removesuffix(f"_{temp_number}")
                    updated_hit_table_data.append([
                        assay,
                        row_data[5],
                        row_data[7]
                    ])
        else:
            updated_hit_table_data = [["No", "Hits", "Found"]]

        updated_purify_table_data = [["Missing", "Coding", "TODO!!!"]]

        # Update overview:
        window["-COMPOUND_INFO_INFO_MP-"].update(value=amount_mp)
        window["-COMPOUND_INFO_INFO_DP-"].update(value=amount_dp)
        window["-COMPOUND_INFO_INFO_ASSAY-"].update(value=amount_assays)
        window["-COMPOUND_INFO_INFO_HITS-"].update(value=amount_assay_hits)
        window["-COMPOUND_INFO_INFO_TRANSFERS-"].update(value=amount_transfers)
        window["-COMPOUND_INFO_INFO_PURITY-"].update(value=amount_purity)

        # Update tables

        window["-COMPOUND_INFO_INFO_MP_TABLE-"].update(values=updated_mp_table_data)
        window["-COMPOUND_INFO_INFO_DP_TABLE-"].update(values=updated_dp_table_data)
        window["-COMPOUND_INFO_INFO_ASSAY_TABLE-"].update(values=updated_assay_table_data)
        window["-COMPOUND_INFO_INFO_HITS_TABLE-"].update(values=updated_hit_table_data)
        # window["-COMPOUND_INFO_INFO_TRANSFERS_TABLE-"].update(values=updated_transfer_table_data)
        window["-COMPOUND_INFO_INFO_PURITY_USED_TABLE-"].update(values=updated_purify_table_data)

    else:
        sg.PopupError(f'No compound data for:\n "{compound_id}"')