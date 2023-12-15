from database_functions import grab_table_data


def table_group_lcms(config, window, values):
    if values["-TABLE_TAB_GRP-"] == "LC Experiment table":
        # print("update listbox with data, if list box is empty")
        lc_exp_data, headlines = grab_table_data(config, "lc_experiment")
        window["-LC_MS_TABLE_BATCH_LIST_BOX-"].update(values=lc_exp_data)


def date_set_update(config, window, values):
    start_date = values["-LC_MS_TABLE_DATE_START_TARGET-"]
    end_date = values["-LC_MS_TABLE_DATE_END_TARGET-"]

    if start_date:
        use_start_date = True
    else:
        use_start_date = False
    if end_date:
        use_end_date = True
    else:
        use_end_date = False

    search_limiter = {
        "start_date": {"value": start_date, "operator": "<", "target_column": "date", "use": use_start_date},
        "end_date": {"value": end_date, "operator": ">", "target_column": "date", "use": use_end_date},
    }

    table_name = "lc_experiment"

    table_data, _ = grab_table_data(config, table_name, search_limiter)
    table_data_2, _ = grab_table_data(config, table_name)

    window["-LC_MS_TABLE_BATCH_LIST_BOX-"].update(values=table_data)


def batch_list_box_update(config, window, values):
    batch_date = values["-LC_MS_TABLE_BATCH_LIST_BOX-"]
    batch = []
    for data in batch_date:
        batch.append(data[0])
    if batch:
        tamp_table_data, _ = grab_table_data(config, "lc_raw")
        window["-LC_MS_SAMPLE_TABLE-"].update(values=tamp_table_data)

