from natsort import natsorted

from database_functions import grab_table_data
from start_up_values import plate_table_table_heading_mp, plate_table_table_heading_dp, all_table_data


def table_group_tables(config, window, values):
    if values["-TABLE_TAB_GRP-"] == "Plate tables":
        temp_mp_plates, _ = grab_table_data(config, "mp_plates")
        mp_plates_list = []
        for rows in temp_mp_plates:
            mp_plates_list.append(rows[0])

        # sortes the table
        mp_plates_list = natsorted(mp_plates_list)

        window["-PLATE_TABLE_BARCODE_LIST_BOX-"].update(values=mp_plates_list)
        return temp_mp_plates


def clear_plate_table_update(config, window, values):
    window["-PLATE_TABLE_TABLE-"].update(values=[])
    window["-PLATE_TABLE_START_DATE_TARGET-"].update(value="")
    window["-PLATE_TABLE_END_DATE_TARGET-"].update(value="")
    window["-PLATE_TABLE_CHOOSER-"].update(value="Mother Plates")

    temp_mp_plates, _ = grab_table_data(config, "mp_plates")
    mp_plates_list = []
    for rows in temp_mp_plates:
        mp_plates_list.append(rows[0])

    # sortes the table
    mp_plates_list = natsorted(mp_plates_list)

    window["-PLATE_TABLE_BARCODE_LIST_BOX-"].update(values=mp_plates_list)
    window.Element("-PLATE_TABLE_TABLE-").Widget.configure(displaycolumns=plate_table_table_heading_mp)

    return temp_mp_plates


def plate_chooser_update(config, window, values):
    if values["-PLATE_TABLE_CHOOSER-"]:

        window["-PLATE_TABLE_TABLE-"].update(values=[])

        table_dict = {"Mother Plates": "mp_plates", "Daughter Plates": "dp_plates"}

        if values["-PLATE_TABLE_START_DATE_TARGET-"]:
            use_start_date = True
        else:
            use_start_date = False
        if values["-PLATE_TABLE_END_DATE_TARGET-"]:
            use_end_date = True
        else:
            use_end_date = False

        search_limiter = {
            "start_date": {"value": values["-PLATE_TABLE_START_DATE_TARGET-"], "operator": "<",
                           "target_column": "date", "use": use_start_date},
            "end_date": {"value": values["-PLATE_TABLE_END_DATE_TARGET-"], "operator": ">",
                         "target_column": "date", "use": use_end_date},
        }
        plate_data, _ = grab_table_data(config, table_dict[values["-PLATE_TABLE_CHOOSER-"]], search_limiter)
        if plate_data:
            plates = []
            for plate in plate_data:
                plates.append(plate[0])

            window["-PLATE_TABLE_BARCODE_LIST_BOX-"].update(values=plates)
        else:
            window["-PLATE_TABLE_BARCODE_LIST_BOX-"].update(values=[[]])

        if values["-PLATE_TABLE_CHOOSER-"] == "Mother Plates":
            window.Element("-PLATE_TABLE_TABLE-").Widget.configure(displaycolumns=plate_table_table_heading_mp)
        else:
            window.Element("-PLATE_TABLE_TABLE-").Widget.configure(displaycolumns=plate_table_table_heading_dp)


def barcode_list_box_update(config, window, values):
    table_dict = {"Mother Plates": {"clm": "mp_barcode", "table": "compound_mp"},
                  "Daughter Plates": {"clm": "dp_barcode", "table": "compound_dp"}}
    search_limiter = {"academic_commercial": {"value": values["-PLATE_TABLE_BARCODE_LIST_BOX-"],
                                              "operator": "IN",
                                              "target_column": table_dict[values["-PLATE_TABLE_CHOOSER-"]][
                                                  "clm"],
                                              "use": True}}

    all_table_data["-PLATE_TABLE_TABLE-"], _ = grab_table_data(config,
                                                               table_dict[values["-PLATE_TABLE_CHOOSER-"]][
                                                                   "table"], search_limiter)

    if values["-PLATE_TABLE_BARCODE_LIST_BOX-"]:
        window["-PLATE_TABLE_TABLE-"].update(values=all_table_data["-PLATE_TABLE_TABLE-"])
    else:
        window["-PLATE_TABLE_TABLE-"].update(values=[])


def table_limiter_update(config, window, values, temp_mp_plates):
    mp_limiter = values["-PLATE_TABLE_TEXT_LIMITER-"]
    mp_plates_list = []
    for rows in temp_mp_plates:
        if mp_limiter in rows[0].casefold():
            mp_plates_list.append(rows[0])

    window["-PLATE_TABLE_BARCODE_LIST_BOX-"].update(values=mp_plates_list)


