from PySimpleGUI import PopupError, Popup

from csv_handler import CSVWriter
from database_functions import get_number_of_rows, update_database, database_to_table
from excel_handler import plate_dilution_write_vol_well_amount, plate_dilution_excel
from plate_dilution import PlateDilution
from plate_formatting import plate_layout_re_formate
from upstarts_values import window_1_extra


def method_do_update(window, values):
    if values["-PD_METHOD_DD-"] == "Generate":
        window["-PD_SAVE_PLATES-"].update(disabled=False)
        window["-PD_WELL_AMOUNT-"].update(disabled=False)
        window["-PD_ADD_SOURCE_WELLS-"].update(value=False, disabled=True)
        window["-DP_PLATE_LAYOUT-"].update(disabled=False)
        window["-DP_SOURCE_FILE_GENERATE-"].update(disabled=False)
        window["-PD_WELL_LAYOUT-"].update(disabled=False)
    elif values["-PD_METHOD_DD-"] == "Calculate":
        window["-PD_SAVE_PLATES-"].update(disabled=True)
        window["-PD_WELL_AMOUNT-"].update(disabled=True)
        window["-PD_ADD_SOURCE_WELLS-"].update(disabled=False)
        window["-PD_SOURCE_WELL_AMOUNT-"].update(disabled=True)
        window["-DP_PLATE_LAYOUT-"].update(disabled=True)
        window["-DP_SOURCE_FILE_GENERATE-"].update(disabled=True)
        window["-PD_WELL_LAYOUT-"].update(disabled=True)


def add_source_wells_update(window, values):
    window["-PD_SOURCE_WELL_AMOUNT-"].update(disabled=not values["-PD_ADD_SOURCE_WELLS-"])
    window["-PD_SOURCE_WELL_AMOUNT_INPUT-"].update(disabled=not values["-PD_ADD_SOURCE_WELLS-"])
    window["-PD_WELL_LAYOUT-"].update(disabled=not values["-PD_ADD_SOURCE_WELLS-"])


def execute_button_pressed(config, window, values, archive_plates_dict):
    if not values["-PD_FILE-"]:
        PopupError("Please select an Excel file")
    elif values["-PD_METHOD_DD-"] == "Generate" and not values["-DP_PLATE_LAYOUT-"]:
        PopupError("Please select a Plate-Layout")
    else:
        function = values["-PD_METHOD_DD-"]
        file = values["-PD_FILE-"]
        if values["-PD_WELL_AMOUNT-"] == "Input":
            dw_amount = int(values["-PD_WELL_AMOUNT_INPUT-"])
        else:
            dw_amount = values["-PD_WELL_AMOUNT-"]
        save_plates = values["-PD_SAVE_PLATES-"]
        well_layout = values["-PD_WELL_LAYOUT-"]
        try:
            plate_layout = archive_plates_dict[values["-DP_PLATE_LAYOUT-"]]
        except KeyError:
            plate_layout = None
        add_source_wells = values["-PD_ADD_SOURCE_WELLS-"]
        if values["-PD_SOURCE_WELL_AMOUNT-"] == "Minimum":
            source_well_amount = values["-PD_SOURCE_WELL_AMOUNT-"]
        else:
            source_well_amount = int(values["-PD_SOURCE_WELL_AMOUNT_INPUT-"])
        pb_source_file = values["-DP_SOURCE_FILE_GENERATE-"]

        state = _plate_dilution(config, function, file, dw_amount, add_source_wells, source_well_amount,
                                save_plates, well_layout, plate_layout, pb_source_file)

        Popup(state)


def database_tab_pressed(config, window, values, db_active):
    if values["-EXTRA_SUB_DATABASE_TABS-"] == "Responsible" and not window_1_extra["responsible_data"]:
        table = "responsible"
        headings = config["Extra_tab_database_headings"]["responsible"].split(",")
        if db_active:
            responsible_data = database_to_table(config, table, headings)
        else:
            responsible_data = []
        window["-EXTRA_DATABASE_RESPONSIBLE_TABLE-"].update(values=[responsible_data])
        print(f"responsible_data: {responsible_data}")
    elif values["-EXTRA_SUB_DATABASE_TABS-"] == "Customers" and not window_1_extra["customers_data"]:
        table = "customers"
        headings = config["Extra_tab_database_headings"]["customers"].split(",")
        customers_data = database_to_table(config, table, headings)
        window["-EXTRA_DATABASE_CUSTOMERS_TABLE-"].update(values=[customers_data])
        print(f"customers_data: {customers_data}")
    elif values["-EXTRA_SUB_DATABASE_TABS-"] == "Vendors" and not window_1_extra["vendors_data"]:
        table = "vendors"
        headings = config["Extra_tab_database_headings"]["vendors"].split(",")
        vendors_data = database_to_table(config, table, headings)
        window["-EXTRA_DATABASE_VENDORS_TABLE-"].update(values=[vendors_data])
        print(f"vendors_data: {vendors_data}")
    elif values["-EXTRA_SUB_DATABASE_TABS-"] == "AC" and not window_1_extra["ac_data"]:
        table = "origin"
        headings = config["Extra_tab_database_headings"]["origin"].split(",")
        ac_data = database_to_table(config, table, headings)
        window["-EXTRA_DATABASE_AC_TABLE-"].update(values=[ac_data])
        print(f"ac_data/origine: {ac_data}")
    elif values["-EXTRA_SUB_DATABASE_TABS-"] == "Plate Types" and not window_1_extra["plate_types_data"]:
        table = "plate_type"
        headings = config["Extra_tab_database_headings"]["plate_type"].split(",")
        plate_types_data = database_to_table(config, table, headings)
        window["-EXTRA_PLATE_TYPE_LISTBOX-"].update(values=[plate_types_data])

        vendors = database_to_table(config, "vendors", "name")
        window["-EXTRA_PLATE_TYPE_VENDOR-"].update(values=vendors)
        print(f"plate_types_data: {plate_types_data}")

    elif values["-EXTRA_SUB_DATABASE_TABS-"] == "Location" and not window_1_extra["location_data"]:
        table = "locations"
        headings = config["Extra_tab_database_headings"]["locations"].split(",")
        location_data = database_to_table(config, table, headings)
        window["-EXTRA_DATABASE_LOCATIONS_TABLE-"].update(values=[location_data])
        print(f"location_data: {location_data}")


def database_responsible_import(config, window, values):
    if not values["-EXTRA_DATABASE_RESPONSIBLE_NAME-"]:
        PopupError("Please fill in name")
    elif not values["-EXTRA_DATABASE_RESPONSIBLE_E_MAIL-"]:
        PopupError("Please fill in E-mail")
    else:
        table = "responsible"
        row_id = get_number_of_rows(config, table) + 1
        name = values["-EXTRA_DATABASE_RESPONSIBLE_NAME-"]
        e_mail = values["-EXTRA_DATABASE_RESPONSIBLE_E_MAIL-"]
        info = values["-EXTRA_DATABASE_RESPONSIBLE_INFO-"]
        data = {"row_id": row_id,
                "name": name,
                "e_mail": e_mail,
                "info": info}
        update_database(config, table, data)


def database_customers_import(config, window, values):
    if not values["-EXTRA_DATABASE_CUSTOMERS_NAME-"]:
        PopupError("Please fill in name")
    elif not values["-EXTRA_DATABASE_CUSTOMERS_E_MAIL-"]:
        PopupError("Please fill in E-mail")
    else:
        table = "customers"
        row_id = get_number_of_rows(config, table) + 1
        name = values["-EXTRA_DATABASE_CUSTOMERS_NAME-"]
        e_mail = values["-EXTRA_DATABASE_CUSTOMERS_E_MAIL-"]
        info = values["-EXTRA_DATABASE_CUSTOMERS_INFO-"]
        data = {"row_id": row_id,
                "name": name,
                "e_mail": e_mail,
                "info": info}
        update_database(config, table, data)


def database_vendors_import(config, window, values):
    if not values["-EXTRA_DATABASE_VENDORS_NAME-"]:
        PopupError("Please fill in name")
    elif not values["-EXTRA_DATABASE_VENDORS_E_MAIL-"]:
        PopupError("Please fill in E-mail")
    else:
        table = "vendors"
        row_id = get_number_of_rows(config, table) + 1
        name = values["-EXTRA_DATABASE_VENDORS_NAME-"]
        e_mail = values["-EXTRA_DATABASE_VENDORS_E_MAIL-"]
        info = values["-EXTRA_DATABASE_VENDORS_INFO-"]
        data = {"row_id": row_id,
                "name": name,
                "e_mail": e_mail,
                "info": info}
        update_database(config, table, data)


def database_academia_company_import(config, window, values):
    if not values["-EXTRA_DATABASE_AC_NAME-"]:
        PopupError("Please fill in name")
    elif not values["-EXTRA_DATABASE_AC_AC-"]:
        PopupError("Please Choose A/C")
    else:
        table = "origin"
        row_id = get_number_of_rows(config, table) + 1
        origin = values["-EXTRA_DATABASE_AC_NAME-"]
        ac = values["-EXTRA_DATABASE_AC_AC-"]
        contact_person = values["-EXTRA_DATABASE_AC_CONTACT_NAME-"]
        e_mail = values["-EXTRA_DATABASE_AC_E_MAIL-"]
        info = values["-EXTRA_DATABASE_AC_INFO-"]
        data = {"ac_id": row_id,
                "origin": origin,
                "ac": ac,
                "contact_person": contact_person,
                "e_mail": e_mail,
                "info": info}
        update_database(config, table, data)


def database_place_type_import(config, window, values):
    if not values["-EXTRA_PLATE_TYPE_NAME-"]:
        PopupError("Please fill in Plate Type")
    elif not values["-EXTRA_PLATE_TYPE_VENDOR-"]:
        PopupError("Please Choose a Vendor")
    elif not values["-EXTRA_PLATE_TYPE_PRODUCT_NUMBER-"]:
        PopupError("Please fill in product number")
    else:
        table = "plate_type"
        try:
            data = {"row_id": "",
                    "plate_type": values["-EXTRA_PLATE_TYPE_NAME-"],
                    "vendor": values["-EXTRA_PLATE_TYPE_VENDOR-"][0],
                    "product_number": values["-EXTRA_PLATE_TYPE_PRODUCT_NUMBER-"],
                    "sterile": values["-EXTRA_PLATE_TYPE_STERILE-"],
                    "info": values["-EXTRA_PLATE_TYPE_INFO-"],
                    "size": int(values["-EXTRA_PLATE_TYPE_SIZE-"]),
                    "well_offset_x": float(values["-EXTRA_PLATE_TYPE_WELL_OFFSET_X-"]),
                    "well_offset_y": float(values["-EXTRA_PLATE_TYPE_WELL_OFFSET_Y-"]),
                    "well_spacing_x": float(values["-EXTRA_PLATE_TYPE_WELL_SPACING_X-"]),
                    "well_spacing_y": float(values["-EXTRA_PLATE_TYPE_WELL_SPACING_Y-"]),
                    "plate_height": float(values["-EXTRA_PLATE_TYPE_PLATE_HEIGHT-"]),
                    "plate_height_lid": float(values["EXTRA_PLATE_TYPE_PLATE_HEIGHT_LID-"]),
                    "flang_height": float(values["-EXTRA_PLATE_TYPE_FLANGE_HEIGHT-"]),
                    "well_depth": (float(values["-EXTRA_PLATE_TYPE_PLATE_HEIGHT-"]) -
                                   float(values["-EXTRA_PLATE_TYPE_FLANGE_HEIGHT-"])),
                    "well_width": float(values["-EXTRA_PLATE_TYPE_WELL_WIDTH-"]),
                    "max_volume": float(values["-EXTRA_PLATE_TYPE_MAX_VOL-"]),
                    "working_volume": float(values["-EXTRA_PLATE_TYPE_WORKING_VOL-"]),
                    "dead_volume": float(values["-EXTRA_PLATE_TYPE_WELL_DEAD_VOL-"])}
        except ValueError:
            data = {"row_id": 0,
                    "plate_type": values["-EXTRA_PLATE_TYPE_NAME-"],
                    "vendor": values["-EXTRA_PLATE_TYPE_VENDOR-"][0],
                    "product_number": values["-EXTRA_PLATE_TYPE_PRODUCT_NUMBER-"],
                    "sterile": values["-EXTRA_PLATE_TYPE_STERILE-"],
                    "info": values["-EXTRA_PLATE_TYPE_INFO-"],
                    "well_offset_x": 0,
                    "well_offset_y": 0,
                    "well_spacing_x": 0,
                    "well_spacing_y": 0,
                    "plate_height": 0,
                    "plate_height_lid": 0,
                    "flang_height": 0,
                    "well_depth": 0,
                    "well_width": 0,
                    "max_volume": 0,
                    "working_volume": 0,
                    "dead_volume": 0}
        update_database(config, table, data)


def database_location_import(config, window, values):
    if not values["-EXTRA_DATABASE_LOCATIONS_ROOM-"]:
        PopupError("Please fill in Room")
    elif not values["-EXTRA_DATABASE_LOCATIONS_BUILDING-"]:
        PopupError("Please fill in Building")
    else:
        table = "locations"
        row_id = get_number_of_rows(config, table) + 1
        room = values["-EXTRA_DATABASE_LOCATIONS_ROOM-"]
        building = values["-EXTRA_DATABASE_LOCATIONS_BUILDING-"]
        spot = values["-EXTRA_DATABASE_LOCATIONS_SPOT-"]
        data = {"loc_id": row_id,
                "room": room,
                "building": building,
                "spot": spot}
        update_database(config, table, data)


def _plate_dilution(config, function, file, dw_amount, add_source_wells, source_well_amount, save_plates, well_layout,
                    plate_layout, pb_source_file):
    output_folder = config["folders"]["main_output_folder"]
    destination_plate_naming_scheme = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]

    if function == "Calculate":
        state = plate_dilution_write_vol_well_amount(config, file, dw_amount, add_source_wells, source_well_amount,
                                                     well_layout)
        return state

    elif function == "Generate":
        plate_layout = plate_layout_re_formate(config, plate_layout["well_layout"])
        sample_info_dict, replicate_samples_max, replicate_plate_sets, dilution_factor, concentration_counter, \
            control_vol, control_conc = plate_dilution_excel(file, save_plates, dw_amount, well_layout)

        pd = PlateDilution(config, save_plates, pb_source_file, control_vol, control_conc)

        sample_dict, dilution_dict = pd.pd_controller(sample_info_dict, replicate_samples_max, replicate_plate_sets,
                                                      dilution_factor, concentration_counter, plate_layout,
                                                      destination_plate_naming_scheme, dw_amount)


        csv_w = CSVWriter()
        if pb_source_file:
            csv_w.source_plate_dilution(dilution_dict, output_folder)

        csv_w.plate_dilution(sample_dict, output_folder)

        return "CSV file created"