
from PySimpleGUI import PopupGetFile, popup_error, Popup, WIN_CLOSED
from helpter_functions import config_update, sort_table
from gui_help_info_layout import info_help_window
from gui_help_info_text import TextForInfo
from start_up_values import all_table_data


def menu_open(config, config_writer, window, gui_layout):
    db = PopupGetFile("Choose database", file_types=(("Database", ".db"),))
    setting_dict = {"Database": {"database": db}}
    config_writer.run(setting_dict, "simple_settings", True)
    config_update(config)
    window.close()
    return gui_layout.full_layout()


def menu_save():
    popup_error("Not working - Will never work... as all data is written to the DB... "
                   "unless I want to add a temp DB...")


def menu_about():
    with open("README.txt") as file:
        info = file.read()

    Popup(info)


def help_info_controller(config):

    window = info_help_window()
    tfi = TextForInfo(window, config)
    tfi.window_1_search_text()

    while True:
        event, values = window.read()
        if event == WIN_CLOSED:
            break

        if event == "-HELP_INFO_TAB_GROUPS-" and values["-HELP_INFO_TAB_GROUPS-"] == "          Search          " \
                and not values["-INFO_HELP_SEARCH-"]:
            tfi.window_1_search_text()
        if event == "-HELP_INFO_TAB_GROUPS-" and values["-HELP_INFO_TAB_GROUPS-"] == "        Bio Data          "\
                and not values["-INFO_HELP_BIO_DATA-"]:
            tfi.window_1_bio_data_text()
        if event == "-HELP_INFO_TAB_GROUPS-" and values["-HELP_INFO_TAB_GROUPS-"] == "       Purity Data       "\
                and not values["-INFO_HELP_PURITY_DATA-"]:
            tfi.window_1_purity_data_text()
        if event == "-HELP_INFO_TAB_GROUPS-" and values["-HELP_INFO_TAB_GROUPS-"] == "      Plate Layout      "\
                and not values["-INFO_HELP_PLATE_LAYOUT-"]:
            tfi.window_1_plate_layout_text()
        if event == "-HELP_INFO_TAB_GROUPS-" and values["-HELP_INFO_TAB_GROUPS-"] == "          Update          "\
                and not values["-INFO_HELP_UPDATE-"]:
            tfi.window_1_update_text()
        if event == "-HELP_INFO_TAB_GROUPS-" and values["-HELP_INFO_TAB_GROUPS-"] == "           Sim             "\
                and not values["-INFO_HELP_SIM-"]:
            tfi.window_1_sim_text()
        if event == "-HELP_INFO_TAB_GROUPS-" and values["-HELP_INFO_TAB_GROUPS-"] == "    Compound Info    "\
                and not values["-INFO_HELP_COMPOUND_INFO-"]:
            tfi.window_2_info_text()
        if event == "-HELP_INFO_TAB_GROUPS-" and values["-HELP_INFO_TAB_GROUPS-"] == "         Bio Info          "\
                and not values["-INFO_HELP_BIO_INFO-"]:
            tfi.window_2_bio_info_text()
        if event == "-HELP_INFO_TAB_GROUPS-" and values["-HELP_INFO_TAB_GROUPS-"] == "       Purity Info        "\
                and not values["-INFO_HELP_PURITY_INFO-"]:
            tfi.window_2_purity_info_text()
        if event == "-HELP_INFO_TAB_GROUPS-" and values["-HELP_INFO_TAB_GROUPS-"] == "   Compound Table   "\
                and not values["-INFO_HELP_COMPOUND_TABLE-"]:
            tfi.window_3_compound_table_text()
        if event == "-HELP_INFO_TAB_GROUPS-" and values["-HELP_INFO_TAB_GROUPS-"] == "         Bio Table        "\
                and not values["-INFO_HELP_EXP_TABLE-"]:
            tfi.window_3_bio_exp_table_text()
        if event == "-HELP_INFO_TAB_GROUPS-" and values["-HELP_INFO_TAB_GROUPS-"] == "       Plate Table       "\
                and not values["-INFO_HELP_PLATE_TABLE-"]:
            tfi.window_3_plate_table_text()
        if event == "-HELP_INFO_TAB_GROUPS-" and values["-HELP_INFO_TAB_GROUPS-"] == "         Glossary         "\
                and not values["-INFO_HELP_GLOSSARY-"]:
            tfi.glossary()


def sorting_the_tables(window, event, search_reverse):
    # TABLE CLICKED Event has value in format ('-TABLE=', '+CLICKED+', (row,col))
    # You can also call Table.get_last_clicked_position to get the cell clicked
    if event[0] in all_table_data and all_table_data[event[0]]:
        clicked_table = event[0]
        try:
            search_reverse[clicked_table]
        except KeyError:
            search_reverse[clicked_table] = {}
        if event[2][0] == -1 and event[2][1] != -1:  # Header was clicked and wasn't the "row" column
            col_num_clicked = event[2][1]
            try:
                search_reverse[clicked_table][col_num_clicked]
            except KeyError:
                search_reverse[clicked_table][col_num_clicked] = False

            new_table, search_reverse[clicked_table][col_num_clicked] = \
                sort_table(all_table_data[clicked_table][0:][:], (col_num_clicked, 0),
                           search_reverse[clicked_table][col_num_clicked])
            window[clicked_table].update(new_table)
            # all_table_data[clicked_table] = [all_table_data[clicked_table][0]] + new_table
            all_table_data[clicked_table] = new_table
