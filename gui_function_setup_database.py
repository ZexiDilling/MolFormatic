from PySimpleGUI import popup_error, popup

from helpter_functions import config_update
from database_functions import update_database


def update_compound(config, window, values, gui_layout):
    if not values["-UPDATE_FOLDER-"]:
        popup_error("Please select a folder containing compound data")
    else:
        update_database(config, "compound_main", values["-UPDATE_FOLDER-"])
        config_update(config)
        window.close()
        return gui_layout.full_layout()


def update_plates(config, window, values, event):
    if not values["-UPDATE_FOLDER-"]:
        popup_error("Please select a folder containing the data")
    else:
        if event == "-UPDATE_MP-":
            update_database(config, "compound_mp", values["-UPDATE_FOLDER-"], "pb_mp_output")
        else:
            update_database(config, "compound_dp", values["-UPDATE_FOLDER-"], "echo_dp_out")
        popup("Done")


