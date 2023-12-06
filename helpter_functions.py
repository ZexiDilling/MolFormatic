from operator import itemgetter
import time
import PySimpleGUI as sg
from natsort import natsorted

from config_writer import ConfigWriter
from database_controller import FetchData


def int_guard(window, value, fall_back_value):
    if value:
        try:
            int(value)
        except ValueError:
            new_value = fall_back_value
            window[value].update(value=fall_back_value)
        else:
            new_value = int(value)

        return new_value
    else:
        return None


def config_update(config):
    fd = FetchData(config)
    cw = ConfigWriter(config)
    # database_specific_commercial
    search_limiter = {
        "ac": {"value": "Commercial",
               "operator": "=",
               "target_column": "ac",
               "use": True}}
    rows = fd.data_search("origin", search_limiter)
    simple_settings = {"database_specific_commercial": {},
                       "database_specific_academic": {}}
    for row in rows:
        simple_settings["database_specific_commercial"][f"vendor_{rows[row]['ac_id']}"] = rows[row]["origin"]

    # database_specific_academia
    search_limiter = {
        "ac": {"value": "Academic",
               "operator": "=",
               "target_column": "ac",
               "use": True}}
    rows = fd.data_search("origin", search_limiter)
    for row in rows:
        simple_settings["database_specific_commercial"][f"academia_{rows[row]['ac_id']}"] = rows[row]["origin"]

    cw.run(simple_settings, "simple_settings", True)


def sort_table(table, cols, reverse):
    """ sort a table by multiple columns
        table: a list of lists (or tuple of tuples) where each inner list
               represents a row
        cols:  a list (or tuple) specifying the column numbers to sort by
               e.g. (1,0) would sort by column 1, then by column 0
    """
    for col in reversed(cols):
        try:
            table = natsorted(table, key=itemgetter(col), reverse=reverse)
        except Exception as e:
            sg.popup_error('Error in sort_table', 'Exception in sort_table', e)
    reverse = not reverse
    return table, reverse


def time_translater(time_testing):
    time_return_value = time.strftime("%Hh%Mm%Ss", time.gmtime(time_testing))
    # time_return_value = time_testing
    return time_return_value


def eval_guard_dict(test_dict):
    """

    :param test_dict:
    :type test_dict: str
    :return:
    :rtype: dict or None
    """
    if test_dict.startswith("{") and test_dict.endswith("}"):
        return eval(test_dict)
    else:
        return None


def eval_guard_list(test_list):
    """

    :param test_dict:
    :type test_dict: str
    :return:
    :rtype: dict or None
    """
    if test_list.startswith("[") and test_list.endswith("]"):
        return eval(test_list)
    else:
        return None


def plate_layout_to_state_dict(plate_layout):
    state_dict = {"states": []}

    for counter in plate_layout:
        current_state = plate_layout[counter]["state"]
        well_id = plate_layout[counter]["well_id"]

        if current_state == "sample":
            use_it = True
        else:
            use_it = False

        if current_state not in state_dict["states"]:
            state_dict["states"].append(current_state)

        try:
            state_dict[current_state]
        except KeyError:
            state_dict[current_state] = {"use": use_it,
                                         "wells": [well_id]}
        else:
            state_dict[current_state]["wells"].append(well_id)

        state_dict[well_id] = current_state

    return state_dict
