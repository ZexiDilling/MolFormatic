from operator import itemgetter

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