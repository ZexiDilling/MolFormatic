import configparser

import config_writer
from database_functions import _get_list_of_names_from_database_double_lookup
from database_handler import DataBaseFunctions
from start_up_values import database_guard

config = configparser.ConfigParser()
config.read("config.ini")

db_active = database_guard(config, config_writer)
dbf = DataBaseFunctions(config)

column_headline = "layout_name"
limiting_value = "single"
limiting_header = "style"
table = "plate_layout"

plate_list = _get_list_of_names_from_database_double_lookup(dbf, table, column_headline, limiting_value,
                                                            limiting_header)


