import configparser

from database_functions import DataBaseFunctions
from draw_tool_handler import launch_draw_tool

if __name__ == "__main__":
    config = configparser.ConfigParser()
    config.read("config.ini")

    dbf = DataBaseFunctions(config)

    row_data = dbf.find_data_single_lookup("plate_layout", "layout_name", "testing")
    print(row_data)





