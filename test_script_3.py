import configparser

from database_functions import DataBaseFunctions


if __name__ == "__main__":
    config = configparser.ConfigParser()
    config.read("config.ini")

    dbf = DataBaseFunctions(config)
    plate_data = dbf.find_data_single_lookup("plate_layout", "testing", "layout_name")[0][5]
    plate_layout = eval(plate_data[5])
    plate_size = int(plate_data[2].split("_")[-1])

    for rows in plate_layout:
        if plate_layout[rows]["state"] == "sample":
            print(plate_layout[rows])

