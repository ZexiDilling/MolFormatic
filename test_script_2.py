from gui_functions import get_plate_layout

if __name__ == "__main__":
    import configparser
    config = configparser.ConfigParser()
    config.read("config.ini")

    plate_list, archive_plates_dict = get_plate_layout(config)
    plate_to_layout = {"plate_1": "Daniells_Alpha_So"}
    barcode = "plate_1"
    print(archive_plates_dict[plate_to_layout[barcode]]["sample"])
