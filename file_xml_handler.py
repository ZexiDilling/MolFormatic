import xml.etree.ElementTree as ET
from datetime import date
from pathlib import Path


def _convert_dict(temp_dict):
    """
    Convert one dict to a new one

    :param temp_dict: A temporary dict
    :type temp_dict: dict
    :return: The same dict in a new formate that fits with the database
    :rtype: dict
    """
    data_dict = {}
    for transferee in temp_dict:
        counter = 0
        for transferee_counter in temp_dict[transferee]["transferees"]:
            temp_name = f"{transferee}_{counter}"
            data_dict[temp_name] = {}

            # for info in data_dict[transferee]["transferees"][transferee_counter]:
            data_dict[temp_name]["SourceBarcode"] = temp_dict[transferee]["source"]
            data_dict[temp_name]["SourceWell"] = temp_dict[transferee]["transferees"][transferee_counter]["source_well"]
            data_dict[temp_name]["DestinationBarcode"] = temp_dict[transferee]["destination"]
            data_dict[temp_name]["DestinationWell"] = temp_dict[transferee]["transferees"][transferee_counter]["destination_well"]
            data_dict[temp_name]["Volume"] = temp_dict[transferee]["transferees"][transferee_counter]["transferee_volume"]
            counter += 1
    return data_dict


def _get_transferee_dict(file_list):
    """
    Translate XML file_lise in to two dict

    :param file_list: A list of files
    :type file_list: list
    :return:
        - transferee: what have been transfereed
        - destination_plates: What plate the transferee have gone to
    :rtype:
        - dict
        - dict

    """

    transferee = {}
    destination_plates = {}
    # source_plate = {}

    for i in file_list:
        doc = ET.parse(i)
        root = doc.getroot()

        for dates in root.iter("transfer"):
            date_running = dates.get("date")
            date_str = f"plate_production_{date_running}"
            transferee[date_str] = {}

        # finds barcode for source and destination
        for plates in root.iter("plate"):
            source_destination = plates.get("type")
            barcode = plates.get("barcode")
            transferee[date_str][source_destination] = barcode


            # if plates.get("type") == "source":
            #     source_plate[barcode] = {}
            #     source_plate[barcode]["SourceBarcode"] = barcode
            #     source_plate[barcode]["date"] = date.today()

            if plates.get("type") == "destination":
                destination_plates[barcode] = {}
                destination_plates[barcode]["DestinationBarcode"] = barcode
                destination_plates[barcode]["date"] = date.today()

        # find source, destination and volume for each transferee
        for wells_t in root.iter("printmap"):
            wells_transferee = int(wells_t.get("total"))
            transferee[date_str]["transferees"] = {}
            for counter in range(wells_transferee):
                temp_str = f"Transferee_{counter + 1}"
                transferee[date_str]["transferees"][temp_str] = {}

                wells_source = wells_t[counter].get("n")
                wells_destination = wells_t[counter].get("dn")
                transferee_volume = float(wells_t[counter].get("vt")) * 10e-6

                transferee[date_str]["transferees"][temp_str]["source_well"] = wells_source
                transferee[date_str]["transferees"][temp_str]["destination_well"] = wells_destination
                transferee[date_str]["transferees"][temp_str]["transferee_volume"] = transferee_volume

        # find source, destination and reason for each skipped well
        for wells in root.iter("skippedwells"):
            wells_skipped = int(wells.get("total"))
            transferee[date_str]["Skipped"] = {}

            # finds destination and source wells data
            for z in range(wells_skipped):
                temp_str = f"Skipped_{z + 1}"
                transferee[date_str]["Skipped"][temp_str] = {}
                wells_destination = wells[z].get("dn")
                wells_source = wells[z].get("n")
                reason = wells[z].get("reason")

                transferee[date_str]["Skipped"][temp_str]["source_well"] = wells_source
                transferee[date_str]["Skipped"][temp_str]["destination_well"] = wells_destination
                transferee[date_str]["Skipped"][temp_str]["reason"] = reason

    return transferee, destination_plates


def xml_controller(file_list):
    """
    Controls the XML reader

    :param file_list: List of files with XML data
    :type file_list: list
    :return:
        - transferee: what have been transfereed
        - destination_plates: What plate the transferee have gone to
    :rtype:
        - dict
        - dict
    """

    transferee_dict, destination_plates = _get_transferee_dict(file_list)
    data_dict = _convert_dict(transferee_dict)

    return data_dict, destination_plates


def convert_echo_to_db(files):

    echo_to_db = {}
    transfer_counter = 0
    for file_index, files in enumerate(files):
        files = Path(files)
        if files.name.startswith("Transfer"):
            doc = ET.parse(files)
            root = doc.getroot()
            # for counting plates and transferees
            for plates in root.iter("plate"):
                barcode = plates.get("barcode")
                source_destination = plates.get("type")

                if source_destination == "destination":
                    temp_d_barcode = barcode
                if source_destination == "source":
                    temp_s_barcode = barcode

            try:
                echo_to_db[temp_d_barcode]
            except KeyError:
                echo_to_db[temp_d_barcode] = {"skipped_wells": {},
                                              "transferred_wells": {}}
            for wells in root.iter("printmap"):
                wells_transferred = wells.get("total")
                if int(wells_transferred) != 0:
                    for z in range(int(wells_transferred)):
                        destination_well = wells[z].get("dn")
                        source_well = wells[z].get("n")
                        vol = wells[z].get("vt")
                        echo_to_db[temp_d_barcode]["transferred_wells"][destination_well] = {
                            "mp_source_plate": temp_s_barcode,
                            "mp_source_well": source_well,
                            "vol": vol}

            for wells in root.iter("skippedwells"):
                wells_skipped = wells.get("total")
                if int(wells_skipped) != 0:
                    transfer_counter += 1

                    for z in range(int(wells_skipped)):
                        destination_well = wells[z].get("dn")
                        source_well = wells[z].get("n")
                        reason = wells[z].get("reason")

                        reason = reason.split(":")[0]
                        echo_to_db[temp_d_barcode]["skipped_wells"][destination_well] = {
                            "mp_source_plate": temp_s_barcode,
                            "mp_source_well": source_well,
                            "reason": reason}

    return echo_to_db


if __name__ == "__main__":
    path = "2022-03-03"
    from file_handler import get_file_list

    file_list = get_file_list(path)

    data, test = xml_controller(file_list)

    print(data)
