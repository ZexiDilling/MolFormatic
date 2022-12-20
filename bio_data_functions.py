import PySimpleGUI
from openpyxl.formatting.rule import ColorScaleRule
from openpyxl.styles import PatternFill, Font
from openpyxl import load_workbook, Workbook
import pandas as pd
import re
from matplotlib import colors
from openpyxl_fix_functions import ex_cell
import datetime
from extra_functions import rgb_to_hex, hex_to_rgb


def org(all_data, well):
    """
    Original reading data from the platereader, for bio data

    :param all_data: All the data for the reading, including the state of the well following the plate layout, and the
        results from different calculations is added as they get to it.
    :type all_data: dict
    :param well: Witch well is being calculated
    :type well: str
    :return: All the readings from the original readings from the platreader in the same formate as the future methodes
    :rtype: dict
    """
    return all_data["plates"]["original"]["wells"][well]


def norm(all_data, well):
    """
    Normalises the data based on the avg of the minimum for bio data

    :param all_data: All the data for the reading, including the state of the well following the plate layout, and the
        results from different calculations is added as they get to it.
    :type all_data: dict
    :param well: Witch well is being calculated
    :type well: str
    :return: Returns a dict over all the values after the have been normalized.
    :rtype: dict
    """
    return all_data["plates"]["original"]["wells"][well] - all_data["calculations"]["original"]["minimum"]["avg"]


def pora(all_data, well):
    """
    Percentage of remaining activity calculation based on the avg of the max of normalised data

    :param all_data: All the data for the reading, including the state of the well following the plate layout, and the
        results from different calculations is added as they get to it.
    :type all_data: dict
    :param well: Witch well is being calculated
    :type well: str
    :return: Returns the Percentage of remaining activity based on the normalized data.
    :rtype: dict
    """
    return ((100 * all_data["plates"]["normalised"]["wells"][well]) /
                all_data["calculations"]["normalised"]["max"]["avg"])


def pora_internal(all_data, well):
    """
    Percentage of remaining activity calculation based on the avg of the max of normalised data
    This should properly be deleted as it is not needed...

    :param all_data: All the data for the reading, including the state of the well following the plate layout, and the
        results from different calculations is added as they get to it.
    :type all_data: dict
    :param well: Witch well is being calculated
    :type well: str
    :return: Returns the Percentage of remaining activity based on the original data.
    :rtype: dict
    """
    return ((100 * all_data["plates"]["normalised"]["wells"][well]) /
            all_data["calculations"]["original"]["max"]["avg"])


def z_prime_calculator(all_data, method):
    """
    Calculate Z-prime

    :param all_data: All the data for the reading, including the state of the well following the plate layout, and the
        results from different calculations is added as they get to it.
    :type all_data: dict
    :param method: The heading for the calculations that is not avg and stdev
    :type method: str
    :return: Returns the Z-Prime value
    :rtype: int
    """
    return 1 - ((3 * (all_data["calculations"][method]["max"]["stdev"] +
                (all_data["calculations"][method]["minimum"]["stdev"]))) /
                abs(all_data["calculations"][method]["max"]["avg"] +
                (all_data["calculations"][method]["minimum"]["avg"])))


def state_mapping(config, ws, translate_wells_to_cells, plate, init_row, free_col, temp_dict, methode):
    """
    Colour in the state of the wells and write a guide in the site to translate
    Might need to re-writes this modul, to exclude temp_dict, as it should not be needed, as all the information should
    come from the plate layout

    :param config: The config file
    :type config: configparser.ConfigParser
    :param ws: Witch worksheet that mappings is conected too
    :type ws: openpyxl.worksheet.worksheet.Worksheet
    :param translate_wells_to_cells: A dict containing what well in the plate belongs to what cell in the excel file
    :type translate_wells_to_cells: dict
    :param plate: Plate layout for what well have what state
    :type plate: dict
    :param init_row: Row to start writing to
    :type init_row: int
    :param free_col: What column should be free / the first column after the last column used for the plate data
    :type free_col: int
    :param temp_dict: The data that is being analysed for the state mapping
    :type temp_dict: dict
    :param methode: The Method that is being used
    :type methode: str
    :return: The colouring of the cells, and a reading guide for the colours, in the excel ark
    """
    # colour the wells
    init_row_start = init_row
    for counter in plate["well_layout"]:
        state = plate["well_layout"][counter]["state"]
        cell_color = config["plate_colouring"][state]
        # cell_color = colors.cnames[colour]
        cell_color = cell_color.replace("#", "")
        temp_cell = translate_wells_to_cells[plate["well_layout"][counter]["well_id"]]
        ws[temp_cell].fill = PatternFill("solid", fgColor=cell_color)
    for state in temp_dict["plates"][methode]:
        # writes the colour guide
        if state != "wells":
            if init_row_start == init_row:
                ws[ex_cell(init_row + 1, free_col)] = "well state"
                ws[ex_cell(init_row + 1, free_col + 1)] = "colour coding"
                ws[ex_cell(init_row + 1, free_col)].font = Font(b=True)
                ws[ex_cell(init_row + 1, free_col + 1)].font = Font(b=True)
            ws[ex_cell(init_row + 2, free_col)] = state
            ws[ex_cell(init_row + 2, free_col + 1)] = config["plate_colouring"][state]
            init_row += 1


def heatmap(config, ws, pw_dict, translate_wells_to_cells, heatmap_colours):
    """
    Colour code based on values.

    :param ws: The worksheet where the heat map is placed
    :type ws: openpyxl.worksheet.worksheet.Worksheet
    :param pw_dict: Dict for each well and what state it is (sample, blank, empty....
    :type pw_dict: dict
    :param translate_wells_to_cells: Dict for each well cells value
    :type translate_wells_to_cells: dict
    :return: A heatmap coloured depending on the options, in the excel file.
    """
    temp_list = []
    for well in pw_dict:
        if pw_dict[well] == "sample":
            temp_list.append(well)

    # 3 colours heat mao
    ws.conditional_formatting.add(f"{translate_wells_to_cells[temp_list[0]]}:"
                                  f"{translate_wells_to_cells[temp_list[-1]]}",
                                  ColorScaleRule(start_type='percentile', start_value=10,
                                                 start_color=config["Settings_bio"]["plate_report_heatmap_colours_low"].removeprefix("#"),
                                                 mid_type='percentile', mid_value=50,
                                                 mid_color=config["Settings_bio"]["plate_report_heatmap_colours_mid"].removeprefix("#"),
                                                 end_type='percentile',  end_value=90,
                                                 end_color=config["Settings_bio"]["plate_report_heatmap_colours_high"].removeprefix("#"))
                                  )
    # 2 colours heat mao
    # ws.conditional_formatting.add(f"{translate_wells_to_cells[temp_list[0]]}:"
    #                               f"{translate_wells_to_cells[temp_list[-1]]}",
    #                               ColorScaleRule(start_type="min",
    #                                              start_color=config["colours to hex"][heatmap_colours["start"]],
    #                                              end_type="max",
    #                                              end_color=config["colours to hex"][heatmap_colours["end"]]))


def hit_mapping(ws, temp_dict, pora_threshold, methode, translate_wells_to_cells, free_col, init_row):
    """
    Colour coding a plate depending on hits (the values are selected before hand), and writes the bounderies for the
    hits.

    :param ws: The worksheet where the hit map is
    :type ws: openpyxl.worksheet.worksheet.Worksheet
    :param temp_dict: A dicts from all_data, that is only containing the values for the plate. These well values are
        the basis for the hits
    :type temp_dict: dict
    :param pora_threshold: The hit map threshold / bounderies
    :type pora_threshold: dict
    :param methode: The method that is being used
    :type methode: str
    :param translate_wells_to_cells: A dict containing what well in the plate belongs to what cell in the excel file
    :type translate_wells_to_cells: dict
    :param free_col: The first free column after the plate.
    :type free_col: int
    :param init_row: The row where the plate can state to be writen in the excel sheet.
    :type init_row: int
    :return: Colouring the cells in the excel file, depending on the boundaries of the hit.
    """
    # Colour the wells:
    for wells in translate_wells_to_cells:
        for split in pora_threshold:
            if split != "colour":
                if float(pora_threshold[split]["min"]) < temp_dict["plates"][methode]["wells"][wells] < \
                                            float(pora_threshold[split]["max"]):
                    temp_colour = pora_threshold["colour"][split]
                    # cell_color = colors.cnames[temp_colour]
                    cell_color = temp_colour.replace("#", "")
                    temp_cell = translate_wells_to_cells[wells]
                    ws[temp_cell].fill = PatternFill("solid", fgColor=cell_color)

    write_row = init_row + 1

    # Write the guide:
    for threshold in pora_threshold:
        ws[ex_cell(write_row, free_col)] = threshold
        ws[ex_cell(write_row, free_col)].font = Font(b=True, underline="single")
        for level in pora_threshold[threshold]:
            ws[ex_cell(write_row, free_col + 1)] = level
            ws[ex_cell(write_row, free_col + 2)] = pora_threshold[threshold][level]
            write_row += 1


def well_row_col_type(plate_layout):
    """
    Makes two dicts. one for what type/state each well/cell is and a dict that translate each well to a column and row
    values

    :param plate_layout: The plate-layout for what state each well is in (sample, min, max...)
    :type plate_layout: dict
    :return: well_col_row: a dict with a list of values for rows and for columns, well_type: a dict over the
        type/state of each well/cell
    :rtype: dict, dict
    """
    well_type = {}
    well_col_row = {"well_col": [], "well_row": []}
    if "well_layout" in plate_layout:
        temp_plate_layout = plate_layout["well_layout"]
    else:
        temp_plate_layout = plate_layout

    for counter in temp_plate_layout:
        for keys in temp_plate_layout[counter]:
            if keys == "well_id":
                temp_well_row = re.sub(r"\d+", "", temp_plate_layout[counter][keys])
                temp_well_col = re.sub(r"\D+", "", temp_plate_layout[counter][keys])
                if temp_well_row not in well_col_row["well_row"]:
                    well_col_row["well_row"].append(temp_well_row)
                if temp_well_col not in well_col_row["well_col"]:
                    well_col_row["well_col"].append(temp_well_col)

            if keys == "state":
                well_type.setdefault(temp_plate_layout[counter]["state"], []).append(
                    temp_plate_layout[counter]["well_id"])

    return well_col_row, well_type


def original_data_dict(file, plate_layout):
    """
    The original data from the plate reader, loaded in from the excel file

    :param file: the excel file, with the platereaders data
    :type file: str
    :param plate_layout: The platelayout to tell witch cell/well is in witch state (sample, blank, empty....)
    :type plate_layout: dict
    :return:
        - all_data: all_data: A dict over the original data, and what each cell/well is
        - well_col_row: a dict with a list of values for rows and for columns,
        - well_type: a dict over the type/state of each well/cell
        - barcode: The barcode for the plate
    :rtype:
        - dict
        - dict
        - dict
        - str
    """
    well_col_row, well_type = well_row_col_type(plate_layout)
    wb = load_workbook(file)        #Todo Make it work for txt file, as platebutle spits txt files out. I could Check for file type, convert txt to excel using test_script_2...
    sheet = wb.sheetnames[0]
    ws = wb[sheet]
    plate_type_384 = False
    plate_type_1536 = False

    all_data = {"plates": {}, "calculations": {}}
    n_rows = len(well_col_row["well_row"])
    for row_index, row in enumerate(ws.values):     #Todo Change this to be a dict, as there can be multiple readings per data.
        for index_row, value in enumerate(row):

            if value == "Date:":
                date = row[4]

            if value == "Name" and row[1]:
                barcode = row[1]

            if value == "<>":
                skipped_rows = row_index
            if value == "I":
                plate_type_384 = True
            if value == "AA":
                plate_type_1536 = True
            try:
                print(skipped_rows)
            except UnboundLocalError:
                pass

    df_plate = pd.read_excel(file, sheet_name=sheet, skiprows=skipped_rows, nrows=n_rows)

    df_plate_dict = df_plate.to_dict()
    if n_rows == 8 and plate_type_384:
        PySimpleGUI.PopupError("Wrong plate size, data is larger than 96-wells")
        return None, None, None
    elif n_rows == 16 and not plate_type_384:
        PySimpleGUI.PopupError("Wrong plate size, data is 96-wells")
        return None, None, None
    elif n_rows == 16 and plate_type_1536:
        PySimpleGUI.PopupError("wrong plate size, data is 1536-wells")
        return None, None, None
    elif n_rows == 32 and not plate_type_1536:
        PySimpleGUI.PopupError("wrong plate size, data is larger smaller than 1536-wells")
        return None, None, None

    temp_reading_dict = {}
    for counter, heading in enumerate(df_plate_dict):
        for index, values in enumerate(df_plate_dict[heading]):
            temp_reading_dict[f"{well_col_row['well_row'][index]}{counter}"] = df_plate_dict[heading][values]

    all_data["plates"]["original"] = {}
    all_data["plates"]["original"]["wells"] = {}
    for state in well_type:
        for well in well_type[state]:
            try:
                all_data["plates"]["original"]["wells"][well] = temp_reading_dict[well]

            # Incase blanked wells are skipped for the reading
            except KeyError:
                all_data["plates"]["original"]["wells"][well] = 1
    try:
        barcode
    except UnboundLocalError:
        barcode = file.split("/")[-1]


    try:
        date
    except UnboundLocalError:
        date = datetime.today()

    return all_data, well_col_row, well_type, barcode, date


if __name__ == "__main__":

    plate_layout = {"well_layout": {"1": {"group": 0, "well_id": "A1", "state": "empty", "colour": "#1e0bc8"}, "2": {"group": 0, "well_id": "B1", "state": "empty", "colour": "#1e0bc8"}, "3": {"group": 0, "well_id": "C1", "state": "empty", "colour": "#1e0bc8"}, "4": {"group": 0, "well_id": "D1", "state": "empty", "colour": "#1e0bc8"}, "5": {"group": 0, "well_id": "E1", "state": "empty", "colour": "#1e0bc8"}, "6": {"group": 0, "well_id": "F1", "state": "empty", "colour": "#1e0bc8"}, "7": {"group": 0, "well_id": "G1", "state": "empty", "colour": "#1e0bc8"}, "8": {"group": 0, "well_id": "H1", "state": "empty", "colour": "#1e0bc8"}, "9": {"group": 0, "well_id": "I1", "state": "empty", "colour": "#1e0bc8"}, "10": {"group": 0, "well_id": "J1", "state": "empty", "colour": "#1e0bc8"}, "11": {"group": 0, "well_id": "K1", "state": "empty", "colour": "#1e0bc8"}, "12": {"group": 0, "well_id": "L1", "state": "empty", "colour": "#1e0bc8"}, "13": {"group": 0, "well_id": "M1", "state": "empty", "colour": "#1e0bc8"}, "14": {"group": 0, "well_id": "N1", "state": "empty", "colour": "#1e0bc8"}, "15": {"group": 0, "well_id": "O1", "state": "empty", "colour": "#1e0bc8"}, "16": {"group": 0, "well_id": "P1", "state": "empty", "colour": "#1e0bc8"}, "17": {"group": 0, "well_id": "A2", "state": "empty", "colour": "#1e0bc8"}, "18": {"group": 0, "well_id": "B2", "state": "minimum", "colour": "#ff8000"}, "19": {"group": 0, "well_id": "C2", "state": "minimum", "colour": "#ff8000"}, "20": {"group": 0, "well_id": "D2", "state": "minimum", "colour": "#ff8000"}, "21": {"group": 0, "well_id": "E2", "state": "minimum", "colour": "#ff8000"}, "22": {"group": 0, "well_id": "F2", "state": "minimum", "colour": "#ff8000"}, "23": {"group": 0, "well_id": "G2", "state": "minimum", "colour": "#ff8000"}, "24": {"group": 0, "well_id": "H2", "state": "minimum", "colour": "#ff8000"}, "25": {"group": 0, "well_id": "I2", "state": "minimum", "colour": "#ff8000"}, "26": {"group": 0, "well_id": "J2", "state": "minimum", "colour": "#ff8000"}, "27": {"group": 0, "well_id": "K2", "state": "minimum", "colour": "#ff8000"}, "28": {"group": 0, "well_id": "L2", "state": "minimum", "colour": "#ff8000"}, "29": {"group": 0, "well_id": "M2", "state": "minimum", "colour": "#ff8000"}, "30": {"group": 0, "well_id": "N2", "state": "minimum", "colour": "#ff8000"}, "31": {"group": 0, "well_id": "O2", "state": "minimum", "colour": "#ff8000"}, "32": {"group": 0, "well_id": "P2", "state": "empty", "colour": "#1e0bc8"}, "33": {"group": 0, "well_id": "A3", "state": "empty", "colour": "#1e0bc8"}, "34": {"group": 0, "well_id": "B3", "state": "max", "colour": "#790dc1"}, "35": {"group": 0, "well_id": "C3", "state": "max", "colour": "#790dc1"}, "36": {"group": 0, "well_id": "D3", "state": "max", "colour": "#790dc1"}, "37": {"group": 0, "well_id": "E3", "state": "max", "colour": "#790dc1"}, "38": {"group": 0, "well_id": "F3", "state": "max", "colour": "#790dc1"}, "39": {"group": 0, "well_id": "G3", "state": "max", "colour": "#790dc1"}, "40": {"group": 0, "well_id": "H3", "state": "max", "colour": "#790dc1"}, "41": {"group": 0, "well_id": "I3", "state": "max", "colour": "#790dc1"}, "42": {"group": 0, "well_id": "J3", "state": "max", "colour": "#790dc1"}, "43": {"group": 0, "well_id": "K3", "state": "max", "colour": "#790dc1"}, "44": {"group": 0, "well_id": "L3", "state": "max", "colour": "#790dc1"}, "45": {"group": 0, "well_id": "M3", "state": "max", "colour": "#790dc1"}, "46": {"group": 0, "well_id": "N3", "state": "max", "colour": "#790dc1"}, "47": {"group": 0, "well_id": "O3", "state": "max", "colour": "#790dc1"}, "48": {"group": 0, "well_id": "P3", "state": "empty", "colour": "#1e0bc8"}, "49": {"group": 0, "well_id": "A4", "state": "empty", "colour": "#1e0bc8"}, "50": {"group": 0, "well_id": "B4", "state": "sample", "colour": "#ff00ff"}, "51": {"group": 0, "well_id": "C4", "state": "sample", "colour": "#ff00ff"}, "52": {"group": 0, "well_id": "D4", "state": "sample", "colour": "#ff00ff"}, "53": {"group": 0, "well_id": "E4", "state": "sample", "colour": "#ff00ff"}, "54": {"group": 0, "well_id": "F4", "state": "sample", "colour": "#ff00ff"}, "55": {"group": 0, "well_id": "G4", "state": "sample", "colour": "#ff00ff"}, "56": {"group": 0, "well_id": "H4", "state": "sample", "colour": "#ff00ff"}, "57": {"group": 0, "well_id": "I4", "state": "sample", "colour": "#ff00ff"}, "58": {"group": 0, "well_id": "J4", "state": "sample", "colour": "#ff00ff"}, "59": {"group": 0, "well_id": "K4", "state": "sample", "colour": "#ff00ff"}, "60": {"group": 0, "well_id": "L4", "state": "sample", "colour": "#ff00ff"}, "61": {"group": 0, "well_id": "M4", "state": "sample", "colour": "#ff00ff"}, "62": {"group": 0, "well_id": "N4", "state": "sample", "colour": "#ff00ff"}, "63": {"group": 0, "well_id": "O4", "state": "sample", "colour": "#ff00ff"}, "64": {"group": 0, "well_id": "P4", "state": "empty", "colour": "#1e0bc8"}, "65": {"group": 0, "well_id": "A5", "state": "empty", "colour": "#1e0bc8"}, "66": {"group": 0, "well_id": "B5", "state": "sample", "colour": "#ff00ff"}, "67": {"group": 0, "well_id": "C5", "state": "sample", "colour": "#ff00ff"}, "68": {"group": 0, "well_id": "D5", "state": "sample", "colour": "#ff00ff"}, "69": {"group": 0, "well_id": "E5", "state": "sample", "colour": "#ff00ff"}, "70": {"group": 0, "well_id": "F5", "state": "sample", "colour": "#ff00ff"}, "71": {"group": 0, "well_id": "G5", "state": "sample", "colour": "#ff00ff"}, "72": {"group": 0, "well_id": "H5", "state": "sample", "colour": "#ff00ff"}, "73": {"group": 0, "well_id": "I5", "state": "sample", "colour": "#ff00ff"}, "74": {"group": 0, "well_id": "J5", "state": "sample", "colour": "#ff00ff"}, "75": {"group": 0, "well_id": "K5", "state": "sample", "colour": "#ff00ff"}, "76": {"group": 0, "well_id": "L5", "state": "sample", "colour": "#ff00ff"}, "77": {"group": 0, "well_id": "M5", "state": "sample", "colour": "#ff00ff"}, "78": {"group": 0, "well_id": "N5", "state": "sample", "colour": "#ff00ff"}, "79": {"group": 0, "well_id": "O5", "state": "sample", "colour": "#ff00ff"}, "80": {"group": 0, "well_id": "P5", "state": "empty", "colour": "#1e0bc8"}, "81": {"group": 0, "well_id": "A6", "state": "empty", "colour": "#1e0bc8"}, "82": {"group": 0, "well_id": "B6", "state": "sample", "colour": "#ff00ff"}, "83": {"group": 0, "well_id": "C6", "state": "sample", "colour": "#ff00ff"}, "84": {"group": 0, "well_id": "D6", "state": "sample", "colour": "#ff00ff"}, "85": {"group": 0, "well_id": "E6", "state": "sample", "colour": "#ff00ff"}, "86": {"group": 0, "well_id": "F6", "state": "sample", "colour": "#ff00ff"}, "87": {"group": 0, "well_id": "G6", "state": "sample", "colour": "#ff00ff"}, "88": {"group": 0, "well_id": "H6", "state": "sample", "colour": "#ff00ff"}, "89": {"group": 0, "well_id": "I6", "state": "sample", "colour": "#ff00ff"}, "90": {"group": 0, "well_id": "J6", "state": "sample", "colour": "#ff00ff"}, "91": {"group": 0, "well_id": "K6", "state": "sample", "colour": "#ff00ff"}, "92": {"group": 0, "well_id": "L6", "state": "sample", "colour": "#ff00ff"}, "93": {"group": 0, "well_id": "M6", "state": "sample", "colour": "#ff00ff"}, "94": {"group": 0, "well_id": "N6", "state": "sample", "colour": "#ff00ff"}, "95": {"group": 0, "well_id": "O6", "state": "sample", "colour": "#ff00ff"}, "96": {"group": 0, "well_id": "P6", "state": "empty", "colour": "#1e0bc8"}, "97": {"group": 0, "well_id": "A7", "state": "empty", "colour": "#1e0bc8"}, "98": {"group": 0, "well_id": "B7", "state": "sample", "colour": "#ff00ff"}, "99": {"group": 0, "well_id": "C7", "state": "sample", "colour": "#ff00ff"}, "100": {"group": 0, "well_id": "D7", "state": "sample", "colour": "#ff00ff"}, "101": {"group": 0, "well_id": "E7", "state": "sample", "colour": "#ff00ff"}, "102": {"group": 0, "well_id": "F7", "state": "sample", "colour": "#ff00ff"}, "103": {"group": 0, "well_id": "G7", "state": "sample", "colour": "#ff00ff"}, "104": {"group": 0, "well_id": "H7", "state": "sample", "colour": "#ff00ff"}, "105": {"group": 0, "well_id": "I7", "state": "sample", "colour": "#ff00ff"}, "106": {"group": 0, "well_id": "J7", "state": "sample", "colour": "#ff00ff"}, "107": {"group": 0, "well_id": "K7", "state": "sample", "colour": "#ff00ff"}, "108": {"group": 0, "well_id": "L7", "state": "sample", "colour": "#ff00ff"}, "109": {"group": 0, "well_id": "M7", "state": "sample", "colour": "#ff00ff"}, "110": {"group": 0, "well_id": "N7", "state": "sample", "colour": "#ff00ff"}, "111": {"group": 0, "well_id": "O7", "state": "sample", "colour": "#ff00ff"}, "112": {"group": 0, "well_id": "P7", "state": "empty", "colour": "#1e0bc8"}, "113": {"group": 0, "well_id": "A8", "state": "empty", "colour": "#1e0bc8"}, "114": {"group": 0, "well_id": "B8", "state": "sample", "colour": "#ff00ff"}, "115": {"group": 0, "well_id": "C8", "state": "sample", "colour": "#ff00ff"}, "116": {"group": 0, "well_id": "D8", "state": "sample", "colour": "#ff00ff"}, "117": {"group": 0, "well_id": "E8", "state": "sample", "colour": "#ff00ff"}, "118": {"group": 0, "well_id": "F8", "state": "sample", "colour": "#ff00ff"}, "119": {"group": 0, "well_id": "G8", "state": "sample", "colour": "#ff00ff"}, "120": {"group": 0, "well_id": "H8", "state": "sample", "colour": "#ff00ff"}, "121": {"group": 0, "well_id": "I8", "state": "sample", "colour": "#ff00ff"}, "122": {"group": 0, "well_id": "J8", "state": "sample", "colour": "#ff00ff"}, "123": {"group": 0, "well_id": "K8", "state": "sample", "colour": "#ff00ff"}, "124": {"group": 0, "well_id": "L8", "state": "sample", "colour": "#ff00ff"}, "125": {"group": 0, "well_id": "M8", "state": "sample", "colour": "#ff00ff"}, "126": {"group": 0, "well_id": "N8", "state": "sample", "colour": "#ff00ff"}, "127": {"group": 0, "well_id": "O8", "state": "sample", "colour": "#ff00ff"}, "128": {"group": 0, "well_id": "P8", "state": "empty", "colour": "#1e0bc8"}, "129": {"group": 0, "well_id": "A9", "state": "empty", "colour": "#1e0bc8"}, "130": {"group": 0, "well_id": "B9", "state": "sample", "colour": "#ff00ff"}, "131": {"group": 0, "well_id": "C9", "state": "sample", "colour": "#ff00ff"}, "132": {"group": 0, "well_id": "D9", "state": "sample", "colour": "#ff00ff"}, "133": {"group": 0, "well_id": "E9", "state": "sample", "colour": "#ff00ff"}, "134": {"group": 0, "well_id": "F9", "state": "sample", "colour": "#ff00ff"}, "135": {"group": 0, "well_id": "G9", "state": "sample", "colour": "#ff00ff"}, "136": {"group": 0, "well_id": "H9", "state": "sample", "colour": "#ff00ff"}, "137": {"group": 0, "well_id": "I9", "state": "sample", "colour": "#ff00ff"}, "138": {"group": 0, "well_id": "J9", "state": "sample", "colour": "#ff00ff"}, "139": {"group": 0, "well_id": "K9", "state": "sample", "colour": "#ff00ff"}, "140": {"group": 0, "well_id": "L9", "state": "sample", "colour": "#ff00ff"}, "141": {"group": 0, "well_id": "M9", "state": "sample", "colour": "#ff00ff"}, "142": {"group": 0, "well_id": "N9", "state": "sample", "colour": "#ff00ff"}, "143": {"group": 0, "well_id": "O9", "state": "sample", "colour": "#ff00ff"}, "144": {"group": 0, "well_id": "P9", "state": "empty", "colour": "#1e0bc8"}, "145": {"group": 0, "well_id": "A10", "state": "empty", "colour": "#1e0bc8"}, "146": {"group": 0, "well_id": "B10", "state": "sample", "colour": "#ff00ff"}, "147": {"group": 0, "well_id": "C10", "state": "sample", "colour": "#ff00ff"}, "148": {"group": 0, "well_id": "D10", "state": "sample", "colour": "#ff00ff"}, "149": {"group": 0, "well_id": "E10", "state": "sample", "colour": "#ff00ff"}, "150": {"group": 0, "well_id": "F10", "state": "sample", "colour": "#ff00ff"}, "151": {"group": 0, "well_id": "G10", "state": "sample", "colour": "#ff00ff"}, "152": {"group": 0, "well_id": "H10", "state": "sample", "colour": "#ff00ff"}, "153": {"group": 0, "well_id": "I10", "state": "sample", "colour": "#ff00ff"}, "154": {"group": 0, "well_id": "J10", "state": "sample", "colour": "#ff00ff"}, "155": {"group": 0, "well_id": "K10", "state": "sample", "colour": "#ff00ff"}, "156": {"group": 0, "well_id": "L10", "state": "sample", "colour": "#ff00ff"}, "157": {"group": 0, "well_id": "M10", "state": "sample", "colour": "#ff00ff"}, "158": {"group": 0, "well_id": "N10", "state": "sample", "colour": "#ff00ff"}, "159": {"group": 0, "well_id": "O10", "state": "sample", "colour": "#ff00ff"}, "160": {"group": 0, "well_id": "P10", "state": "empty", "colour": "#1e0bc8"}, "161": {"group": 0, "well_id": "A11", "state": "empty", "colour": "#1e0bc8"}, "162": {"group": 0, "well_id": "B11", "state": "sample", "colour": "#ff00ff"}, "163": {"group": 0, "well_id": "C11", "state": "sample", "colour": "#ff00ff"}, "164": {"group": 0, "well_id": "D11", "state": "sample", "colour": "#ff00ff"}, "165": {"group": 0, "well_id": "E11", "state": "sample", "colour": "#ff00ff"}, "166": {"group": 0, "well_id": "F11", "state": "sample", "colour": "#ff00ff"}, "167": {"group": 0, "well_id": "G11", "state": "sample", "colour": "#ff00ff"}, "168": {"group": 0, "well_id": "H11", "state": "sample", "colour": "#ff00ff"}, "169": {"group": 0, "well_id": "I11", "state": "sample", "colour": "#ff00ff"}, "170": {"group": 0, "well_id": "J11", "state": "sample", "colour": "#ff00ff"}, "171": {"group": 0, "well_id": "K11", "state": "sample", "colour": "#ff00ff"}, "172": {"group": 0, "well_id": "L11", "state": "sample", "colour": "#ff00ff"}, "173": {"group": 0, "well_id": "M11", "state": "sample", "colour": "#ff00ff"}, "174": {"group": 0, "well_id": "N11", "state": "sample", "colour": "#ff00ff"}, "175": {"group": 0, "well_id": "O11", "state": "sample", "colour": "#ff00ff"}, "176": {"group": 0, "well_id": "P11", "state": "empty", "colour": "#1e0bc8"}, "177": {"group": 0, "well_id": "A12", "state": "empty", "colour": "#1e0bc8"}, "178": {"group": 0, "well_id": "B12", "state": "sample", "colour": "#ff00ff"}, "179": {"group": 0, "well_id": "C12", "state": "sample", "colour": "#ff00ff"}, "180": {"group": 0, "well_id": "D12", "state": "sample", "colour": "#ff00ff"}, "181": {"group": 0, "well_id": "E12", "state": "sample", "colour": "#ff00ff"}, "182": {"group": 0, "well_id": "F12", "state": "sample", "colour": "#ff00ff"}, "183": {"group": 0, "well_id": "G12", "state": "sample", "colour": "#ff00ff"}, "184": {"group": 0, "well_id": "H12", "state": "sample", "colour": "#ff00ff"}, "185": {"group": 0, "well_id": "I12", "state": "sample", "colour": "#ff00ff"}, "186": {"group": 0, "well_id": "J12", "state": "sample", "colour": "#ff00ff"}, "187": {"group": 0, "well_id": "K12", "state": "sample", "colour": "#ff00ff"}, "188": {"group": 0, "well_id": "L12", "state": "sample", "colour": "#ff00ff"}, "189": {"group": 0, "well_id": "M12", "state": "sample", "colour": "#ff00ff"}, "190": {"group": 0, "well_id": "N12", "state": "sample", "colour": "#ff00ff"}, "191": {"group": 0, "well_id": "O12", "state": "sample", "colour": "#ff00ff"}, "192": {"group": 0, "well_id": "P12", "state": "empty", "colour": "#1e0bc8"}, "193": {"group": 0, "well_id": "A13", "state": "empty", "colour": "#1e0bc8"}, "194": {"group": 0, "well_id": "B13", "state": "sample", "colour": "#ff00ff"}, "195": {"group": 0, "well_id": "C13", "state": "sample", "colour": "#ff00ff"}, "196": {"group": 0, "well_id": "D13", "state": "sample", "colour": "#ff00ff"}, "197": {"group": 0, "well_id": "E13", "state": "sample", "colour": "#ff00ff"}, "198": {"group": 0, "well_id": "F13", "state": "sample", "colour": "#ff00ff"}, "199": {"group": 0, "well_id": "G13", "state": "sample", "colour": "#ff00ff"}, "200": {"group": 0, "well_id": "H13", "state": "sample", "colour": "#ff00ff"}, "201": {"group": 0, "well_id": "I13", "state": "sample", "colour": "#ff00ff"}, "202": {"group": 0, "well_id": "J13", "state": "sample", "colour": "#ff00ff"}, "203": {"group": 0, "well_id": "K13", "state": "sample", "colour": "#ff00ff"}, "204": {"group": 0, "well_id": "L13", "state": "sample", "colour": "#ff00ff"}, "205": {"group": 0, "well_id": "M13", "state": "sample", "colour": "#ff00ff"}, "206": {"group": 0, "well_id": "N13", "state": "sample", "colour": "#ff00ff"}, "207": {"group": 0, "well_id": "O13", "state": "sample", "colour": "#ff00ff"}, "208": {"group": 0, "well_id": "P13", "state": "empty", "colour": "#1e0bc8"}, "209": {"group": 0, "well_id": "A14", "state": "empty", "colour": "#1e0bc8"}, "210": {"group": 0, "well_id": "B14", "state": "sample", "colour": "#ff00ff"}, "211": {"group": 0, "well_id": "C14", "state": "sample", "colour": "#ff00ff"}, "212": {"group": 0, "well_id": "D14", "state": "sample", "colour": "#ff00ff"}, "213": {"group": 0, "well_id": "E14", "state": "sample", "colour": "#ff00ff"}, "214": {"group": 0, "well_id": "F14", "state": "sample", "colour": "#ff00ff"}, "215": {"group": 0, "well_id": "G14", "state": "sample", "colour": "#ff00ff"}, "216": {"group": 0, "well_id": "H14", "state": "sample", "colour": "#ff00ff"}, "217": {"group": 0, "well_id": "I14", "state": "sample", "colour": "#ff00ff"}, "218": {"group": 0, "well_id": "J14", "state": "sample", "colour": "#ff00ff"}, "219": {"group": 0, "well_id": "K14", "state": "sample", "colour": "#ff00ff"}, "220": {"group": 0, "well_id": "L14", "state": "sample", "colour": "#ff00ff"}, "221": {"group": 0, "well_id": "M14", "state": "sample", "colour": "#ff00ff"}, "222": {"group": 0, "well_id": "N14", "state": "sample", "colour": "#ff00ff"}, "223": {"group": 0, "well_id": "O14", "state": "sample", "colour": "#ff00ff"}, "224": {"group": 0, "well_id": "P14", "state": "empty", "colour": "#1e0bc8"}, "225": {"group": 0, "well_id": "A15", "state": "empty", "colour": "#1e0bc8"}, "226": {"group": 0, "well_id": "B15", "state": "sample", "colour": "#ff00ff"}, "227": {"group": 0, "well_id": "C15", "state": "sample", "colour": "#ff00ff"}, "228": {"group": 0, "well_id": "D15", "state": "sample", "colour": "#ff00ff"}, "229": {"group": 0, "well_id": "E15", "state": "sample", "colour": "#ff00ff"}, "230": {"group": 0, "well_id": "F15", "state": "sample", "colour": "#ff00ff"}, "231": {"group": 0, "well_id": "G15", "state": "sample", "colour": "#ff00ff"}, "232": {"group": 0, "well_id": "H15", "state": "sample", "colour": "#ff00ff"}, "233": {"group": 0, "well_id": "I15", "state": "sample", "colour": "#ff00ff"}, "234": {"group": 0, "well_id": "J15", "state": "sample", "colour": "#ff00ff"}, "235": {"group": 0, "well_id": "K15", "state": "sample", "colour": "#ff00ff"}, "236": {"group": 0, "well_id": "L15", "state": "sample", "colour": "#ff00ff"}, "237": {"group": 0, "well_id": "M15", "state": "sample", "colour": "#ff00ff"}, "238": {"group": 0, "well_id": "N15", "state": "sample", "colour": "#ff00ff"}, "239": {"group": 0, "well_id": "O15", "state": "sample", "colour": "#ff00ff"}, "240": {"group": 0, "well_id": "P15", "state": "empty", "colour": "#1e0bc8"}, "241": {"group": 0, "well_id": "A16", "state": "empty", "colour": "#1e0bc8"}, "242": {"group": 0, "well_id": "B16", "state": "sample", "colour": "#ff00ff"}, "243": {"group": 0, "well_id": "C16", "state": "sample", "colour": "#ff00ff"}, "244": {"group": 0, "well_id": "D16", "state": "sample", "colour": "#ff00ff"}, "245": {"group": 0, "well_id": "E16", "state": "sample", "colour": "#ff00ff"}, "246": {"group": 0, "well_id": "F16", "state": "sample", "colour": "#ff00ff"}, "247": {"group": 0, "well_id": "G16", "state": "sample", "colour": "#ff00ff"}, "248": {"group": 0, "well_id": "H16", "state": "sample", "colour": "#ff00ff"}, "249": {"group": 0, "well_id": "I16", "state": "sample", "colour": "#ff00ff"}, "250": {"group": 0, "well_id": "J16", "state": "sample", "colour": "#ff00ff"}, "251": {"group": 0, "well_id": "K16", "state": "sample", "colour": "#ff00ff"}, "252": {"group": 0, "well_id": "L16", "state": "sample", "colour": "#ff00ff"}, "253": {"group": 0, "well_id": "M16", "state": "sample", "colour": "#ff00ff"}, "254": {"group": 0, "well_id": "N16", "state": "sample", "colour": "#ff00ff"}, "255": {"group": 0, "well_id": "O16", "state": "sample", "colour": "#ff00ff"}, "256": {"group": 0, "well_id": "P16", "state": "empty", "colour": "#1e0bc8"}, "257": {"group": 0, "well_id": "A17", "state": "empty", "colour": "#1e0bc8"}, "258": {"group": 0, "well_id": "B17", "state": "sample", "colour": "#ff00ff"}, "259": {"group": 0, "well_id": "C17", "state": "sample", "colour": "#ff00ff"}, "260": {"group": 0, "well_id": "D17", "state": "sample", "colour": "#ff00ff"}, "261": {"group": 0, "well_id": "E17", "state": "sample", "colour": "#ff00ff"}, "262": {"group": 0, "well_id": "F17", "state": "sample", "colour": "#ff00ff"}, "263": {"group": 0, "well_id": "G17", "state": "sample", "colour": "#ff00ff"}, "264": {"group": 0, "well_id": "H17", "state": "sample", "colour": "#ff00ff"}, "265": {"group": 0, "well_id": "I17", "state": "sample", "colour": "#ff00ff"}, "266": {"group": 0, "well_id": "J17", "state": "sample", "colour": "#ff00ff"}, "267": {"group": 0, "well_id": "K17", "state": "sample", "colour": "#ff00ff"}, "268": {"group": 0, "well_id": "L17", "state": "sample", "colour": "#ff00ff"}, "269": {"group": 0, "well_id": "M17", "state": "sample", "colour": "#ff00ff"}, "270": {"group": 0, "well_id": "N17", "state": "sample", "colour": "#ff00ff"}, "271": {"group": 0, "well_id": "O17", "state": "sample", "colour": "#ff00ff"}, "272": {"group": 0, "well_id": "P17", "state": "empty", "colour": "#1e0bc8"}, "273": {"group": 0, "well_id": "A18", "state": "empty", "colour": "#1e0bc8"}, "274": {"group": 0, "well_id": "B18", "state": "sample", "colour": "#ff00ff"}, "275": {"group": 0, "well_id": "C18", "state": "sample", "colour": "#ff00ff"}, "276": {"group": 0, "well_id": "D18", "state": "sample", "colour": "#ff00ff"}, "277": {"group": 0, "well_id": "E18", "state": "sample", "colour": "#ff00ff"}, "278": {"group": 0, "well_id": "F18", "state": "sample", "colour": "#ff00ff"}, "279": {"group": 0, "well_id": "G18", "state": "sample", "colour": "#ff00ff"}, "280": {"group": 0, "well_id": "H18", "state": "sample", "colour": "#ff00ff"}, "281": {"group": 0, "well_id": "I18", "state": "sample", "colour": "#ff00ff"}, "282": {"group": 0, "well_id": "J18", "state": "sample", "colour": "#ff00ff"}, "283": {"group": 0, "well_id": "K18", "state": "sample", "colour": "#ff00ff"}, "284": {"group": 0, "well_id": "L18", "state": "sample", "colour": "#ff00ff"}, "285": {"group": 0, "well_id": "M18", "state": "sample", "colour": "#ff00ff"}, "286": {"group": 0, "well_id": "N18", "state": "sample", "colour": "#ff00ff"}, "287": {"group": 0, "well_id": "O18", "state": "sample", "colour": "#ff00ff"}, "288": {"group": 0, "well_id": "P18", "state": "empty", "colour": "#1e0bc8"}, "289": {"group": 0, "well_id": "A19", "state": "empty", "colour": "#1e0bc8"}, "290": {"group": 0, "well_id": "B19", "state": "sample", "colour": "#ff00ff"}, "291": {"group": 0, "well_id": "C19", "state": "sample", "colour": "#ff00ff"}, "292": {"group": 0, "well_id": "D19", "state": "sample", "colour": "#ff00ff"}, "293": {"group": 0, "well_id": "E19", "state": "sample", "colour": "#ff00ff"}, "294": {"group": 0, "well_id": "F19", "state": "sample", "colour": "#ff00ff"}, "295": {"group": 0, "well_id": "G19", "state": "sample", "colour": "#ff00ff"}, "296": {"group": 0, "well_id": "H19", "state": "sample", "colour": "#ff00ff"}, "297": {"group": 0, "well_id": "I19", "state": "sample", "colour": "#ff00ff"}, "298": {"group": 0, "well_id": "J19", "state": "sample", "colour": "#ff00ff"}, "299": {"group": 0, "well_id": "K19", "state": "sample", "colour": "#ff00ff"}, "300": {"group": 0, "well_id": "L19", "state": "sample", "colour": "#ff00ff"}, "301": {"group": 0, "well_id": "M19", "state": "sample", "colour": "#ff00ff"}, "302": {"group": 0, "well_id": "N19", "state": "sample", "colour": "#ff00ff"}, "303": {"group": 0, "well_id": "O19", "state": "sample", "colour": "#ff00ff"}, "304": {"group": 0, "well_id": "P19", "state": "empty", "colour": "#1e0bc8"}, "305": {"group": 0, "well_id": "A20", "state": "empty", "colour": "#1e0bc8"}, "306": {"group": 0, "well_id": "B20", "state": "sample", "colour": "#ff00ff"}, "307": {"group": 0, "well_id": "C20", "state": "sample", "colour": "#ff00ff"}, "308": {"group": 0, "well_id": "D20", "state": "sample", "colour": "#ff00ff"}, "309": {"group": 0, "well_id": "E20", "state": "sample", "colour": "#ff00ff"}, "310": {"group": 0, "well_id": "F20", "state": "sample", "colour": "#ff00ff"}, "311": {"group": 0, "well_id": "G20", "state": "sample", "colour": "#ff00ff"}, "312": {"group": 0, "well_id": "H20", "state": "sample", "colour": "#ff00ff"}, "313": {"group": 0, "well_id": "I20", "state": "sample", "colour": "#ff00ff"}, "314": {"group": 0, "well_id": "J20", "state": "sample", "colour": "#ff00ff"}, "315": {"group": 0, "well_id": "K20", "state": "sample", "colour": "#ff00ff"}, "316": {"group": 0, "well_id": "L20", "state": "sample", "colour": "#ff00ff"}, "317": {"group": 0, "well_id": "M20", "state": "sample", "colour": "#ff00ff"}, "318": {"group": 0, "well_id": "N20", "state": "sample", "colour": "#ff00ff"}, "319": {"group": 0, "well_id": "O20", "state": "sample", "colour": "#ff00ff"}, "320": {"group": 0, "well_id": "P20", "state": "empty", "colour": "#1e0bc8"}, "321": {"group": 0, "well_id": "A21", "state": "empty", "colour": "#1e0bc8"}, "322": {"group": 0, "well_id": "B21", "state": "sample", "colour": "#ff00ff"}, "323": {"group": 0, "well_id": "C21", "state": "sample", "colour": "#ff00ff"}, "324": {"group": 0, "well_id": "D21", "state": "sample", "colour": "#ff00ff"}, "325": {"group": 0, "well_id": "E21", "state": "sample", "colour": "#ff00ff"}, "326": {"group": 0, "well_id": "F21", "state": "sample", "colour": "#ff00ff"}, "327": {"group": 0, "well_id": "G21", "state": "sample", "colour": "#ff00ff"}, "328": {"group": 0, "well_id": "H21", "state": "sample", "colour": "#ff00ff"}, "329": {"group": 0, "well_id": "I21", "state": "sample", "colour": "#ff00ff"}, "330": {"group": 0, "well_id": "J21", "state": "sample", "colour": "#ff00ff"}, "331": {"group": 0, "well_id": "K21", "state": "sample", "colour": "#ff00ff"}, "332": {"group": 0, "well_id": "L21", "state": "sample", "colour": "#ff00ff"}, "333": {"group": 0, "well_id": "M21", "state": "sample", "colour": "#ff00ff"}, "334": {"group": 0, "well_id": "N21", "state": "sample", "colour": "#ff00ff"}, "335": {"group": 0, "well_id": "O21", "state": "sample", "colour": "#ff00ff"}, "336": {"group": 0, "well_id": "P21", "state": "empty", "colour": "#1e0bc8"}, "337": {"group": 0, "well_id": "A22", "state": "empty", "colour": "#1e0bc8"}, "338": {"group": 0, "well_id": "B22", "state": "sample", "colour": "#ff00ff"}, "339": {"group": 0, "well_id": "C22", "state": "sample", "colour": "#ff00ff"}, "340": {"group": 0, "well_id": "D22", "state": "sample", "colour": "#ff00ff"}, "341": {"group": 0, "well_id": "E22", "state": "sample", "colour": "#ff00ff"}, "342": {"group": 0, "well_id": "F22", "state": "sample", "colour": "#ff00ff"}, "343": {"group": 0, "well_id": "G22", "state": "sample", "colour": "#ff00ff"}, "344": {"group": 0, "well_id": "H22", "state": "sample", "colour": "#ff00ff"}, "345": {"group": 0, "well_id": "I22", "state": "sample", "colour": "#ff00ff"}, "346": {"group": 0, "well_id": "J22", "state": "sample", "colour": "#ff00ff"}, "347": {"group": 0, "well_id": "K22", "state": "sample", "colour": "#ff00ff"}, "348": {"group": 0, "well_id": "L22", "state": "sample", "colour": "#ff00ff"}, "349": {"group": 0, "well_id": "M22", "state": "sample", "colour": "#ff00ff"}, "350": {"group": 0, "well_id": "N22", "state": "sample", "colour": "#ff00ff"}, "351": {"group": 0, "well_id": "O22", "state": "sample", "colour": "#ff00ff"}, "352": {"group": 0, "well_id": "P22", "state": "empty", "colour": "#1e0bc8"}, "353": {"group": 0, "well_id": "A23", "state": "empty", "colour": "#1e0bc8"}, "354": {"group": 0, "well_id": "B23", "state": "empty", "colour": "#1e0bc8"}, "355": {"group": 0, "well_id": "C23", "state": "empty", "colour": "#1e0bc8"}, "356": {"group": 0, "well_id": "D23", "state": "empty", "colour": "#1e0bc8"}, "357": {"group": 0, "well_id": "E23", "state": "empty", "colour": "#1e0bc8"}, "358": {"group": 0, "well_id": "F23", "state": "empty", "colour": "#1e0bc8"}, "359": {"group": 0, "well_id": "G23", "state": "empty", "colour": "#1e0bc8"}, "360": {"group": 0, "well_id": "H23", "state": "empty", "colour": "#1e0bc8"}, "361": {"group": 0, "well_id": "I23", "state": "empty", "colour": "#1e0bc8"}, "362": {"group": 0, "well_id": "J23", "state": "empty", "colour": "#1e0bc8"}, "363": {"group": 0, "well_id": "K23", "state": "empty", "colour": "#1e0bc8"}, "364": {"group": 0, "well_id": "L23", "state": "empty", "colour": "#1e0bc8"}, "365": {"group": 0, "well_id": "M23", "state": "empty", "colour": "#1e0bc8"}, "366": {"group": 0, "well_id": "N23", "state": "empty", "colour": "#1e0bc8"}, "367": {"group": 0, "well_id": "O23", "state": "empty", "colour": "#1e0bc8"}, "368": {"group": 0, "well_id": "P23", "state": "empty", "colour": "#1e0bc8"}, "369": {"group": 0, "well_id": "A24", "state": "empty", "colour": "#1e0bc8"}, "370": {"group": 0, "well_id": "B24", "state": "empty", "colour": "#1e0bc8"}, "371": {"group": 0, "well_id": "C24", "state": "empty", "colour": "#1e0bc8"}, "372": {"group": 0, "well_id": "D24", "state": "empty", "colour": "#1e0bc8"}, "373": {"group": 0, "well_id": "E24", "state": "empty", "colour": "#1e0bc8"}, "374": {"group": 0, "well_id": "F24", "state": "empty", "colour": "#1e0bc8"}, "375": {"group": 0, "well_id": "G24", "state": "empty", "colour": "#1e0bc8"}, "376": {"group": 0, "well_id": "H24", "state": "empty", "colour": "#1e0bc8"}, "377": {"group": 0, "well_id": "I24", "state": "empty", "colour": "#1e0bc8"}, "378": {"group": 0, "well_id": "J24", "state": "empty", "colour": "#1e0bc8"}, "379": {"group": 0, "well_id": "K24", "state": "empty", "colour": "#1e0bc8"}, "380": {"group": 0, "well_id": "L24", "state": "empty", "colour": "#1e0bc8"}, "381": {"group": 0, "well_id": "M24", "state": "empty", "colour": "#1e0bc8"}, "382": {"group": 0, "well_id": "N24", "state": "empty", "colour": "#1e0bc8"}, "383": {"group": 0, "well_id": "O24", "state": "empty", "colour": "#1e0bc8"}, "384": {"group": 0, "well_id": "P24", "state": "empty", "colour": "#1e0bc8"}}, "plate_type": "plate_384"}
    folder ="C:/Users/phch/Desktop/Daniels data/alphaSO-ELISA_20221216_151031/Export/xlsx"
    file = "alphaSO-ELISA_20221216_151223-plate-2.xlsx"
    full_file = f"{folder}/{file}"
    original_data_dict(full_file, plate_layout)