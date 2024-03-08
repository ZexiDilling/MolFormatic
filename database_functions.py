from database_controller import AddData, FetchData
from database_handler import DataBaseFunctions


def get_number_of_rows(config, table):
    dbf = DataBaseFunctions(config)
    return dbf.number_of_rows(table)


def update_database(config, table, data, file_type=None):
    """
    Update the database with data. Both for adding data to the different tables, but updating tables with new values

    :param data: Output file from the plate_butler system, tube files for the comPOUND freezer, or other data.
    :type data: str or dict
    :param file_type: If the file is a CSV file, it needs a specific file_type to know how to handle the file. as the
        CSV files are different depending on where they are coming from.
    :type file_type: str
    :param table: What table needs to be updated/where the data is going
    :type table: str
    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :return: Updated database with values
    """

    ad_db = AddData(config)
    ad_db.add_controller(table, data, file_type)


def grab_compound_table_data(config, table, sample):
    dbf = DataBaseFunctions(config)
    Search_limiter = {table:
                          {"value": [sample],
                           "operator": "IN",
                           "target_column": "compound_id",
                           "use": True}
                      }

    return dbf.return_table_data(table, Search_limiter)


def database_to_table(config, table, headings):
    dbf = DataBaseFunctions(config)
    row_data = dbf.return_table_data(table, None)

    table_data = []
    for row in row_data:
        temp_table_row = []
        for data in row_data[row]:
            if data in headings:
                temp_table_row.append(row_data[row][data])
        table_data.append(temp_table_row)

    return table_data


def delete_records_from_database(config, table, headline, data_value):
    """
    Deletes a record from the database
    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :param table: What table are the plates in
    :type table: str
    :param data_value: The value of the thing you are looking for
    :type data_value: str
    :param headline: Headline for the column where  the data is, in the table
    :type headline: str
    :return:
    """
    dbf = DataBaseFunctions(config)
    dbf.delete_records(table, headline, data_value)


def rename_record_in_the_database(config, table, headline_for_changing_value, headline_for_indicator_value, indicator_value,
                                  new_value):
    """
    Rename a record from the database
   :param table: The table where the data is located
    :type table: str
    :param headline_for_changing_value: the column headline for the data that needs to be changed
    :type headline_for_changing_value: str
    :param headline_for_indicator_value: The headling for the indicator value
    :type headline_for_indicator_value: str
    :param indicator_value: A value to find the right row from
    :type indicator_value: str
    :param new_value: The new value, that eeds to be changed to
    :type new_value: str
    :return:
    """
    dbf = DataBaseFunctions(config)
    dbf.rename_record_value( table, headline_for_changing_value, headline_for_indicator_value, indicator_value,
                            new_value)


def grab_table_data(config, table_name, search_limiter=None, single_row=None, data_value=None, headline=None,
                    search_list_clm=None, specific_rows=None):
    """
    Grabs data from tables
    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :param table_name: What table are the plates in
    :type table_name: str
    :param search_limiter: A dict to limit the search
    :type search_limiter: dict or list or None
    :param single_row: bool statement if it is a sling row look up
    :type single_row: bool
    :param data_value: The value of the thing you are looking for
    :type data_value: any or None
    :param headline: Headline for the column where  the data is, in the table
    :type headline: str or None
    :return:
    """

    if single_row:
        dbf = DataBaseFunctions(config)
        row_data = dbf.find_data_single_lookup(table_name, data_value, headline)
        return row_data
    else:
        fd = FetchData(config)
        rows = fd.data_search(table_name, search_limiter, search_list_clm, specific_rows)
        all_table_data = []
        headlines = []
        if rows:
            for row in rows:
                temp_data = []
                try:
                    rows[row]
                except TypeError:
                    print(f"table_name: {table_name}")
                    print(f"search_limiter: {search_limiter}")
                    print(f"search_list_clm: {search_list_clm}")
                    print(f"specific_rows: {specific_rows}")
                else:
                    for data in rows[row]:
                        headlines.append(data)
                        temp_data.append(rows[row][data])

                all_table_data.append(temp_data)

            return all_table_data, headlines
        else:
            return None, None


def _get_list_of_names_from_database_double_lookup(dbf, table, column_headline, limiting_value, limiting_header):

    temp_rows = dbf.find_data_single_lookup(table, limiting_value, limiting_header)
    table_header = dbf.grab_table_headers(table)
    try:
        table_index = table_header.index(column_headline)
    except ValueError:
        return None
    else:
        name_list = []
        for row_data in temp_rows:
            for counter, data in enumerate(row_data):
                if counter == table_index + 1:
                    name_list.append(data)

    return name_list


def _get_list_of_names_from_database(dbf, table, column_headline):
    """
    Gets a list of the plate names from the database
    :param dbf: The DataBaseFunction
    :type dbf: class
    :param table: The table where the data is storage
    :type table: str
    :param column_headline: The headline of the column that is used for finding data
    :type column_headline: str
    :return:
    """
    temp_rows = dbf.find_column_data(table, column_headline)
    temp_list = [names for names in temp_rows]
    return temp_list


def _get_plate_archive(dbf, table, column_headline, plate_names):
    """
    Gets a dict of the plate layouts from the database, based on a list of plate names
    :param dbf: The DataBaseFunction
    :type dbf: class
    :param plate_names: A list of plate names
    :type plate_names: list
    :param table: The table where the data is storage
    :type table: str
    :param column_headline: The headline of the column that is used for finding data
    :type column_headline: str
    :return: a dict of the plate layouts
    :rtype: dict
    """
    archive_plates_dict = {}

    # loops over all the plates
    for plates in plate_names:
        temp_row_data = dbf.find_data_single_lookup(table, plates, column_headline)

        # grabs data from the return rown
        plate_name = temp_row_data[0][1]
        plate_type = temp_row_data[0][2]

        temp_sub_row_data = dbf.find_data_single_lookup("plate_layout", plate_name, "layout_name")
        try:
            well_layout = eval(temp_sub_row_data[0][5])
        except IndexError:
            well_layout = None

        try:
            sample_type = temp_sub_row_data[0][4]
        except IndexError:
            sample_type = None

        # Generates the dict
        archive_plates_dict[plate_name] = {"well_layout": well_layout,
                                           "plate_type": plate_type,
                                           "sample_type": sample_type,
                                           "sample": [],
                                           "blank": [],
                                           "max": [],
                                           "minimum": [],
                                           "positive": [],
                                           "negative": [],
                                           "empty": []}

        # Makes list of the different well-types and adds them
        if well_layout is not None:
            temp_well_dict = well_layout
            for counter in temp_well_dict:
                well_id = temp_well_dict[counter]["well_id"]
                if temp_well_dict[counter]["state"] == "sample":
                    archive_plates_dict[plate_name]["sample"].append(well_id)
                elif temp_well_dict[counter]["state"] == "blank":
                    archive_plates_dict[plate_name]["blank"].append(well_id)
                elif temp_well_dict[counter]["state"] == "max":
                    archive_plates_dict[plate_name]["max"].append(well_id)
                elif temp_well_dict[counter]["state"] == "minimum":
                    archive_plates_dict[plate_name]["minimum"].append(well_id)
                elif temp_well_dict[counter]["state"] == "positive":
                    archive_plates_dict[plate_name]["positive"].append(well_id)
                elif temp_well_dict[counter]["state"] == "negative":
                    archive_plates_dict[plate_name]["negative"].append(well_id)
                elif temp_well_dict[counter]["state"] == "empty":
                    archive_plates_dict[plate_name]["empty"].append(well_id)

    return archive_plates_dict


def get_plate_layout(config):
    """
    Gets a list of plate names and their layout from the database
    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :return: plate_names, archive_plates_dict
    :rtype: list, dict
    """

    # Connects to the database and setting up standard values
    dbf = DataBaseFunctions(config)
    table = "plate_layout"
    column_headline = "layout_name"

    # gets a list of the plate names
    plate_names = _get_list_of_names_from_database(dbf, table, column_headline)

    # Gets a dict over all the plate_layouts
    archive_plates_dict = _get_plate_archive(dbf, table, column_headline, plate_names)

    return plate_names, archive_plates_dict


def _get_assay_list(config):
    assay_table_data = grab_table_data(config, "assay")
    assay_list = []
    for row in assay_table_data[0]:
        assay_list.append(row[1])

    return assay_list

#
# def gui_data_fetcher(db_active, config):
#     if db_active:
#         plate_list, archive_plates_dict = get_plate_layout(config)
#         assay_list = _get_assay_list(config)
#     else:
#         archive_plates_dict = {}
#         plate_list = []
#         assay_list = []
#
#     return plate_list, assay_list, archive_plates_dict

