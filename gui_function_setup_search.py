from chem_operators import structure_search
from database_controller import FetchData
import PySimpleGUI as sg

from gui_popup import export_chooser_popup

def search_compound(window, values, config):
    ac = values["-SEARCH_AC-"]
    origin_values = []
    if ac:
        for acs in ac:
            acs = acs.casefold()
            for values in config[f"database_specific_{acs}"]:
                origin_values.append(config[f"database_specific_{acs}"][values])
    window["-SEARCH_ORIGIN-"].update(values=origin_values)


def sub_search_method_update_values(window, values, event):
    if values[event] == "Finger Print":
        window["-SUB_SEARCH_THRESHOLD-"].update(disabled=True)
        window["-SUB_SEARCH_MORGAN_VALUES-"].update(disabled=True)
    elif values[event] == "Morgan":
        window["-SUB_SEARCH_THRESHOLD-"].update(disabled=False)
        window["-SUB_SEARCH_MORGAN_VALUES-"].update(disabled=False)
    elif values[event] == "Skeleton":
        window["-SUB_SEARCH_THRESHOLD-"].update(disabled=True)
        window["-SUB_SEARCH_MORGAN_VALUES-"].update(disabled=True)

    else:
        window["-SUB_SEARCH_THRESHOLD-"].update(disabled=False)
        window["-SUB_SEARCH_MORGAN_VALUES-"].update(disabled=True)


def from_assay_updater(window, values, event):

    window["-SUB_SEARCH_ASSAY-"].update(disabled=not values[event])
    window["-SUB_SEARCH_APPROVED_ONLY-"].update(disabled=not values[event])
    window["-SUB_SEARCH_HITS_ONLY-"].update(disabled=not values[event])


def search_daughter_plates_update_values(window):
    window["-SEARCH_PLATE_LAYOUT-"].update(disabled=False)
    window["-SEARCH_MP_MINIMIZED-"].update(disabled=False)
    window["-SEARCH_IGNORE_PLATED_COMPOUNDS-"].update(disabled=True)


def search_mother_plates_update_values(window):
    window["-SEARCH_PLATE_LAYOUT-"].update(disabled=True)
    window["-SEARCH_MP_MINIMIZED-"].update(disabled=True)
    window["-SEARCH_IGNORE_PLATED_COMPOUNDS-"].update(disabled=False)
    window["-SEARCH_PLATE_LAYOUT_SAMPLE_AMOUNT-"].update(value=384)
    window["-SEARCH_PLATE_LAYOUT-"].update(value="")


def search_sample_counter_update(dbf, window, values):
    temp_layout = eval(dbf.find_data_single_lookup("plate_layout", values["-SEARCH_PLATE_LAYOUT-"],
                                                   "layout_name")[0][5])
    temp_counter = 0
    for counter in temp_layout:
        if temp_layout[counter]["state"] == "sample":
            temp_counter += 1
    window["-SEARCH_PLATE_LAYOUT_SAMPLE_AMOUNT-"].update(value=temp_counter)


def list_box_update(window, values):
    new_smiles = values["-SUB_SEARCH_SMILES-"].strip()

    smiles_list = window["-SUB_SEARCH_SMILES_LIST-"].get()
    new_smiles_list = [[new_smiles]]
    for row in smiles_list:
        for smiles in row:
            temp_smiles = [smiles]
        new_smiles_list.append(temp_smiles)
    # smiles_list.append(smiles)
    window["-SUB_SEARCH_SMILES_LIST-"].update(values=new_smiles_list)
    window["-SUB_SEARCH_SMILES-"].update(value="")


def _sub_search_table_data(dbf, compound_data, assay_name, approved_only, hits_only):
    table_data = []
    counter = 0
    for compound in compound_data:
        compound_id = compound
        smiles = compound_data[compound]["smiles"]
        if assay_name:
            assay_scores = dbf.find_data_double_lookup("biological_compound_data", compound_id, assay_name,
                                                       "compound_id", "assay_name")
        else:
            assay_scores = [""]

        try:
            compound_data[compound]["match_score"]
        except KeyError:
            match_score = ""
        else:
            match_score = compound_data[compound]["match_score"]
        for score in assay_scores:
            try:
                score[6]
            except IndexError:
                assay_score = score
            else:
                if approved_only:
                    if score[10] == "1":
                        assay_score = round(float(score[6]), 2)
                    else:
                        continue
                elif hits_only:
                    if score[7] == "1":
                        assay_score = round(float(score[6]), 2)
                    else:
                        continue
                else:
                    assay_score = round(float(score[6]), 2)

            temp_row = [compound_id, match_score, assay_score, smiles]

            table_data.append(temp_row)
            counter += 1
    return table_data, counter


def sub_search(dbf, config, window, values, sub_search_info):
    smiles_list = window["-SUB_SEARCH_SMILES_LIST-"].get()
    if smiles_list:
        smiles_target_list = []
        for smiles in smiles_list:
            smiles_target_list.append(smiles[0])
    else:
        smiles_target_list = [values["-SUB_SEARCH_SMILES-"]]
    print(smiles_target_list)
    sub_search_methode = values["-SUB_SEARCH_METHOD-"]
    if sub_search_methode.casefold() == "skeleton":
        threshold = None
    else:
        threshold = float(values["-SUB_SEARCH_THRESHOLD-"])
    table = config["Tables"]["compound_main"]
    assay_name = values["-SUB_SEARCH_ASSAY-"]
    approved_only = values["-SUB_SEARCH_APPROVED_ONLY-"]
    hits_only = values["-SUB_SEARCH_HITS_ONLY-"]
    fd = FetchData(config)
    rows = fd.sub_structure_search_compound_list(None, table)

    compound_data, hit_list = structure_search(sub_search_methode, threshold, rows, smiles_target_list)

    if not compound_data:
        return sub_search_info

    temp_table_data, sample_counter = _sub_search_table_data(dbf, compound_data, assay_name, approved_only, hits_only)
    sub_search_info.clear()
    sub_search_info = {
        "method": values["-SUB_SEARCH_METHOD-"],
        "smiles_list": smiles_target_list,
        "sample_amount": sample_counter,
        "threshold": threshold,
        "morgan_values": {
            "chirality": values["-SUB_SEARCH_MORGAN_CHIRALITY-"],
            "Features": values["-SUB_SEARCH_MORGAN_FEATURES-"],
            "n bits": values["-SUB_SEARCH_MORGAN_BITS-"],
            "bound range": values["-SUB_SEARCH_MORGAN_RANGE-"],
        },
        "assay_name": assay_name,
        "approved_only": approved_only,
        "hits_only": hits_only,
    }
    window["-SUB_SEARCH_TABLE-"].update(values=temp_table_data)
    window["-SUB_SEARCH_SAMPLE_AMOUNT-"].update(value=sample_counter)
    return sub_search_info


def sub_search_export_table(window, values, sub_search_info):
    export_excel, export_csv = export_chooser_popup()
    print(f"export_excel - {export_excel}")
    print(f"export_csv - {export_csv}")
    if not export_excel and not export_csv:
        print("returning None")
        return None
    else:
        search_table_data = window["-SUB_SEARCH_TABLE-"].get()

        if export_excel:
            sg.PopupError("Not set up yet")
            print("exporting excel")
            # export_excel_search_table(search_table_data, sub_search_info)

        if export_csv:
            sg.PopupError("Not set up yet")
            print("exporting csv")
            # export_csv_search_table(search_table_data, sub_search_info)


if __name__ == "__main__":
    import configparser
    from database_controller import DataBaseFunctions

    config = configparser.ConfigParser()
    config.read("config.ini")
    dbf = DataBaseFunctions(config)
    # sub_search(config, dbf, None, None)














