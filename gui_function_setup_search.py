def search_compound(window, values, config):
    ac = values["-SEARCH_AC-"]
    origin_values = []
    if ac:
        for acs in ac:
            acs = acs.casefold()
            for values in config[f"database_specific_{acs}"]:
                origin_values.append(config[f"database_specific_{acs}"][values])
    window["-SEARCH_ORIGIN-"].update(values=origin_values)


def sub_search_method_update_values(window, values):
    if values["-SUB_SEARCH_METHOD-"] == "morgan":
        window["-SUB_SEARCH_MORGAN_OPTIONS-"].update(visible=True)
        window["-SUB_SEARCH_MORGAN_CHIRALITY-"].update(visible=True)
        window["-SUB_SEARCH_MORGAN_FEATURES-"].update(visible=True)
        window["-SUB_SEARCH_BITS_TEXT-"].update(visible=True)
        window["-SUB_SEARCH_MORGAN_BITS-"].update(visible=True)
        window["-SUB_SEARCH_BOUND_TEXT-"].update(visible=True)
        window["-SUB_SEARCH_MORGAN_RANGE-"].update(visible=True)
    else:
        window["-SUB_SEARCH_MORGAN_OPTIONS-"].update(visible=False)
        window["-SUB_SEARCH_MORGAN_CHIRALITY-"].update(visible=False)
        window["-SUB_SEARCH_MORGAN_FEATURES-"].update(visible=False)
        window["-SUB_SEARCH_BITS_TEXT-"].update(visible=False)
        window["-SUB_SEARCH_MORGAN_BITS-"].update(visible=False)
        window["-SUB_SEARCH_BOUND_TEXT-"].update(visible=False)
        window["-SUB_SEARCH_MORGAN_RANGE-"].update(visible=False)


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












