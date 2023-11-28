
def heatmap_low_colour_update(window, values):
    if values["-BIO_INFO_HEATMAP_LOW_COLOUR_TARGET-"] != "None":
        window["-BIO_INFO_HEATMAP_LOW_COLOUR_BOX-"]. \
            update(background_color=values["-BIO_INFO_HEATMAP_LOW_COLOUR_TARGET-"])


def heatmap_mid_colour_update(window, values):
    if values["-BIO_INFO_HEATMAP_MID_COLOUR_TARGET-"] != "None":
        window["-BIO_INFO_HEATMAP_MID_COLOUR_BOX-"]. \
            update(background_color=values["-BIO_INFO_HEATMAP_MID_COLOUR_TARGET-"])


def heatmap_high_colour_update(window, values):
    if values["-BIO_INFO_HEATMAP_HIGH_COLOUR_TARGET-"] != "None":
        window["-BIO_INFO_HEATMAP_HIGH_COLOUR_BOX-"]. \
            update(background_color=values["-BIO_INFO_HEATMAP_HIGH_COLOUR_TARGET-"])


def hit_low_colour_update(window, values):
    if values["-BIO_INFO_HIT_MAP_LOW_COLOUR_TARGET-"] != "None":
        window["-BIO_INFO_HIT_MAP_LOW_COLOUR_BOX-"].\
            update(background_color=values["-BIO_INFO_HIT_MAP_LOW_COLOUR_TARGET-"])


def hit_mid_colour_update(window, values):
    if values["-BIO_INFO_HIT_MAP_MID_COLOUR_TARGET-"] != "None":
        window["-BIO_INFO_HIT_MAP_MID_COLOUR_BOX-"]. \
            update(background_color=values["-BIO_INFO_HIT_MAP_MID_COLOUR_TARGET-"])


def hit_high_colour_update(window, values):
    if values["-BIO_INFO_HIT_MAP_HIGH_COLOUR_TARGET-"] != "None":
        window["-BIO_INFO_HIT_MAP_HIGH_COLOUR_BOX-"]. \
            update(background_color=values["-BIO_INFO_HIT_MAP_HIGH_COLOUR_TARGET-"])


def update_bio_info_values(window, values, plate_bio_info):
    if values["-BIO_INFO_STATES-"]:

        temp_plate_name = values["-BIO_INFO_PLATES-"]
        temp_analyse_method = values["-BIO_INFO_ANALYSE_METHOD-"]
        temp_state = values["-BIO_INFO_STATES-"]

        window["-INFO_BIO_AVG-"].update(
            value=plate_bio_info[temp_plate_name]["calculations"][temp_analyse_method][temp_state]["avg"])
        window["-INFO_BIO_STDEV-"].update(
            value=plate_bio_info[temp_plate_name]["calculations"][temp_analyse_method][temp_state]["stdev"])
        window["-INFO_BIO_Z_PRIME-"].update(value=plate_bio_info[temp_plate_name]["calculations"]["other"]["z_prime"])


