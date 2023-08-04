import re
from math import floor
from info import unit_converter_dict

def hex_to_rgb(hex_colour):
    """
    Changes Hex-colour to RBG-colour.
    :param hex_colour: Code for hex-colour
    :type hex_colour: str
    :return: RBG colour as a list
    :rtype: list
    """
    # Pass 16 to the integer function for change of base
    return [int(hex_colour[i:i + 2], 16) for i in range(1, 6, 2)]


def rgb_to_hex(rgb):
    """
    Changes RBG-colour to Hex-colour .
    :param rgb: A list of 3 RGB values
    :type rgb: list
    :return: The colour in hex format, starting with '#'
    :rtype: str
    """
    # Components need to be integers for hex to make sense
    rgb = [int(x) for x in rgb]
    return "#" + "".join(["0{0:x}".format(v) if v < 16 else "{0:x}".format(v) for v in rgb])


def row_col_to_cell(row, col):
    col_names = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J",
                 "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T",
                 "U", "V", "W", "X", "Y", "Z"]
    col -= 1
    if col < len(col_names):
        cell_name = f"{col_names[col]}{row}"
    else:
        stacking_letter = floor(col/len(col_names))
        temp_col = col - (len(col_names) * stacking_letter )
        stacking_letter -= 1
        cell_name = f"{col_names[stacking_letter]}{col_names[temp_col]}{row}"
    return cell_name


def increment_text_string(txt):
    head = txt.rstrip('0123456789')
    tail = txt[len(head):]
    tail = int(tail) + 1
    incremented_text = f"{head}{tail}"
    return incremented_text


def unit_converter(input_value, new_unit_out=None, old_unit_out=False, as_list=False):
    # Validate input_value format
    pattern = r'^(\d+(\.\d+)?(?:e[-+]?\d+)?)([a-zA-Z]+)([a-zA-Z]*)$'
    try:
        re.match(pattern, input_value)
    except TypeError:
        if new_unit_out:
            new_unit = [letter for letter in new_unit_out]
            if len(new_unit) > 1:
                new_unit_out = new_unit[0]
                unit_type = new_unit[1]
            else:
                unit_type = new_unit[0]
        else:
            unit_type = ""
        number_str = input_value

    else:
        match = re.match(pattern, input_value)
        if not match:
            print(input_value)
            raise ValueError("Invalid input format. Input should be in the format '<number><unit>' or '<number><unit><type>'")
        number_str, _, unit_type, _ = match.groups()

    unit_type = [letter for letter in unit_type]

    if len(unit_type) > 1:
        original_unit = unit_type[0]
        temp_type = unit_type[1]
    else:
        original_unit = ""
        try:
            unit_type[0]
        except IndexError:
            temp_type = ""
        else:
            temp_type = unit_type[0]
    # Make the unit case-insensitive
    original_unit = original_unit.lower()

    # Check if the unit is valid
    if original_unit not in unit_converter_dict:
        raise ValueError(f"Invalid unit '{original_unit}'. Available units are: {', '.join(unit_converter_dict.keys())}")

    number = float(number_str)
    base_number = unit_converter_dict[original_unit] * number

    if old_unit_out:
        unit_out = original_unit
    elif new_unit_out:
        unit_out = new_unit_out
    else:
        unit_out = None

    if unit_out:
        new_unit_out = unit_out.lower()  # Make the unit_out case-insensitive

        # Check if the unit_out is valid
        if new_unit_out not in unit_converter_dict:
            raise ValueError(f"Invalid unit_out '{new_unit_out}'. Available units are: {', '.join(unit_converter_dict.keys())}")

        converted_number = base_number / unit_converter_dict[new_unit_out]

    else:
        converted_number = base_number
        new_unit_out = ""

    if as_list:
        return converted_number, new_unit_out, temp_type, original_unit
    else:
        return f"{converted_number}{new_unit_out}{temp_type}"


if __name__ == "__main__":
    print(row_col_to_cell(3, 53))
