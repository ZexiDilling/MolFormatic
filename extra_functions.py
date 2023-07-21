from math import floor


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


if __name__ == "__main__":
    print(row_col_to_cell(3, 53))
