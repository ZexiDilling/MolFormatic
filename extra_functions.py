
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

