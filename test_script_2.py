import PySimpleGUI as sg
import random
import string

BLANK_BOX = '☐'
CHECKED_BOX = '☑'

# ------ Some functions to help generate data for the table ------
def word():
    return ''.join(random.choice(string.ascii_lowercase) for i in range(10))
def number(max_val=1000):
    return random.randint(0, max_val)

def make_table(num_rows, num_cols):
    data = [[j for j in range(num_cols)] for i in range(num_rows)]
    data[0] = [word() for __ in range(num_cols)]
    for i in range(1, num_rows):
        data[i] = [BLANK_BOX if i % 2 else CHECKED_BOX] + [word(), *[number() for i in range(num_cols - 1)]]
    return data

# ------ Make the Table Data ------
data = make_table(num_rows=15, num_cols=6)
headings = [str(data[0][x])+' ..' for x in range(len(data[0]))]
selected = {i for i, row in enumerate(data[1:][:]) if row[0] == CHECKED_BOX}

# ------ Window Layout ------
layout = [[sg.Table(values=data[1:][:], headings=headings, max_col_width=25,
                    auto_size_columns=False,
                    col_widths=[10, 10, 20, 20 ,30, 5],
                    display_row_numbers=True,
                    justification='center',
                    num_rows=20,
                    key='-TABLE-',
                    # selected_row_colors='red on yellow',
                    # enable_events=True,
                    expand_x=False,
                    expand_y=True,
                    vertical_scroll_only=False,
                    enable_click_events=True,
                    tooltip='This is a table', font='_ 14'),
                    sg.Sizegrip()]]
# ------ Create Window ------
window = sg.Window('The Table Element', layout, resizable=True, finalize=True)

window['-TABLE-'].update(values=data[1:][:])

# ------ Event Loop ------
while True:
    event, values = window.read()

    if event == sg.WIN_CLOSED:
        break
    elif event[0] == '-TABLE-' and event[2][1] == 0:

        row = event[2][0]
        print(row)
        print(data)
        if data[row+1][0] == CHECKED_BOX:
            selected.remove(row)
            data[row + 1][0] = BLANK_BOX
        else:
            selected.add(row)
            data[row + 1][0] = CHECKED_BOX

        window['-TABLE-'].update(values=data[1:][:])

window.close()

mapping = {
    "mapping": "Hit Map",
    "bins": {"th_1": {
        "use": True,
        "min": 0.0,
        "max": 10.0,
        "colour": "#a756ad"},
        "th_2": {
            "use": True,
            "min": 10.0,
            "max": 20.0,
            "colour": "#cf347a"},
        "th_3": {
            "use": True,
            "min": 20.0,
            "max": 30.0,
            "colour": "#53d033"},
        "th_4": {
            "use": True,
            "min": 30.0,
            "max": 40.0,
            "colour": "#4993ba"},
        "th_5": {
            "use": True,
            "min": 40.0,
            "max": 50.0,
            "colour": "#57bb48"},
        "th_6": {
            "use": True,
            "min": 50.0,
            "max": 60.0,
            "colour": "#99b053"},
        "th_7": {
            "use": True,
            "min": 60.0,
            "max": 70.0,
            "colour": "#0b8b07"},
        "th_8": {
            "use": True,
            "min": 70.0,
            "max": 80.0,
            "colour": "#b922f2"},
        "th_9": {
            "use": True,
            "min": 80.0,
            "max": 90.0,
            "colour": "#e43086"},
        "th_10": {
            "use": True,
            "min": 90.0,
            "max": 100.0,
            "colour": "#6d83a7"}
    }
}




