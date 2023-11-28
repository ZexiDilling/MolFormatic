from database_startup import DatabaseSetUp
from helpter_functions import config_update
from PySimpleGUI import PopupGetText


def start_up_database(config, db_active, window):
    if db_active:
        dsu = DatabaseSetUp(config, config["Database"]["database"])
    else:
        db = PopupGetText("Choose database name")
        if db is not None:
            db += ".db"
            dsu = DatabaseSetUp(config, db)

    try:
        dsu.controller()
    except UnboundLocalError:
        print("No dp selected")
    else:
        config_update(config)
        window.close()




