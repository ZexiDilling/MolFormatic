import configparser
import multiprocessing as mp
from multiprocessing import Queue

from gui_controller import main
from draw_tool_handler import launch_draw_tool


def main_controller():

    config = configparser.ConfigParser()
    config.read("config.ini")
    queue_gui = Queue()
    queue_mol = Queue()
    process_gui = mp.Process(target=main, args=(config, queue_gui, queue_mol))
    process_gui.start()

    while True:
        msg, smiles = queue_mol.get()
        print(smiles)
        if msg is None:
            break
        elif msg == "start_draw_tool":
            # Start the PySide6 process
            process_side6 = mp.Process(target=launch_draw_tool, args=(queue_gui, queue_mol, smiles))
            process_side6.start()
            process_side6.join()
        elif msg == "close":
            break


if __name__ == "__main__":
    main_controller()
