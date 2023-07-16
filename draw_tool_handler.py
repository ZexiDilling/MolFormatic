from PySide6.QtWidgets import QApplication

from draw_tool_gui import MainWindow


def launch_draw_tool(queue_gui, queue_mol, molecule=None):
    print("starting to draw")
    app = QApplication([])
    window = MainWindow(queue_gui, queue_mol, molecule=molecule)
    app.exec()


if __name__ == "__main__":
    # print(Path("icons").absolute())
    # mol = Chem.MolFromSmiles('c1ccccc1')
    mol = None
    launch_draw_tool(queue_gui=None, queue_mol=None, molecule=mol)
