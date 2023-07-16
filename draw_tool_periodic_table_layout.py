from __future__ import print_function
import sys

from PySide6 import QtCore, QtWidgets
from PySide6.QtGui import QActionGroup, QAction

from draw_tool_periodic_table_info import periodic_table


class PeriodicTable(QtWidgets.QWidget):

    def __init__(self):
        super(PeriodicTable, self).__init__()
        self.periodic_table = periodic_table
        self.init_ui()

    # noinspection PyAttributeOutsideInit
    def init_ui(self):
        grid = QtWidgets.QGridLayout()

        # Create actions dictionary and group dictionary
        self.atom_action_group = QActionGroup(self)
        self.atom_actions = {}

        # Set up the layout for the table
        for key in self.periodic_table.keys():
            atom_name = self.periodic_table[key]["Symbol"]
            temp_atom = QAction(f"{atom_name}", self, statusTip=f"Set atom type to {atom_name}",
                                triggered=self.atom_type_push, objectName=atom_name, checkable=True)

            # Add the atoms and the "action" to the group and dictionary
            self.atom_action_group.addAction(temp_atom)
            self.atom_actions[atom_name] = temp_atom

            if temp_atom.objectName() == "C":
                temp_atom.setChecked(True)

            # Set up the buttons
            button = QtWidgets.QToolButton()
            button.setDefaultAction(temp_atom)
            button.setFocusPolicy(QtCore.Qt.NoFocus)
            button.setMaximumWidth(40)

            if self.periodic_table[key]["Group"] is not None:
                grid.addWidget(button, self.periodic_table[key]["Period"], self.periodic_table[key]["Group"])
            else:
                if key < 72:
                    grid.addWidget(button, 9, key - 54)
                else:
                    grid.addWidget(button, 10, key - 86)

        # Ensure spacing between main table and actinides/lathanides
        grid.addWidget(QtWidgets.QLabel(''), 8, 1)

        self.setLayout(grid)

        self.move(300, 150)
        self.setWindowTitle('Periodic Table')

    atom_type_changed = QtCore.Signal(str, name="atom_type_changed")

    def atom_type_push(self):
        sender = self.sender()
        self.atom_type_changed.emit(sender.objectName())

    # For setting the new atomtype
    def select_atom_type(self, atom_name):
        if atom_name in self.atom_actions.keys():
            self.atom_actions[atom_name].setChecked(True)
        else:
            print(f"Unknown atom type or key error: {atom_name}")


def main():
    """
	Start up the periodic table
	:return:
	"""
    app = QtWidgets.QApplication(sys.argv)
    pt = periodic_table()
    pt.selectAtomtype("N")
    pt.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
