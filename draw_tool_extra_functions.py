from math import sqrt
from PySide6.QtWidgets import QDialog, QVBoxLayout, QLabel, QPushButton, QLineEdit
from rdkit.Chem import MolToSmiles
from rdkit.Chem.Descriptors import ExactMolWt


class CustomDialog(QDialog):
    def __init__(self, mol, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout()
        self.setWindowTitle('Molecule Info')
        self.setLayout(layout)
        self.mol = mol

        self.label = QLabel("Smiles")
        layout.addWidget(self.label)

        self.smiles_field = QLineEdit(str(MolToSmiles(self.mol)))  # Added QLineEdit widget for mass field
        layout.addWidget(self.smiles_field)

        self.label = QLabel("Exact Mass")
        layout.addWidget(self.label)

        self.mass_field = QLineEdit(str(ExactMolWt(self.mol)))  # Added QLineEdit widget for mass field
        layout.addWidget(self.mass_field)

        update_button = QPushButton("Update")
        update_button.clicked.connect(self.update_data)
        layout.addWidget(update_button)

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.close)
        layout.addWidget(close_button)

    def update_data(self):
        # Perform your mass update logic here
        # For demonstration purposes, let's just double the mass value

        # self.smiles_field.setText(str(MolToSmiles(mol)))
        # self.mass_field.setText(str(ExactMolWt(mol)))
        self.smiles_field.setText("Dose not work!")
        self.mass_field.setText("Open a new window")


def calculate_distance(coord1, coord2):
    x1, y1 = coord1
    x2, y2 = coord2
    distance = sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
    return distance