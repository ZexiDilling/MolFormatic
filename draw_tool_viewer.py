import numpy as np
from PySide6 import QtSvgWidgets, QtCore
from PySide6.QtCore import Qt, QRect
from PySide6.QtGui import QPainter, QPen
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point2D


class MoleculeViewer(QtSvgWidgets.QSvgWidget):
    def __init__(self, mol=None, parent=None):

        # Also init the super class
        super(MoleculeViewer, self).__init__(parent)

        # Private Properties
        self._mol = None                # The molecule
        self._draw_mol = None           # Molecule for drawing
        self.drawer = None              # Drawing object for producing SVG
        self._selected_atoms = []       # List of selected atoms
        self._closest_object = None     # Getting the closest object for drawing location circles
        self.coordination_list = None                         # SVG coordinates of the current molecules atoms
        self.distance_check = 50                              # Sets the distance for tracking the nearest atom or bond

        # Bind signales to slots for automatic actions
        self.signal_molecule_changed.connect(self.sanitize_draw)
        self.signal_selection_style_changed.connect(self.draw)


    # Getter and setter for mol
    signal_molecule_changed = QtCore.Signal(name="signal_molecule_changed")

    @property
    def mol(self):
        return self._mol

    @mol.setter
    def mol(self, mol):
        if mol is None:
            mol = Chem.MolFromSmiles('')
        if mol != self._mol:

            if self._mol is not None:
                self._previous_mol = Chem.Mol(self._mol.ToBinary())  # Copy
            self._mol = mol
            self.signal_molecule_changed.emit()

    def setMol(self, mol):
        self.mol = mol

    # Handling of selections
    signal_selection_style_changed = QtCore.Signal(name="signal_selection_style_changed")

    def clear_atom_selection(self):
        if self._selected_atoms:
            self._selected_atoms = []
            self.signal_selection_style_changed.emit()

    @property
    def selected_atoms(self):
        return self._selected_atoms

    @selected_atoms.setter
    def selected_atoms(self, atom_list):
        if atom_list != self._selected_atoms:
            assert type(atom_list) == list, "selected_atoms should be a list of integers"
            assert all(isinstance(item, int) for item in atom_list), "selected_atoms should be a list of integers"
            self._selected_atoms = atom_list
            self.signal_selection_style_changed.emit()

    def set_selected_atoms(self, atom_list):
        self.selected_atoms = atom_list

    signal_circle_requested = QtCore.Signal(name="signal_circle_requested")

    @property
    def closest_object(self):
        return self._closest_object

    @closest_object.setter
    def closest_object(self, closest_object):
        if closest_object != self._closest_object:
            self._closest_object = closest_object
            self.signal_circle_requested.emit()

    def set_closest_object(self, closest_object):
        self.closest_object = closest_object

    # Actions and functions
    @QtCore.Slot()
    def draw(self):
        svg = self.get_mol_svg()
        self.load(QtCore.QByteArray(svg.encode('utf-8')))

    @QtCore.Slot()
    def sanitize_draw(self):
        self.sanitizeMol()
        self.draw()

    def compute_new_coordinates(self, ignoreExisting=False, canonOrient=False):
        """Computes new coordinates for the molecule taking into account all
        existing positions (feeding these to the rdkit coordinate generation as
        prev_coords).
        """
        # This code is buggy when you are not using the CoordGen coordinate
        # generation system, so we enable it here
        rdDepictor.SetPreferCoordGen(True)
        prev_coords = {}
        if self._mol.GetNumConformers() == 0:
            pass
            # print("No Conformers found, computing all 2D coords")
        elif ignoreExisting:
            pass
            # print("Ignoring existing conformers, computing all 2D coords")
        else:
            assert self._mol.GetNumConformers() == 1
            # print("1 Conformer found, computing 2D coords not in found conformer")
            conf = self._mol.GetConformer(0)
            for a in self._mol.GetAtoms():
                pos3d = conf.GetAtomPosition(a.GetIdx())
                if (pos3d.x, pos3d.y) == (0, 0):
                    continue
                prev_coords[a.GetIdx()] = Point2D(pos3d.x, pos3d.y)

        rdDepictor.Compute2DCoords(self._mol, coordMap=prev_coords, canonOrient=canonOrient)

    signal_sanitize = QtCore.Signal(str, name="signal_sanitize")

    def clean_up_structure(self):
        self.compute_new_coordinates(canonOrient=True, ignoreExisting=True)
        self._draw_mol = Chem.Mol(self._mol.ToBinary())
        self.draw()

    @QtCore.Slot()
    def sanitizeMol(self, kekulize=False, draw_kekulize=False):

        self.compute_new_coordinates()
        self._draw_mol = Chem.Mol(self._mol.ToBinary()) #Is this necessary?
        try:
            Chem.SanitizeMol(self._draw_mol)
            self.signal_sanitize.emit("Sanitizable")
        except:
            print("hit exception - UNSANITIZABLE")
            self.signal_sanitize.emit("UNSANITIZABLE")
            try:
                self._draw_mol.UpdatePropertyCache(strict=False)
            except:
                self.signal_sanitize.emit("UpdatePropertyCache FAIL")

        # Kekulize
        if kekulize:
            try:
                Chem.Kekulize(self._draw_mol)
            except:
                print("hit exception - Unkekulizable")
        try:
            self._draw_mol = rdMolDraw2D.PrepareMolForDrawing(self._draw_mol, kekulize=draw_kekulize)
        except ValueError:  # <- can happen on a kekulization failure
            self._draw_mol = rdMolDraw2D.PrepareMolForDrawing(self._draw_mol, kekulize=False)

    signal_finished_drawing = QtCore.Signal(name="signal_finished_drawing")

    def get_mol_svg(self):
        self.drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)

        if self._draw_mol is not None:
            # Chiral tags on R/S
            chiral_tags = Chem.FindMolChiralCenters(self._draw_mol)
            opts = self.drawer.drawOptions()

            for tag in chiral_tags:
                tag_index = tag[0]
                opts.atomLabels[tag_index] = self._draw_mol.GetAtomWithIdx(tag_index).GetSymbol() + ':' + tag[1]

            if len(self._selected_atoms) > 0:
                colors = {self._selected_atoms[-1]: (1, 0.2, 0.2)}          # Color lastly selected a different color
                self.drawer.DrawMolecule(self._draw_mol, highlightAtoms=self._selected_atoms,
                                         highlightAtomColors=colors, )
            else:
                self.drawer.DrawMolecule(self._draw_mol)

        self.drawer.FinishDrawing()
        self.signal_finished_drawing.emit()                                 # Signal that drawer has finished
        svg = self.drawer.GetDrawingText().replace('svg:', '')
        return svg

    def paint_circle(self, painter=QPainter):
        # pass
        # Call the base class implementation first to ensure any default drawing is performed
        # ToDo Make this work
        if self.closest_object is not None:

            # Get the coordinates of the closest atom or bond and re-shape the array to make it easier to get the middle
            data = (self.coordination_list[self.closest_object]).reshape(-1, 2)

            # Calculate the center of the closest object
            center = np.mean(data, axis=0)

            # Convert center coordinates to integers
            center_int = center.astype(int)
            print(f"I am at {center_int}")
            # Sets the radius to the trigger distance of adding a new bond
            radius = self.distance_check

            # Calculate the bounding rectangle coordinates
            x = int(center_int[0] - radius)
            y = int(center_int[1] - radius)
            width = height = radius * 2

            # Set the pen properties for drawing the circle
            pen = QPen(Qt.green)
            pen.setWidth(10)
            painter.setPen(pen)

            # Draw the green circle around the closest object
            rect = QRect(x, y, width, height)
            painter.drawEllipse(rect)

    # def paintEvent(self, event):
        # ToDo This is triggered every time I click on something... figure out what triggers it
    #     if isinstance(event, QtGui.QPaintEvent):
    #         super().paintEvent(event)
    #     painter = QPainter(self)
    #     self.paint_circle(painter)
    #
    #