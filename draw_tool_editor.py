from math import sqrt
import PySide6
import numpy as np
from PySide6 import QtCore, QtGui
from PySide6.QtGui import QPainter, QPen, Qt, QPixmap, QKeyEvent
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Geometry import Point2D, Point3D

from draw_tool_viewer import MoleculeViewer
from draw_tool_extra_functions import calculate_distance
from draw_tool_periodic_table_info import symbol_to_int


class MoleculeEditor(MoleculeViewer):
    def __init__(self, molecule=None, parent=None,):
        # Init the super class
        super(MoleculeEditor, self).__init__(parent)

        # Activating tracking of the mouse as default -
        # used for seeing if the mouse is close to an atom.
        # If it is the case, mark it with a circle...
        # ToDo
        self.setMouseTracking(True)
        self.setAcceptDrops(True)
        self.setFocusPolicy(Qt.StrongFocus)
        self.closest_object = None
        self.start_coordinates = None
        self.end_coordinates = None
        self.drag_atom = None
        self.dragged_atom = None

        # Creates a stamp dict for adding standard molecules to the canvas
        self.stamp_dict = {
            "Cyclopentane": "C1CCCC1",
            "Benzene": "C1=CC=CC=C1",
            "glucose": "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O",
            "mannose": "C",
            "galactose": "O[C@H]1[C@@H](O)[C@H](O[C@H](O)[C@@H]1O)CO",
            "xylose": "C1[C@H]([C@@H]([C@H](C(O1)O)O)O)O",
            "arabinose": "C([C@H]([C@H]([C@@H](C=O)O)O)O)O",
            "ribose": "C([C@H]([C@H]([C@H](C=O)O)O)O)O",
            "deoxyribose": "C(C=O)[C@@H]([C@@H](CO)O)O",
            "lyxose": "C([C@@H]([C@H]([C@H](C=O)O)O)O)O",
            "allose": "OCC(O)[C@@H](O)[C@@H](O)[C@@H](O)C=O",
            "altrose": "OCC(O)[C@@H](O)[C@@H](O)[C@H](O)C=O",

            "talose": "O=C[C@@H](O)[C@@H](O)[C@@H](O)[C@H](O)CO",
            "idose": "O[C@@H]1[C@@H](O)[C@H](OC(O)[C@H]1O)CO",
        }

        # This is to check if an atom should be replaced or added. Make the value higher to make "click" more likely
        self.click_replace_check = 10

        # Connect the signal to the paint_circle slot
        # self.signal_circle_requested.connect(self.paintEvent)

        # Set the canvas to delete itself when it's closed.
        # noinspection PyUnresolvedReferences
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)

        # Properties
        self._previous_mol = None                             # For undo
        self.coordination_list = None                   # SVG coords of the current mols atoms

        self.view_box = self.renderer().viewBox()
        self.canvas_size = self.size()

        # Set up the atmos and bonds
        self.symbol_to_int = symbol_to_int
        # noinspection PyUnresolvedReferences
        self.bond_types = Chem.rdchem.BondType.names          # A dictionary with all available rdkit bond types

        # Default actions (add, single_bond and Carbon atom
        self._action_tool = "add"                             # Add atoms or bonds to the drawing
        self._bond_type = self.bond_types["SINGLE"]           # single bond - Should always be single bond !!
        self._atom_type = 6                                   # Atom is Carbon.

        # Points to calculate the SVG to coordinates scaling
        self.points = [Point2D(0, 0), Point2D(1, 1)]

        # When drawing finished, update coordination list of SVG atoms.
        self.signal_finished_drawing.connect(self.update_coordlist)

        # Init with a molecule if passed into the constructor on launch
        self.mol = molecule

    # Getters and Setters for properties
    signal_action_tool_changed = QtCore.Signal(name="signal_action_tool_changed")

    @property
    def action_tool(self):
        return self._action_tool

    @action_tool.setter
    def action_tool(self, action_tool_name):
        if action_tool_name != self.action_tool:
            self._action_tool = action_tool_name
            self.signal_action_tool_changed.emit()

    def set_action_tool(self, action_tool_name):
        self.action_tool = action_tool_name

    signal_bond_type_Changed = QtCore.Signal(name="signal_bond_type_Changed")

    @property
    def bond_type(self):
        return self._bond_type

    @bond_type.setter
    def bond_type(self, bond_type):
        if bond_type != self.bond_type:
            self._bond_type = bond_type
            self.signal_bond_type_Changed.emit()

    def set_bond_type(self, bond_type):
        if type(bond_type) == Chem.rdchem.BondType:
            self.bond_type = bond_type
        elif bond_type == "BENZENE":
            print("benzene")
        elif isinstance(bond_type, str):
            assert bond_type in self.bond_types.keys(), f"Bond type {bond_type} not known"
            self.bond_type = self.bond_types[bond_type]
        else:
            print(f"Bond type must be string or rdchem.BondType, not {type(bond_type)}")

    atom_type_changed = QtCore.Signal(name="atom_type_changed")

    @property
    def atom_type(self):
        return self._atom_type

    @atom_type.setter
    def atom_type(self, atom_type):
        if atom_type != self.atom_type:
            self._atom_type = atom_type
            self.atom_type_changed.emit()

    def set_atom_type(self, atom_type):
        if atom_type in self.symbol_to_int.keys():
            self.atom_type = self.symbol_to_int[atom_type]
        elif type(atom_type) == int:

            self.atom_type = atom_type
        else:
            print(f"Atom type must be string or integer, not {atom_type}")

    def update_coordlist(self):
        if self.mol is not None:
            self.coordination_list = np.array([list(self.drawer.GetDrawCoords(atom_index)) for
                                               atom_index in range(self.mol.GetNumAtoms())])
        else:
            self.coordination_list = None

    # Check mouse events
    def mouseMoveEvent(self, event: PySide6.QtGui.QMouseEvent) -> None:
        # Get x- and y-coordinates from the mouse position.
        x = event.pos().x()
        y = event.pos().y()

        self.current_coordinates = (x, y)
        # Gets the canvas to draw on and the size
        view_box = self.renderer().viewBox()
        size = self.size()

        # Rescale, divide by the size of the widget, multiply by the size of the view_box + offset.
        x_svg = float(x)/size.width() * view_box.width() + view_box.left()
        y_svg = float(y)/size.height() * view_box.height() + view_box.top()

        self.closest_object = self.get_mol_object(event)

        if self.drag_atom is not None:
            dragged_atom = self.get_mol_object(event)
            if type(dragged_atom) == Chem.rdchem.Atom and dragged_atom.GetIdx() != self.drag_atom.GetIdx():
                self.dragged_atom = dragged_atom
            if self.action_tool == "Move":
                self.move_atom(self.drag_atom, x_svg, y_svg)
        # self.paintEvent(event)

    def mousePressEvent(self, event):
        """
        Check for mouse button being pressed. both right and left.
        :param event:
        :return:
        """
        clicked = self.get_mol_object(event)
        # noinspection PyUnresolvedReferences
        if event.button() == QtCore.Qt.LeftButton:
            self.start_coordinates = (event.pos().x(), event.pos().y())

            # noinspection PyUnresolvedReferences
            if type(clicked) == Chem.rdchem.Atom:
                self.drag_atom = clicked
            elif type(clicked) == Chem.rdchem.Bond:
                self.bond_click(clicked)
            elif type(clicked) == Point2D:
                self.canvas_click(clicked)

        elif event.button() == QtCore.Qt.RightButton:
            if type(clicked) == Chem.rdchem.Atom:
                self.remove_atom(clicked)
            elif type(clicked) == Chem.rdchem.Bond:
                self.remove_bond(clicked)

    def mouseReleaseEvent(self, event):
        self.end_coordinates = (event.pos().x(), event.pos().y())

        if self.drag_atom is not None and self.dragged_atom is not None:
            # Update the bond between the drag_atom and dragged_atom
            self.update_bond(self.drag_atom, self.dragged_atom)
        elif self.drag_atom is not None and self.dragged_atom is None:
            # This is a check for making the interface nicer to work with. for how to add or replacing atoms...
            clicker_check = calculate_distance(self.start_coordinates, self.end_coordinates)
            if self.action_tool == "add" and clicker_check > self.click_replace_check:
                self.add_atom_to(self.drag_atom)
            else:
                self.atom_click(self.drag_atom)

        self.drag_atom = None
        self.dragged_atom = None

    def keyPressEvent(self, event: PySide6.QtGui.QKeyEvent) -> None:
        if event.key() == Qt.Key_Delete:
            if type(self.closest_object) == Chem.rdchem.Atom:
                self.remove_atom(self.closest_object)
            elif type(self.closest_object) == Chem.rdchem.Bond:
                self.remove_bond(self.closest_object)

        elif event.key() == Qt.Key_H and type(self.closest_object) == Chem.rdchem.Atom:
            self.replace_atom(self.closest_object, 1)

        elif event.key() == Qt.Key_C and type(self.closest_object) == Chem.rdchem.Atom:
            self.replace_atom(self.closest_object, 6)

        elif event.key() == Qt.Key_N and type(self.closest_object) == Chem.rdchem.Atom:
            self.replace_atom(self.closest_object, 7)

        elif event.key() == Qt.Key_O and type(self.closest_object) == Chem.rdchem.Atom:
            self.replace_atom(self.closest_object, 8)

        elif event.key() == Qt.Key_F and type(self.closest_object) == Chem.rdchem.Atom:
            self.replace_atom(self.closest_object, 9)

        super().keyPressEvent(event)

    def canvas_click(self, point):

        if self.action_tool == "add":
            self.add_canvas_atom(point)
        elif self.action_tool in self.stamp_dict:
            smiles = self.stamp_dict[self.action_tool]
            self.add_canvas_stamp(smiles)
        else:
            print(f"The combination of Canvas click and Action %{self.action_tool} undefined")

    def get_mol_object(self, event):
        """
        Gets the molecule object
        :param event:
        :return:
        """

        # Get x- and y-coordinates from the mouse position.
        x = event.pos().x()
        y = event.pos().y()

        # Gets the canvas to draw on and the size
        view_box = self.renderer().viewBox()
        size = self.size()

        # Rescale, divide by the size of the widget, multiply by the size of the view_box + offset.
        x_svg = float(x)/size.width() * view_box.width() + view_box.left()
        y_svg = float(y)/size.height() * view_box.height() + view_box.top()

        # Check distant to the nearest object (atom or bond)
        atom_index, atom_distance = self.get_nearest_atom(x_svg, y_svg)
        bond_index, bond_distance = self.get_nearest_bond(x_svg, y_svg)

        # Check if the click is near enough to an atom or bond, to bind to it or change it
        if min([atom_distance, bond_distance]) < self.distance_check:
            # Check what is closer, the atom or the bond
            if atom_distance < bond_distance:
                return self.mol.GetAtomWithIdx(int(atom_index))
            else:
                return self.mol.GetBondWithIdx(int(bond_index))

        # If nothing is close enough, the software will draw a new thingy
        else:
            return self.svg_to_2d_coordinates(x_svg, y_svg)

    def atom_click(self, atom):
        if self.action_tool == "add":
            if atom.GetAtomicNum() != self.atom_type:
                self.replace_atom(atom)
            else:
                self.add_atom_to(atom)
        elif self.action_tool == "benzene":
            self.add_benzene_to(atom)
        elif self.action_tool == "cyclopentane":
            self.add_cyclopentane_to(atom)
        elif self.action_tool == "eraser":
            self.remove_atom(atom)
        # elif self.action == "Select":
        #     self.select_atom_add(atom)
        # elif self.action == "Replace":
        #     self.replace_atom(atom)
        elif self.action_tool == "increase_charge":
            self.increase_charge(atom)
        elif self.action_tool == "decrease_charge":
            self.decrease_charge(atom)
        # elif self.action == "RStoggle":
        #     self.toogleRS(atom)

        elif self.action_tool == "down":
            # ToDO figure out the right behavior of these...
            self.chiral_type_down_bond(atom)
        elif self.action_tool == "up":
            # ToDO figure out the right behavior of these...
            self.chiral_type_up_bond(atom)
        elif self.action_tool == "empty":
            pass

        else:
            print(f"ATOM_CLICK - The combination of Atom click and Action {self.action_tool} undefined")

    def bond_click(self, bond):
        # if a bond is clicked, and it is on add-mode change between bond types (single, double, triple)
        if self.action_tool == "add":
            self.toggle_bond(bond)
        # elif self.action_tool == "Single Bond":
        #     self.set_to_single_bond(bond)
        elif self.action_tool == "eraser":
            self.remove_bond(bond)

        # elif self.action_tool == "EZtoggle":
        #     self.toogleEZ(bond)
        else:
            print(f"The combination of Bond click and Action {self.action_tool} undefined")

    def get_nearest_atom(self, x_svg, y_svg):
        if self.mol is not None and self.mol.GetNumAtoms() > 0:
            atom_svg_coordinates = np.array([x_svg, y_svg])
            deltas = self.coordination_list - atom_svg_coordinates
            distance_2 = np.einsum('ij,ij->i', deltas, deltas)
            min_index = np.argmin(distance_2)
            return min_index, sqrt(distance_2[min_index])
        else:
            return None, 1e10       # Return a hugh number. if set to None, other code needs to be change

    def get_nearest_bond(self, x_svg, y_svg):
        if self.mol is not None and self.mol.GetNumAtoms() > 2:
            bond_list = []
            for bond in self.mol.GetBonds():
                start_index = bond.GetBeginAtomIdx()
                end_index = bond.GetEndAtomIdx()
                svg_coordinates = np.mean(self.coordination_list[[start_index, end_index]], axis=0)
                bond_list.append(svg_coordinates)

            bond_list = np.array(bond_list)

            atom_svg_coordinates = np.array([x_svg, y_svg])
            deltas = bond_list - atom_svg_coordinates
            distance_2 = np.einsum('ij,ij->i', deltas, deltas)
            min_index = np.argmin(distance_2)
            return min_index, sqrt(distance_2[min_index])
        else:
            return None, 1e10       # Return a hugh number. if set to None, other code needs to be change #ToDo time it

    def svg_to_2d_coordinates(self, x_svg, y_svg):
        """
        Function to translate from SVG coordinates to atom coordinates to generate svg-file,
        using scaling calculated from atom coordinates (0,0) and (1,1)

        :param x_svg:
        :type x_svg: int
        :param y_svg:
        :type y_svg: int
        :return: rdkit Point2D
        :rtype: rdkit.Geometry.rdGeometry.Point2D
        """
        if self.drawer is not None:
            scale0 = self.drawer.GetDrawCoords(self.points[0])
            scale1 = self.drawer.GetDrawCoords(self.points[1])

            ax = scale1.x - scale0.x
            bx = scale0.x

            ay = scale1.y - scale0.y
            by = scale0.y

            return Point2D((x_svg - bx) / ax, (y_svg - by) / ay)
        else:
            return Point2D(0., 0.)

    def move_atom(self, atom, x_svg, y_svg):
        """
        Move the specified atom to the given SVG coordinates.

        :param atom: The atom to be moved.
        :type atom: Chem.rdchem.Atom
        :param x_svg: The x-coordinate in SVG coordinates.
        :type x_svg: float
        :param y_svg: The y-coordinate in SVG coordinates.
        :type y_svg: float
        """
        pass
        # ToDO Make this work
        # rwmol = Chem.rdchem.RWMol(self.mol)
        #
        # atom_idx = atom.GetIdx()
        # conformer = self.mol.GetConformer()
        # conformer.SetAtomPosition(atom_idx, (x_svg, y_svg, 0.0))
        # self.signal_molecule_changed.emit()
        # # self.mol = rwmol

    def add_canvas_atom(self, point):
        """
        Add an atom to the canvas at the given SVG coordinates.

        :param point: The SVG coordinates where the atom will be added.
        :type point: rdkit.Geometry.rdGeometry.Point2D
        """
        rwmol = Chem.rdchem.RWMol(self.mol)
        if rwmol.GetNumAtoms() == 0:
            point.x = 0.0
            point.y = 0.0
        newatom = Chem.rdchem.Atom(self.atom_type)
        newidx = rwmol.AddAtom(newatom)

        # This only triggers if the canvas is empty
        if not rwmol.GetNumConformers():
            rdDepictor.Compute2DCoords(rwmol)
        conf = rwmol.GetConformer(0)
        p3 = Point3D(point.x, point.y, 0)
        conf.SetAtomPosition(newidx, p3)
        self.mol = rwmol
        # print(f"conf: {conf}, p3: {p3}: conf: {conf}")

    def add_canvas_stamp(self, smiles):
        """
        Add a molecule stamp to the canvas based on the provided SMILES.

        :param smiles: The SMILES representation of the molecule to be added.
        :type smiles: str
        """
        benzene = Chem.MolFromSmiles(smiles)  # Create a molecule from smiles

        rwmol = Chem.rdchem.RWMol(self.mol)
        if rwmol.GetNumAtoms() == 0:
            self.mol = benzene
        else:
            rwmol = Chem.CombineMols(rwmol, benzene)
            self.mol = rwmol

    def add_benzene_to(self, atom):
        """
        Add a benzene ring to the molecule, connected to the specified atom.

        :param atom: The atom to which the benzene ring will be connected.
        :type atom: Chem.rdchem.Atom
        """
        rwmol = Chem.rdchem.RWMol(self.mol)
        benzene = Chem.MolFromSmiles('c1ccccc1')

        # Add the benzene ring atoms to the current molecule
        benzene_atom_indices = []
        for benzene_atom in benzene.GetAtoms():
            new_atom_index = rwmol.AddAtom(benzene_atom)
            benzene_atom_indices.append(new_atom_index)

        # Create bonds between the atoms of the benzene ring
        bond_order = self.bond_type  # Set the desired bond order

        # Add single bonds between adjacent ring atoms
        for carbon in range(len(benzene_atom_indices) - 1):
            rwmol.AddBond(benzene_atom_indices[carbon], benzene_atom_indices[carbon + 1], order=Chem.BondType.SINGLE)
        rwmol.AddBond(benzene_atom_indices[-1], benzene_atom_indices[0], order=Chem.BondType.SINGLE)

        # Create a bond between the clicked atom and one of the ring atoms
        rwmol.AddBond(atom.GetIdx(), benzene_atom_indices[0], order=bond_order)

        # Set double bond between adjacent ring atoms
        for carbons in range(len(benzene_atom_indices) - 1):
            if carbons % 2:
                atom_idx1 = benzene_atom_indices[carbons]
                atom_idx2 = benzene_atom_indices[carbons + 1]
                bond = rwmol.GetBondBetweenAtoms(atom_idx1, atom_idx2)
                if bond is not None:
                    bond.SetBondType(Chem.BondType.DOUBLE)

        # Set double bond between the first and last ring atoms
        atom_idx1 = benzene_atom_indices[-1]
        atom_idx2 = benzene_atom_indices[0]
        bond = rwmol.GetBondBetweenAtoms(atom_idx1, atom_idx2)
        if bond is not None:
            bond.SetBondType(Chem.BondType.DOUBLE)

        self.mol = rwmol

    def add_cyclopentane_to(self, atom):
        """
        Add a benzene ring to the molecule, connected to the specified atom.

        :param atom: The atom to which the benzene ring will be connected.
        :type atom: Chem.rdchem.Atom
        """
        rwmol = Chem.rdchem.RWMol(self.mol)
        cyclopentane = Chem.MolFromSmiles('C1CCCC1')

        # Add the benzene ring atoms to the current molecule
        cyclopentane_atom_indices = []
        for temp_atom in cyclopentane.GetAtoms():
            new_atom_index = rwmol.AddAtom(temp_atom)
            cyclopentane_atom_indices.append(new_atom_index)

        # Create bonds between the atoms of the benzene ring
        bond_order = self.bond_type  # Set the desired bond order

        # Add single bonds between adjacent ring atoms
        for carbon in range(len(cyclopentane_atom_indices) - 1):
            rwmol.AddBond(cyclopentane_atom_indices[carbon],
                          cyclopentane_atom_indices[carbon + 1], order=Chem.BondType.SINGLE)
        rwmol.AddBond(cyclopentane_atom_indices[-1], cyclopentane_atom_indices[0], order=Chem.BondType.SINGLE)

        # Create a bond between the clicked atom and one of the ring atoms
        rwmol.AddBond(atom.GetIdx(), cyclopentane_atom_indices[0], order=bond_order)

        self.mol = rwmol

    def add_atom_to(self, atom):
        """
        Add an atom to the molecule, connected to the specified atom.

        :param atom: The atom to which the new atom will be connected.
        :type atom: Chem.rdchem.Atom
        """
        rwmol = Chem.rdchem.RWMol(self.mol)
        newatom = Chem.rdchem.Atom(self.atom_type)

        newidx = rwmol.AddAtom(newatom)
        rwmol.AddBond(atom.GetIdx(), newidx, order=self.bond_type)

        self.mol = rwmol

    def add_bond(self, atom):
        """
        Add a bond between the previously selected atom and the specified atom.

        :param atom: The atom to which the bond will be added.
        :type atom: Chem.rdchem.Atom
        """
        if len(self.selected_atoms) > 0:
            selected = self.selected_atoms[-1]
            rwmol = Chem.rdchem.RWMol(self.mol)
            rwmol.AddBond(selected, atom.GetIdx(), order=self.bond_type)
            self.mol = rwmol
        else:
            print("add_bond - failed")

    def remove_atom(self, atom):
        """
        Remove the specified atom from the molecule.

        :param atom: The atom to be removed.
        :type atom: Chem.rdchem.Atom
        """
        rwmol = Chem.rdchem.RWMol(self.mol)
        rwmol.RemoveAtom(atom.GetIdx())
        self.clear_atom_selection()             # Removing atoms updates Idx'es
        self.mol = rwmol

    def update_bond(self, atom_start, atom_end):
        """
        Update the bond between the specified start and end atoms.

        :param atom_start: The start atom of the bond.
        :type atom_start: Chem.rdchem.Atom
        :param atom_end: The end atom of the bond.
        :type atom_end: Chem.rdchem.Atom
        """
        atom_start_idx = atom_start.GetIdx()
        atom_end_idx = atom_end.GetIdx()

        rwmol = Chem.RWMol(self.mol)
        rwmol.AddBond(atom_start_idx, atom_end_idx, order=self.bond_type)
        self.mol = rwmol

    def increase_charge(self, atom):
        """
        Increase the formal charge of the specified atom by 1.

        :param atom: The atom whose charge will be increased.
        :type atom: Chem.rdchem.Atom
        """
        self.backupMol()
        atom.SetFormalCharge(atom.GetFormalCharge()+1)
        self.signal_molecule_changed.emit()

    def decrease_charge(self, atom):
        """
        Decrease the formal charge of the specified atom by 1.

        :param atom: The atom whose charge will be decreased.
        :type atom: Chem.rdchem.Atom
        """
        self.backupMol()
        atom.SetFormalCharge(atom.GetFormalCharge()-1)
        self.signal_molecule_changed.emit()

    def replace_atom(self, atom, new_atom=None):
        """
        Replace the specified atom with a new atom of the selected type.

        :param atom: The atom to be replaced.
        :type atom: Chem.rdchem.Atom
        """
        rwmol = Chem.rdchem.RWMol(self.mol)
        if new_atom is None:
            newatom = Chem.rdchem.Atom(self.atom_type)
        else:
            newatom = Chem.rdchem.Atom(new_atom)
        rwmol.ReplaceAtom(atom.GetIdx(), newatom)
        self.mol = rwmol

    # noinspection PyUnresolvedReferences
    def toggle_bond(self, bond):
        """
        Toggle the bond type of the specified bond.

        :param bond: The bond whose type will be toggled.
        :type bond: Chem.rdchem.Bond
        """
        self.backupMol()
        bond_type = bond.GetBondType()
        bond_types = [Chem.rdchem.BondType.TRIPLE,
                      Chem.rdchem.BondType.SINGLE,
                      Chem.rdchem.BondType.DOUBLE,
                      Chem.rdchem.BondType.TRIPLE]

        # Find the next type in the list based on current
        # If current is not in list? Then it selects the first and add 1 => SINGLE
        new_index = np.argmax(np.array(bond_types) == bond_type)+1
        new_type = bond_types[new_index]
        bond.SetBondType(new_type)
        self.signal_molecule_changed.emit()

    def remove_bond(self, bond):
        """
        Remove the specified bond from the molecule.

        :param bond: The bond to be removed.
        :type bond: Chem.rdchem.Bond
        """
        rwmol = Chem.rdchem.RWMol(self.mol)
        rwmol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        self.mol = rwmol

    def undo(self):
        """
        Undo the last modification to the molecule.
        """
        self.mol = self._previous_mol

    def chiral_type_down_bond(self, atom):
        """
        Set the chiral type of the specified atom to "down".

        :param atom: The atom to be modified.
        :type atom: Chem.rdchem.Atom
        """
        self.backupMol()

        stereotype = Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW
        atom.SetChiralTag(stereotype)
        self._mol.UpdatePropertyCache()
        rdDepictor.Compute2DCoords(self._mol)
        self.signal_molecule_changed.emit()

    def chiral_type_up_bond(self, atom):
        """
        Set the chiral type of the specified atom to "up".

        :param atom: The atom to be modified.
        :type atom: Chem.rdchem.Atom
        """
        self.backupMol()

        stereotype = Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW
        atom.SetChiralTag(stereotype)
        self._mol.UpdatePropertyCache()
        rdDepictor.Compute2DCoords(self._mol)
        self.signal_molecule_changed.emit()

    def backupMol(self):
        """
        Create a backup of the current molecule for undo functionality.
        """
        self._previous_mol = Chem.Mol(self.mol.ToBinary())


if __name__ == "__main__":
    pass
