from PySide6 import QtWidgets, QtCore
from PySide6.QtGui import QIcon, QAction, QActionGroup
from PySide6.QtWidgets import QStatusBar, QMessageBox, QLineEdit, QMenu, QToolButton, QApplication
from pathlib import Path
from rdkit import Chem

from draw_tool_editor import MoleculeEditor
from draw_tool_extra_functions import CustomDialog
from draw_tool_periodic_table_layout import PeriodicTable


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, queue_gui=None, queue_mol=None, molecule=None):
        super(MainWindow, self).__init__()

        # Gets the path to all the icons for the menus
        self.icon_path = Path("icons").absolute()

        self.queue_gui = queue_gui
        self.queue_mol = queue_mol

        # Initialing the editor:
        self.editor = MoleculeEditor()

        # Initialize the toolbars
        self.main_toolbar = self.addToolBar("Main")
        self.left_side_toolbar = QtWidgets.QToolBar(self)

        # Creates a group of actions, to makes sure that you can only choose one action at the time
        self.actionActionGroup = QActionGroup(self)

        # This is the Periodic table for getting compounds
        self.periodic_table = PeriodicTable()

        self.init_gui()
        print(f"Mol: {molecule}")
        self.periodic_table.atom_type_changed.connect(self.set_atom_type_name)

        if molecule is not None:
            mol = Chem.MolFromSmiles(molecule)
            self.load_mol(mol)

    def load_mol(self, mol):
        self.editor.mol = mol

    def init_gui(self):
        # Set up the main window - Title - Icon - Size
        self.setWindowTitle("A simple mol editor")
        self.setWindowIcon(QIcon(str(self.icon_path/"appicon.svg.png")))
        self.setGeometry(100, 100, 200, 150)

        # Set up toolbar
        self.set_up_components()

        self.show()

    def set_up_components(self):
        """
        Starts up the components to the main window, like the central editor canvas, the menu, the toolbars...
        :return:
        """
        self.central_window()                       # editor canvas
        self.status_bar()                           # The statusbar
        self.create_menu_actions()                  # Sets the action for the menu buttons
        self.menu_layout()                          # Set up the menu
        self.create_top_toolbar_actions()           # Creates the different actions that the toolbar can use
        self.top_toolbar_layout()                   # Top toolbar controls editor tools, like:
        self.create_left_side_toolbar_actions()     # Creates the different actions that the toolbar can use
        self.left_side_toolbar_layout()             # Left toolbar controls the drawing tool,

    def central_window(self):
        """
        Set up the center of the window
        :return:
        """
        self.center = self.editor
        self.center.setFixedSize(600, 600)
        self.setCentralWidget(self.center)

    def status_bar(self):
        """
        The button bar in the window that shows the status of stuff
        :return:
        """
        # noinspection PyAttributeOutsideInit
        self.my_status_bar = QStatusBar()
        self.setStatusBar(self.my_status_bar)
        self.my_status_bar.showMessage("ARG!!!!!!!!!", 10000)

    # noinspection PyAttributeOutsideInit
    def menu_layout(self):
        """
        The menu at the top of the window, with file, info and other settings
        :return:
        """
        self.file_menu = self.menuBar().addMenu("File")
        self.info_menu = self.menuBar().addMenu("Info")

        self.file_menu.addAction(self.action_menu_exit)

        self.info_menu.addAction(self.action_menu_about)
        self.info_menu.addAction(self.action_menu_about_qt)

    def top_toolbar_layout(self):
        """
        The Main toolbar at the top of the window
        :return:
        """
        self.main_toolbar.addAction(self.action_top_undo)
        # self.main_toolbar.addAction(self.action_top_move)
        self.main_toolbar.addSeparator()
        self.main_toolbar.addAction(self.action_top_clear_canvas)
        self.main_toolbar.addAction(self.action_top_clean_structure)
        self.main_toolbar.addSeparator()
        self.main_toolbar.addAction(self.action_top_show_mol_info)
        self.main_toolbar.addAction(self.action_top_show_smiles)
        self.main_toolbar.addAction(self.action_top_export_mol)
        self.main_toolbar.addAction(self.action_top_import_mol)
        self.main_toolbar.addWidget(self.smiles_field)

    def left_side_toolbar_layout(self):

        self.addToolBar(QtCore.Qt.LeftToolBarArea, self.left_side_toolbar)
        self.main_toolbar.addSeparator()
        # Set up the buttons
        # self.left_side_toolbar.addAction(self.action_left_side_lasso)
        # self.left_side_toolbar.addAction(self.action_left_side_arrow_circle)
        self.left_side_toolbar.addAction(self.action_left_side_single_standard_bond)
        # self.left_side_toolbar.addAction(self.action_left_side_line_daz)
        # self.left_side_toolbar.addAction(self.action_left_side_squar_daze)
        self.left_side_toolbar.addAction(self.action_left_side_down_bond)
        self.left_side_toolbar.addAction(self.action_left_side_up_bond)
        # self.left_side_toolbar.addAction(self.action_left_side_empty_triangle)
        # self.left_side_toolbar.addAction(self.action_left_side_squvikeliqy)
        # self.left_side_toolbar.addAction(self.action_left_side_multiple_bond)
        # self.left_side_toolbar.addAction(self.action_left_side_selection_tool)
        # self.left_side_toolbar.addAction(self.action_left_side_bond_cutter)
        self.left_side_toolbar.addAction(self.action_left_side_eraser)
        # self.left_side_toolbar.addAction(self.action_left_side_Letters)
        self.left_side_toolbar.addAction(self.action_left_side_increase_charge)
        self.left_side_toolbar.addAction(self.action_left_side_decrease_charge)
        self.main_toolbar.addSeparator()
        for action in self.atom_actions:
            self.left_side_toolbar.addAction(action)
        self.main_toolbar.addSeparator()
        self.left_side_toolbar.addAction(self.action_left_side_benzene)
        self.left_side_toolbar.addAction(self.action_left_side_cyclopentane)
        self.main_toolbar.addSeparator()
        self.left_side_toolbar.addAction(self.action_left_side_open_periodic_table)

        # Add the tool button to the toolbar
        stamps_button = self.create_left_side_toolbar_action_stamps_button()
        self.left_side_toolbar.addWidget(stamps_button)

    # noinspection PyAttributeOutsideInit
    def create_menu_actions(self):
        self.action_menu_exit = QAction("Exit",
                                        self, shortcut="Ctrl+Q",
                                        statusTip="Exit the Application",
                                        triggered=self.menu_close_drawing_tool)

        self.action_menu_about = QAction("About",
                                         self, statusTip="Displays info about text editor",
                                         triggered=self.menu_about)

        self.action_menu_about_qt = QAction("Help", self,
                                            statusTip="Show the Qt library's About box",
                                            triggered=self.menu_help)

    # noinspection PyAttributeOutsideInit
    def create_top_toolbar_actions(self):
        """
        Creates the buttons for the top menu, adds text and/or icons and sets the function for the user
        :return:
        """

        # Extra buttons that can trigger no matter the action selected:
        self.action_top_undo = QAction(QIcon(str(self.icon_path/"undo.png")), "Undo",
                                       self, shortcut="Ctrl+Z",
                                       statusTip="Undo/Redo changes to molecule Ctrl+Z",
                                       triggered=self.editor.undo, objectName="undo")

        # self.action_top_move = QAction(QIcon(str(self.icon_path/"move.png")), "Move",
        #                                self,
        #                                statusTip="Move single Atom around",
        #                                triggered=self.set_tool, checkable=True, objectName="Move")
        # self.actionActionGroup.addAction(self.action_top_move)

        self.action_top_clear_canvas = QAction(QIcon(str(self.icon_path/"trash.png")), "Clear Canvas",
                                               self, shortcut="Ctrl+X",
                                               statusTip="Clear Canvas (no warning)",
                                               triggered=self.toolbar_top_clear_canvas, objectName="clear_canvas")

        self.action_top_clean_structure = QAction(QIcon(str(self.icon_path/'broom.png')), 'Recalculate coordinates &F',
                                                  self, shortcut="Ctrl+F",
                                                  statusTip="Re-calculates coordinates and redraw",
                                                  triggered=self.editor.clean_up_structure,
                                                  objectName="recalculate_coordinates")

        self.action_top_show_smiles = QAction(QIcon(str(self.icon_path/"smiles.png")), "Show Smiles",
                                              self, statusTip="Shows the smiles for the current mol",
                                              triggered=self.toolbar_top_show_smiles, objectName="show_smiles")

        self.action_top_show_mol_info = QAction(QIcon(str(self.icon_path/"molecule.png")), "Show Mol Info",
                                                self, statusTip="Shows the Information for the current mol",
                                                triggered=self.toolbar_top_show_info, objectName="show_mol_info")

        self.action_top_export_mol = QAction(QIcon(str(self.icon_path/"molecule_export.png")), "Export Mol",
                                             self, statusTip="Export the current molecule",
                                             triggered=self.toolbar_top_export_mol, objectName="export_mol")

        self.action_top_import_mol = QAction(QIcon(str(self.icon_path/"molecule_import.png")), "Import Mol",
                                             self, statusTip="Import a new molecule from the smiles string",
                                             triggered=self.toolbar_top_import_mol, objectName="import_mol")

        self.smiles_field = QLineEdit()

    # noinspection PyAttributeOutsideInit
    def create_left_side_toolbar_actions(self):

        # self.action_left_side_lasso = QAction(QIcon(str(self.icon_path/"lasso.png")), 'lasso',
        #                                               self, shortcut="Ctrl+1",
        #                                               statusTip="Selection_tool",
        #                                               triggered=self.set_tool, objectName="LASSO",
        #                                               checkable=True)
        # self.actionActionGroup.addAction(self.action_left_side_lasso)
        #
        # self.action_left_side_arrow_circle = QAction(QIcon(str(self.icon_path/"arrow_circle.png")), 'arrow_circle',
        #                                               self, shortcut="Ctrl+1",
        #                                               statusTip="arrow_circle",
        #                                               triggered=self.set_tool, objectName="arrow_circle",
        #                                               checkable=True)
        # self.actionActionGroup.addAction(self.action_left_side_arrow_circle)

        self.action_left_side_single_standard_bond = QAction(QIcon(str(self.icon_path / "bond.png")), 'Add',
                                                             self, shortcut="Ctrl+1",
                                                             statusTip="Set to single Standard Bond",
                                                             triggered=self.set_tool, objectName="add",
                                                             checkable=True)
        self.actionActionGroup.addAction(self.action_left_side_single_standard_bond)

        # self.action_left_side_line_daz = QAction(QIcon(str(self.icon_path / "line_daz.png")), 'line_daz',
        #                                      self, shortcut="Ctrl+1",
        #                                      statusTip="line_daz",
        #                                      triggered=self.set_tool, objectName="line_daz",
        #                                      checkable=True)
        # self.actionActionGroup.addAction(self.action_left_side_line_daz)

        # self.action_left_side_squar_daze = QAction(QIcon(str(self.icon_path / "squar_daze.png")), 'squar_daze',
        #                                          self, shortcut="Ctrl+1",
        #                                          statusTip="squar_daze",
        #                                          triggered=self.set_tool, objectName="squar_daze",
        #                                          checkable=True)
        # self.actionActionGroup.addAction(self.action_left_side_squar_daze)

        self.action_left_side_down_bond = QAction(QIcon(str(self.icon_path / "triangle_daz.png")), 'Down Bond',
                                                  self, shortcut="Ctrl+2",
                                                  statusTip="Down Bond",
                                                  triggered=self.set_tool, objectName="down",
                                                  checkable=True)
        self.actionActionGroup.addAction(self.action_left_side_down_bond)

        self.action_left_side_up_bond = QAction(QIcon(str(self.icon_path / "triangle_filled.png")), 'Up Bond',
                                                self, shortcut="Ctrl+3",
                                                statusTip="Up Bond",
                                                triggered=self.set_tool, objectName="up",
                                                checkable=True)
        self.actionActionGroup.addAction(self.action_left_side_up_bond)

        # self.action_left_side_empty_triangle = QAction(QIcon(str(self.icon_path / "triangle_empty.png")),
        #                                                'empty_triangle',
        #                                            self, shortcut="Ctrl+4",
        #                                            statusTip="empty_triangle",
        #                                            triggered=self.set_tool, objectName="empty",
        #                                            checkable=True)
        # self.actionActionGroup.addAction(self.action_left_side_empty_triangle)

        # self.action_left_side_squvikeliqy = QAction(QIcon(str(self.icon_path / "squvikeliqy.png")), 'squvikeliqy',
        #                                            self, shortcut="Ctrl+1",
        #                                            statusTip="squvikeliqy",
        #                                            triggered=self.set_tool, objectName="squvikeliqy",
        #                                            checkable=True)
        # self.actionActionGroup.addAction(self.action_left_side_squvikeliqy)

        # self.action_left_side_multiple_bond = QAction(QIcon(str(self.icon_path / "multiple_bond.png")), 'multiple_bond',
        #                                             self, shortcut="Ctrl+1",
        #                                             statusTip="multiple_bond",
        #                                             triggered=self.set_tool, objectName="multiple_bond",
        #                                             checkable=True)
        # self.actionActionGroup.addAction(self.action_left_side_multiple_bond)

        # self.action_left_side_selection_tool = QAction(QIcon(str(self.icon_path / "selection_tool_2.png")), 'Selection tool',
        #                                             self, shortcut="Ctrl+1",
        #                                             statusTip="Selection Tool",
        #                                             triggered=self.set_tool, objectName="Selection tool",
        #                                             checkable=True)
        # self.actionActionGroup.addAction(self.action_left_side_selection_tool)

        # self.action_left_side_bond_cutter = QAction(QIcon(str(self.icon_path / "bond_cutter.png")), 'bond_cutter',
        #                                             self, shortcut="Ctrl+1",
        #                                             statusTip="bond_cutter",
        #                                             triggered=self.set_tool, objectName="bond_cutter",
        #                                             checkable=True)
        # self.actionActionGroup.addAction(self.action_left_side_bond_cutter)

        self.action_left_side_eraser = QAction(QIcon(str(self.icon_path / "eraser.png")), 'Eraser',
                                                    self, shortcut="Ctrl+E",
                                                    statusTip="Delete Atom or Bond",
                                                    triggered=self.set_tool, objectName="eraser",
                                                    checkable=True)
        self.actionActionGroup.addAction(self.action_left_side_eraser)

        # self.action_left_side_Letters = QAction("A",
        #                                             self, shortcut="Ctrl+1",
        #                                             statusTip="Letters",
        #                                             triggered=self.set_tool, objectName="Letters",
        #                                             checkable=True)

        self.action_left_side_increase_charge = QAction(QIcon(str(self.icon_path / "plus.png")), 'Increase Charge',
                                                        self, shortcut="Ctrl++",
                                                        statusTip="Increase Charge",
                                                        triggered=self.set_tool, objectName="increase_charge",
                                                        checkable=True)
        self.actionActionGroup.addAction(self.action_left_side_increase_charge)

        self.action_left_side_decrease_charge = QAction(QIcon(str(self.icon_path / "minus.png")), 'Decrease Charge',
                                                        self, shortcut="Ctrl+-",
                                                        statusTip="Decrease Charge",
                                                        triggered=self.set_tool, objectName="decrease_charge",
                                                        checkable=True)
        self.actionActionGroup.addAction(self.action_left_side_decrease_charge)

        # Standard Atoms
        self.atom_actions = []
        for atomname in ["C", "N", "O", "F"]:
            temp_atom_action = self.periodic_table.atom_actions[atomname]
            if temp_atom_action.objectName() == "C":
                temp_atom_action.setChecked(True)
                temp_atom_action.setIcon(QIcon(str(self.icon_path / "atom_c.png")))
            elif temp_atom_action.objectName() == "N":
                temp_atom_action.setIcon(QIcon(str(self.icon_path / "atom_n.png")))
            elif temp_atom_action.objectName() == "O":
                temp_atom_action.setIcon(QIcon(str(self.icon_path / "atom_o.png")))
            elif temp_atom_action.objectName() == "F":
                temp_atom_action.setIcon(QIcon(str(self.icon_path / "atom_f.png")))
            self.atom_actions.append(temp_atom_action)

        # Mol Structure buttons:
        self.action_left_side_benzene = QAction(QIcon(str(self.icon_path/"benzene.png")), "Benzene",
                                                self, shortcut="ctrl+B",
                                                statusTip="Add Benzene",
                                                triggered=self.set_tool, objectName="benzene",
                                                checkable=True)
        self.actionActionGroup.addAction(self.action_left_side_benzene)

        self.action_left_side_cyclopentane = QAction(QIcon(str(self.icon_path/"cyclopentane.png")), "Cyclopentane",
                                                self, shortcut="ctrl+P",
                                                statusTip="Add Cyclopentane",
                                                triggered=self.set_tool, objectName="cyclopentane",
                                                checkable=True)
        self.actionActionGroup.addAction(self.action_left_side_cyclopentane)


        # Misc buttons:
        self.action_left_side_open_periodic_table = QAction(QIcon(str(self.icon_path / "periodic_table.png")),
                                                            "O&pen Periodic Table",
                                                            self, shortcut="ctrl+P",
                                                            statusTip="Open the periodic table for atom type selection",
                                                            triggered=self.menu_open_periodic_table)

    def create_left_side_toolbar_action_stamps_button(self):

        stamps_menu = self.stamps_monosaccharides()
        stamps_button = QToolButton(self)
        # Create a tool button for the submenu
        stamps_button.setIcon(QIcon(str(self.icon_path / "stamp.png")))
        stamps_button.setMenu(stamps_menu)
        stamps_button.setPopupMode(QToolButton.InstantPopup)  # Set the popup mode for the tool button
        return stamps_button

    def stamps_monosaccharides(self):
        # Create a submenu for the "Stamps" button
        stamps_menu = QMenu("Stamps", self)

        # Create actions for the submenu
        stamp_glucose = QAction("Glucose",
                                self,
                                statusTip="Add Glucose",
                                triggered=self.set_tool, objectName="glucose",
                                checkable=True)
        self.actionActionGroup.addAction(stamp_glucose)

        stamp_mannose = QAction("Mannose",
                                self,
                                statusTip="Add Mannose",
                                triggered=self.set_tool, objectName="mannose",
                                checkable=True)
        self.actionActionGroup.addAction(stamp_mannose)

        stamp_galactose = QAction("Galactose",
                                self,
                                statusTip="Add Galactose",
                                triggered=self.set_tool, objectName="galactose",
                                checkable=True)
        self.actionActionGroup.addAction(stamp_galactose)

        stamp_xylose = QAction("Xylose",
                                self,
                                statusTip="Add Xylose",
                                triggered=self.set_tool, objectName="xylose",
                                checkable=True)
        self.actionActionGroup.addAction(stamp_xylose)

        stamp_arabinose = QAction("Arabinose",
                                self,
                                statusTip="Add Arabinose",
                                triggered=self.set_tool, objectName="arabinose",
                                checkable=True)
        self.actionActionGroup.addAction(stamp_arabinose)

        stamp_ribose = QAction("Ribose",
                                self,
                                statusTip="Add Ribose",
                                triggered=self.set_tool, objectName="ribose",
                                checkable=True)
        self.actionActionGroup.addAction(stamp_ribose)

        stamp_deoxyribose = QAction("Deoxyribose",
                                self,
                                statusTip="Add Deoxyribose",
                                triggered=self.set_tool, objectName="deoxyribose",
                                checkable=True)
        self.actionActionGroup.addAction(stamp_deoxyribose)

        stamp_lyxose = QAction("Lyxose",
                                self,
                                statusTip="Add Lyxose",
                                triggered=self.set_tool, objectName="lyxose",
                                checkable=True)
        self.actionActionGroup.addAction(stamp_lyxose)

        stamp_allose = QAction("Allose",
                                self,
                                statusTip="Add Allose",
                                triggered=self.set_tool, objectName="allose",
                                checkable=True)
        self.actionActionGroup.addAction(stamp_allose)

        stamp_altrose = QAction("Altrose",
                                self,
                                statusTip="Add Altrose",
                                triggered=self.set_tool, objectName="altrose",
                                checkable=True)
        self.actionActionGroup.addAction(stamp_altrose)

        stamp_talose = QAction("Talose",
                                self,
                                statusTip="Add Talose",
                                triggered=self.set_tool, objectName="talose",
                                checkable=True)
        self.actionActionGroup.addAction(stamp_talose)

        stamp_idose = QAction("Idose",
                                self,
                                statusTip="Add Idose",
                                triggered=self.set_tool, objectName="idose",
                                checkable=True)
        self.actionActionGroup.addAction(stamp_idose)

        # add actions to the menu:
        stamps_menu.addAction(stamp_glucose)
        stamps_menu.addAction(stamp_mannose)
        stamps_menu.addAction(stamp_galactose)
        stamps_menu.addAction(stamp_xylose)
        stamps_menu.addAction(stamp_arabinose)
        stamps_menu.addAction(stamp_ribose)
        stamps_menu.addAction(stamp_deoxyribose)
        stamps_menu.addAction(stamp_lyxose)
        stamps_menu.addAction(stamp_allose)
        stamps_menu.addAction(stamp_altrose)
        stamps_menu.addAction(stamp_talose)
        stamps_menu.addAction(stamp_idose)

        return stamps_menu

    @staticmethod
    def menu_close_drawing_tool():
        exit(0)

    def menu_about(self):
        QMessageBox.about(self, "Drawing and editor tool",
                          "A Simple Molecule Drawing and Editor where you can draw and edit molecules\n\n"
                          "Made By Charlie")

    def menu_help(self):
        QMessageBox.about(self, "Information and help",
                          "<b>As default:</b> Left click is adding an atom or bond, right click is "
                          "<u>deleting</u> atom or bond<br>"
                          "<b>All short-cuts</b> are listed in the status bar when hovering an icon<br>"
                          "<b>The Single bond</b> is the standard tool and can be used to add molecules and bonds<br>"
                          "<b>If the Atom selected</b> is the same as the Atom you click on, "
                          "a new atom of that type is added with a bond between the two.<br>"
                          "<b>If the selected Atom</b> is different than the one that is clicked, "
                          "it will replace it, unless you click and drag away from the atom.<br>"
                          "Then a new atom will be added and a bond between the two will be created.<br>"
                          "<b>If you click a bond</b>, it will switch between single, double and triple bond.<br>"
                          "<b>If you click an atom</b>, and drag to a different atom,"
                          " a bond between the two will be created.<br>"
                          "<b>You can add any smiles</b> code in the smile's text field, "
                          "and then click the 'import smiles'"
                          "button. This will clear the canvas and add the molecule from the smiles to the canvas.<br>"
                          "<b>Stamps can only be added</b> next to the current molecule, or on an empty canvas. <br>"
                          "You can't add stamps to a clicked atom.\n\n")

    def menu_open_periodic_table(self):
        """
        Opens the Periodic table, to make it possible to choose different atoms to add the molecule
        :return:
        """
        self.periodic_table.show()

    def toolbar_top_clear_canvas(self):
        self.editor.clear_atom_selection()
        self.editor.mol = None
        self.my_status_bar.showMessage("Canvas Cleared")

    def toolbar_top_show_smiles(self):
        temp_mol = self.editor.mol
        smiles = Chem.MolToSmiles(temp_mol)
        self.smiles_field.setText(smiles)

    def toolbar_top_show_info(self):
        temp_mol = self.editor.mol
        dialog = CustomDialog(temp_mol, self)
        dialog.show()

    def toolbar_top_export_mol(self):
        temp_mol = self.editor.mol
        self.queue_gui.put(('smiles', temp_mol))
        print(self.queue_gui)

    def toolbar_top_import_mol(self):
        self.toolbar_top_clear_canvas()
        import_mol = self.smiles_field.text()
        self.editor.mol = Chem.MolFromSmiles(import_mol)

    def set_atom_type_name(self, atom_name):
        self.editor.set_atom_type(str(atom_name))
        self.my_status_bar.showMessage(f"Atom-type {atom_name} selected")

    def set_tool(self):
        """
        Changes the selected tool through the draw_tool_editor
        :return:
        """
        sender = self.sender()
        self.editor.set_action_tool(sender.objectName())
        self.my_status_bar.showMessage(f"Action {sender.objectName()} selected")

    def set_bond_type(self):
        """
        Sets the bond type for the drawing tool.
        This should be able to set normal 2D as well as 3D bonds.
        :return:
        """
        sender = self.sender()
        self.editor.set_bond_type(sender.objectName())
        self.my_status_bar.showMessage(f"Bond type {sender.objectName()} selected")


