from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import (
    QComboBox,
    QLabel,
    QLineEdit,
    QGridLayout,
    QFileDialog,
    QCheckBox,
)
from utils.GUI import gui_pages, gui_width, gui_height, check_black_listed_words


class Window(QtWidgets.QWidget):
    """This class defines the main page of the GUI and ask the details for the definition of the new distribuited port-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self, session):
        QtWidgets.QWidget.__init__(self)
        self.session = session
        self.session["cursor_on_page"] = "welcome_page"
        self.setWindowTitle("Welcome in the GUI of SCRIMP")
        self.setFixedWidth(gui_width)
        self.setFixedHeight(gui_height)

        self.layout = QGridLayout()

        self.button_create = QtWidgets.QPushButton("Create NEW Script", self)
        self.button_create.clicked.connect(self.go_to_create_page)
        self.button_load = QtWidgets.QPushButton("Load a Script")
        self.button_load.clicked.connect(self.go_to_load_page)

        self.layout.addWidget(self.button_create, 1, 0)
        self.layout.addWidget(self.button_load, 1, 1)

        self.setLayout(self.layout)

    def update_page(self):
        """This function manages the update of the current page."""
        self.session["cursor_on_page"] = "welcome_page"

    def go_to_load_page(self):
        """This function emits the signal to navigate to next page."""
        self.switch_window.emit("load_page")
        self.hide()

    def go_to_create_page(self):
        """This function emits the signal to navigate to next page."""
        self.switch_window.emit("create_dphs_page")
        self.hide()