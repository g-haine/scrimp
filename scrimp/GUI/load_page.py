from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QComboBox, QLabel, QLineEdit, QGridLayout, QFileDialog
from utils.GUI import gui_pages, gui_width, gui_height, check_black_listed_words
import json


class Window(QtWidgets.QWidget):
    """This class defines the main page of the GUI and ask the details for the definition of the new distribuited port-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self, session):
        QtWidgets.QWidget.__init__(self)
        self.session = session
        self.setWindowTitle("Load a former GUI session")
        self.setFixedWidth(gui_width)
        self.setFixedHeight(gui_height)

        self.layout = QGridLayout()

        self.button_file_dialog = QtWidgets.QPushButton("Select File", self)
        self.button_file_dialog.clicked.connect(self.get_path)

        label_directory = QLabel('<font size="4"> The selected directory is: </font>')
        self.line_edit_directory = QLineEdit()
        self.line_edit_directory.setPlaceholderText(
            "Insert manually your filepath or click the button on the right to select the directory"
        )

        self.file_path = ""

        self.button_next = QtWidgets.QPushButton("Next >")
        self.button_next.clicked.connect(self.next_page)
        self.button_prev = QtWidgets.QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        self.layout.addWidget(label_directory, 1, 0)
        self.layout.addWidget(self.line_edit_directory, 1, 1)
        self.layout.addWidget(self.button_file_dialog, 1, 2)
        self.layout.addWidget(self.button_prev, 5, 2)
        self.layout.addWidget(self.button_next, 5, 3)

        self.setLayout(self.layout)

    def get_path(self):
        self.file_path = QFileDialog.getOpenFileName(
            self, "Open file", "c:\\", "session files (*.session)"
        )

        print(self.file_path)
        self.line_edit_directory.setText(self.file_path[0])

    def text_changed(self, page):  # s is a str
        """This function allows the navigation trhough the navigation list.
        After checking the presence of black listed words, the function hides the current page for showing the selected one.

        Args:
            page (str): the name of the page.
        """
        self.comboBox.setCurrentText("create_dphs_page")
        if not check_black_listed_words(
            self, self.line_edit_dphs_name, "Name for your dpHs"
        ):
            self.switch_window.emit(page)
            self.hide()

    def update_page(self):
        """This function manages the update of the current page."""
        self.session["cursor_on_page"] = "load_page"

    def next_page(self):
        """This function emits the signal to navigate to next page."""
        self.file_path = [self.line_edit_directory.text()]

        import os

        if os.path.isfile(self.file_path[0]):
            self.session["read_from_file"] = {"filepath": self.file_path}
            self.load_session_from_file()
            self.switch_window.emit("create_dphs_page")
            self.hide()
        # TODO alert message file not valid

    def previous_page(self):
        """This funcion emits the signal to navigate to the prvious page."""
        self.switch_window.emit("welcome_page")
        self.hide()

    def load_session_from_file(self):
        """This function load the session object from a json file."""

        with open(self.session["read_from_file"]["filepath"][0], "r") as f:
            self.session["read_from_file"]["dict"] = json.load(f)
