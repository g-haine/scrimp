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
        self.setWindowTitle("Definition of new distributed port-Hamiltonian system")
        self.setFixedWidth(gui_width)
        self.setFixedHeight(gui_height)

        self.layout = QGridLayout()

        label_filename = QLabel(
            '<font size="4"> Name for your script based on SCRIMP: </font>'
        )
        self.line_edit_filname = QLineEdit()
        self.line_edit_filname.setPlaceholderText(
            "Please insert the name for your script"
        )

        self.button_file_dialog = QtWidgets.QPushButton("Select Folder", self)
        self.button_file_dialog.clicked.connect(self.get_path)

        label_directory = QLabel('<font size="4"> The selected directory is: </font>')
        self.line_edit_directory = QLineEdit()
        self.line_edit_directory.setPlaceholderText(
            "Insert manually your filepath or click the button on the right to select the directory"
        )

        self.file_path = ""

        label_dphs_name = QLabel('<font size="4"> Name for your dpHs: </font>')
        self.line_edit_dphs_name = QLineEdit()
        self.line_edit_dphs_name.setPlaceholderText(
            "Please enter the name of your Discrete Port Hamiltonian System."
        )

        self.label_question = QLabel(
            "!! Check to enable the auto save mode to store your session on the GUI !!:"
        )

        self.checkBox_answer = QCheckBox()
        self.checkBox_answer.setChecked(self.session["auto_save"])
        self.checkBox_answer.toggled.connect(self.check_state)

        linlabel_dphs_type = QLabel('<font size="4"> Type of dpHS: </font>')
        self.comboBox_dphs_type = QComboBox()
        self.comboBox_dphs_type.addItems(["real", "complex"])

        self.button_next = QtWidgets.QPushButton("Next >")
        self.button_next.clicked.connect(self.next_page)
        self.button_prev = QtWidgets.QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        self.layout.addWidget(label_filename, 1, 0)
        self.layout.addWidget(self.line_edit_filname, 1, 1)
        self.layout.addWidget(label_directory, 2, 0)
        self.layout.addWidget(self.line_edit_directory, 2, 1)
        self.layout.addWidget(self.button_file_dialog, 2, 2)

        self.layout.addWidget(label_dphs_name, 3, 0)
        self.layout.addWidget(self.line_edit_dphs_name, 3, 1)
        self.layout.addWidget(linlabel_dphs_type, 4, 0)
        self.layout.addWidget(self.comboBox_dphs_type, 4, 1)
        self.layout.addWidget(self.label_question, 5, 0)
        self.layout.addWidget(self.checkBox_answer, 5, 1)
        self.layout.addWidget(self.button_prev, 6, 2)
        self.layout.addWidget(self.button_next, 6, 3)

        # create navigation list
        self.comboBox = QComboBox()
        self.comboBox.addItems(gui_pages)
        self.comboBox.setCurrentText("create_dphs_page")

        # There is an alternate signal to send the text.
        self.comboBox.currentTextChanged.connect(self.text_changed)
        self.layout.addWidget(self.comboBox, 6, 1)

        self.setLayout(self.layout)

    def check_state(self):
        """This function checks wether or not the checkbox has been checked and update the page accordingly."""
        print("clicked")
        if not self.checkBox_answer.isChecked():
            self.session["auto_save"] = False
        else:
            self.session["auto_save"] = True

    def get_path(self):
        self.file_path = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        print(self.file_path)
        self.line_edit_directory.setText(self.file_path)

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
            self.update_session()
            self.switch_window.emit(page)
            self.hide()

    def update_page(self):
        """This function manages the update of the current page."""
        self.line_edit_filname.setText(self.session["filename"])
        self.line_edit_directory.setText(self.session["filepath"])

    def update_session(self):
        self.session["filename"] = self.line_edit_filname.text()
        self.session["filepath"] = self.line_edit_directory.text()

    def next_page(self):
        """This function emits the signal to navigate to next page."""
        if not check_black_listed_words(
            self, self.line_edit_dphs_name, "Name for your dpHs"
        ):
            self.update_session()
            self.switch_window.emit("set_domain_page")
            self.hide()

    def previous_page(self):
        """This funcion emits the signal to navigate to the prvious page."""
        if self.session["cursor_on_page"] == "load_page":
            self.switch_window.emit("load_page")
        else:
            self.switch_window.emit("welcome_page")
        self.hide()
