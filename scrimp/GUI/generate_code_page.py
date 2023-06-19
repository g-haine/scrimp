from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QFileDialog, QComboBox, QLabel, QLineEdit, QGridLayout
from utils.GUI import gui_pages


class Window(QtWidgets.QWidget):
    """This class defines the main page of the GUI and ask the details for the definition of the new distribuited port-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self):
        QtWidgets.QWidget.__init__(self)
        self.setWindowTitle(
            "Generate the Python script for the defined Discrete Port Hamiltonian System"
        )
        self.setFixedWidth(1700)
        self.setFixedHeight(600)

        layout = QGridLayout()

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

        self.button_generate = QtWidgets.QPushButton("Generate")
        self.button_generate.clicked.connect(self.generate_code)

        self.button_prev = QtWidgets.QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        layout.addWidget(label_filename, 1, 0)
        layout.addWidget(self.line_edit_filname, 1, 1)
        layout.addWidget(label_directory, 2, 0)
        layout.addWidget(self.line_edit_directory, 2, 1)
        layout.addWidget(self.button_file_dialog, 2, 2)
        layout.addWidget(self.button_generate, 3, 4)
        layout.addWidget(self.button_prev, 3, 3)

        self.comboBox = QComboBox()
        self.comboBox.addItems(gui_pages)
        self.comboBox.setCurrentText("generate_code_page")

        # There is an alternate signal to send the text.
        self.comboBox.currentTextChanged.connect(self.text_changed)
        layout.addWidget(self.comboBox, 3, 2)

        self.setLayout(layout)

    def get_path(self):
        self.file_path = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        print(self.file_path)
        self.line_edit_directory.setText(self.file_path)

    def text_changed(self, page):  # s is a str
        self.comboBox.setCurrentText("generate_code_page")
        self.switch_window.emit(page)
        self.hide()

    def generate_code(self):
        print("generate code...")

    def previous_page(self):
        """This funcion emits the signal to navigate to the prvious page."""
        self.switch_window.emit("set_time_scheme_page")
        self.hide()
