from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QComboBox, QLabel, QLineEdit, QGridLayout
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

        label_filename = QLabel('<font size="4"> Name for your dpHs: </font>')
        self.line_edit_filname = QLineEdit()
        self.line_edit_filname.setPlaceholderText("Please enter the name")

        self.button_generate = QtWidgets.QPushButton("Generate")
        self.button_generate.clicked.connect(self.generate_code)

        self.button_prev = QtWidgets.QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        layout.addWidget(label_filename, 1, 0)
        layout.addWidget(self.line_edit_filname, 1, 1)
        layout.addWidget(self.button_generate, 3, 4)
        layout.addWidget(self.button_prev, 3, 3)

        self.comboBox = QComboBox()
        self.comboBox.addItems(gui_pages)
        self.comboBox.setCurrentText("generate_code_page")

        # There is an alternate signal to send the text.
        self.comboBox.currentTextChanged.connect(self.text_changed)
        layout.addWidget(self.comboBox, 3, 2)

        self.setLayout(layout)

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
