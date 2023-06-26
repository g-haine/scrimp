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
        self.setWindowTitle("Definition of the Hamiltonian: ")
        self.setFixedWidth(1700)
        self.setFixedHeight(600)

        layout = QGridLayout()

        label_hamiltonian_name = QLabel(
            '<font size="4"> Name for the Hamiltonian: </font>'
        )
        self.line_edit_hamiltonian_name = QLineEdit()
        self.line_edit_hamiltonian_name.setPlaceholderText("Please enter the name")

        self.button_next = QtWidgets.QPushButton("Next >")
        self.button_next.clicked.connect(self.next_page)

        self.button_prev = QtWidgets.QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        layout.addWidget(label_hamiltonian_name, 1, 0)
        layout.addWidget(self.line_edit_hamiltonian_name, 1, 1)

        layout.addWidget(self.button_next, 4, 4)
        layout.addWidget(self.button_prev, 4, 3)

        # create navigation list
        self.comboBox = QComboBox()
        self.comboBox.addItems(gui_pages)
        self.comboBox.setCurrentText("set_hamiltonian_page")

        # There is an alternate signal to send the text.
        self.comboBox.currentTextChanged.connect(self.text_changed)
        layout.addWidget(self.comboBox, 4, 2)

        self.setLayout(layout)

    def text_changed(self, page):  # s is a str
        self.comboBox.setCurrentText("set_hamiltonian_page")
        self.switch_window.emit(page)
        self.hide()

    def next_page(self):
        """This funciont emit the signal to navigate to the next page."""
        self.switch_window.emit("add_term_page")
        self.hide()

    def previous_page(self):
        """This funcion emits the signal to navigate to the prvious page."""
        self.switch_window.emit("add_fem_page")
        self.hide()
