from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import (
    QPushButton,
    QLabel,
    QGridLayout,
)


class Window(QtWidgets.QWidget):
    """This class defines the add port page of the GUI and asks to insert the port for the states and costate of the new distribuited port-Hamiltonian system.


    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self):
        QtWidgets.QWidget.__init__(self)

        self.setWindowTitle("Window Two")

        layout = QGridLayout()

        self.label = QLabel()
        layout.addWidget(self.label)

        self.button = QPushButton("Close")
        self.button.clicked.connect(self.close)

        layout.addWidget(self.button)

        self.button_next = QPushButton("Next >")
        self.button_next.clicked.connect(self.next_page)

        layout.addWidget(self.button_next, 3, 2)

        self.button_prev = QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        layout.addWidget(self.button_prev, 3, 1)

        self.setLayout(layout)

    def next_page(self):
        """This funciont emit the signal to navigate to the next page."""
        self.switch_window.emit("next_page")
        self.hide()

    def previous_page(self):
        """This funciont emit the signal to navigate to the prvious page."""
        self.switch_window.emit("add_state_costate_page")
        self.hide()
