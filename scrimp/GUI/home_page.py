from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QLabel, QLineEdit, QGridLayout


class Window(QtWidgets.QWidget):
    """This class defines the main page of the GUI and ask the details for the definition of the new distribuited port-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self):
        QtWidgets.QWidget.__init__(self)
        self.setWindowTitle("Definition of new distributed port-Hamiltonian system")
        self.setFixedWidth(600)
        self.setFixedHeight(250)

        layout = QGridLayout()

        label_dphs_name = QLabel('<font size="4"> Name for your dpHs</font>')
        self.line_edit_dphs_name = QLineEdit()
        self.line_edit_dphs_name.setPlaceholderText("Please enter the name")
        layout.addWidget(label_dphs_name, 1, 0)
        layout.addWidget(self.line_edit_dphs_name, 1, 1)

        linlabel_dphs_type = QLabel('<font size="4"> Type of dpHS </font>')
        self.line_edit_dphs_type = QLineEdit()
        self.line_edit_dphs_type.setPlaceholderText(
            "Please enter the type of the system"
        )
        layout.addWidget(linlabel_dphs_type, 2, 0)
        layout.addWidget(self.line_edit_dphs_type, 2, 1)

        self.button_next = QtWidgets.QPushButton("Next >")
        self.button_next.clicked.connect(self.next_page)

        layout.addWidget(self.button_next, 3, 2)

        self.setLayout(layout)

    def next_page(self):
        """This funciont emit the signal to navigate to next page."""
        self.switch_window.emit("set_domain_page")
        self.hide()
