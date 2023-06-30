from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QComboBox, QLabel, QLineEdit, QGridLayout
from utils.GUI import gui_pages, gui_width, gui_height


class Window(QtWidgets.QWidget):
    """This class defines the main page of the GUI and ask the details for the definition of the new distribuited port-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self):
        QtWidgets.QWidget.__init__(self)
        self.setWindowTitle("Definition of new distributed port-Hamiltonian system")
        self.setFixedWidth(gui_width)
        self.setFixedHeight(gui_height)

        layout = QGridLayout()

        label_dphs_name = QLabel('<font size="4"> Name for your dpHs: </font>')
        self.line_edit_dphs_name = QLineEdit()
        self.line_edit_dphs_name.setPlaceholderText(
            "Please enter the name of your Discrete Port Hamiltonian System."
        )

        linlabel_dphs_type = QLabel('<font size="4"> Type of dpHS: </font>')
        self.comboBox_dphs_type = QComboBox()
        self.comboBox_dphs_type.addItems(["real", "complex"])

        self.button_next = QtWidgets.QPushButton("Next >")
        self.button_next.clicked.connect(self.next_page)

        layout.addWidget(label_dphs_name, 1, 0)
        layout.addWidget(self.line_edit_dphs_name, 1, 1)
        layout.addWidget(linlabel_dphs_type, 2, 0)
        layout.addWidget(self.comboBox_dphs_type, 2, 1)
        layout.addWidget(self.button_next, 3, 3)

        # create navigation list
        self.comboBox = QComboBox()
        self.comboBox.addItems(gui_pages)
        self.comboBox.setCurrentText("create_dphs_page")

        # There is an alternate signal to send the text.
        self.comboBox.currentTextChanged.connect(self.text_changed)
        layout.addWidget(self.comboBox, 3, 2)

        self.setLayout(layout)

    def text_changed(self, page):  # s is a str
        self.comboBox.setCurrentText("create_dphs_page")
        self.switch_window.emit(page)
        self.hide()

    def next_page(self):
        """This funciont emit the signal to navigate to next page."""
        self.switch_window.emit("set_domain_page")
        self.hide()