from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import (
    QHBoxLayout,
    QPushButton,
    QLineEdit,
    QGridLayout,
    QTableWidget,
    QComboBox,
)
from PyQt5.QtCore import Qt
from utils.GUI import gui_pages


class Window(QtWidgets.QWidget):
    """This class defines the add port and coport page of the GUI and asks to insert
    the ports and co-port realted to the new distribuited port-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self):
        QtWidgets.QWidget.__init__(self)

        self.setWindowTitle("Definition of Port/s")
        self.setFixedWidth(1700)
        self.setFixedHeight(600)
        # self.setGeometry(100, 100, 600, 300)

        layout = QGridLayout()

        # self.line_edit = QLineEdit()
        # layout.addWidget(self.line_edit)

        # create a QTableWidget ports
        self.table_ports = QTableWidget()
        self.table_ports.setRowCount(1)

        # adding header to the table
        header_horizontal_ports = [
            "Name",
            "Flow",
            "Effort",
            "Kind",
            "Mesh ID",
            "Algebraic",
            "Substituted",
            "Region",
        ]
        self.table_ports.setColumnCount(len(header_horizontal_ports))

        self.header_vertical_ports = ["port"]
        self.table_ports.setHorizontalHeaderLabels(header_horizontal_ports)
        self.table_ports.setVerticalHeaderLabels(self.header_vertical_ports)

        # adjust size columns of horizontal header
        for i, _ in enumerate(header_horizontal_ports):
            self.table_ports.setColumnWidth(i, 150)

        self.button_add_port = QPushButton("Add port")
        self.button_add_port.clicked.connect(self.new_port)

        self.button_delete_port = QPushButton("Remove port")
        self.button_delete_port.clicked.connect(self.delete_port)

        self.button_clear_all = QPushButton("Clear All")
        self.button_clear_all.clicked.connect(self.clear_all)

        # layout_buttons_port = QHBoxLayout()

        # layout_buttons_port.addWidget(self.button_add_port)
        # layout_buttons_port.addWidget(self.button_delete_port)

        # cell_double = QTableWidget(layout_buttons_port)

        self.button_next = QPushButton("Next >")
        self.button_next.clicked.connect(self.next_page)

        self.button_prev = QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        layout.addWidget(self.table_ports, 1, 0, 1, 3)
        # layout.addWidget(cell_double, 1, 3)
        layout.addWidget(self.button_clear_all, 0, 1)
        layout.addWidget(self.button_add_port, 0, 2, Qt.AlignTop)
        layout.addWidget(self.button_delete_port, 0, 3, Qt.AlignTop)

        layout.addWidget(self.button_next, 4, 3)
        layout.addWidget(self.button_prev, 4, 2)

        self.comboBox = QComboBox()
        self.comboBox.addItems(gui_pages)
        self.comboBox.setCurrentText("add_port_page")

        # There is an alternate signal to send the text.
        self.comboBox.currentTextChanged.connect(self.text_changed)
        layout.addWidget(self.comboBox, 4, 1)

        self.setLayout(layout)

    def text_changed(self, page):  # s is a str
        self.comboBox.setCurrentText("add_port_page")
        self.switch_window.emit(page)
        self.hide()

    def next_page(self):
        """This funciont emit the signal to navigate to the next page."""
        self.switch_window.emit("add_parameter_page")
        self.hide()

    def previous_page(self):
        """This funcion emits the signal to navigate to the prvious page."""
        self.switch_window.emit("add_state_costate_page")
        self.hide()

    def new_port(self):
        """This function adds 2 rows in the table (1 for port, 1 for co-port)"""
        count = self.table_ports.rowCount()
        self.table_ports.insertRow(count)
        self.header_vertical_ports += ["port"]
        self.table_ports.setVerticalHeaderLabels(self.header_vertical_ports)

    def delete_port(self):
        """This function removes 2 rows in the table (1 for port, 1 for co-port)"""
        if len(self.header_vertical_ports) > 1:
            self.header_vertical_ports.pop()
            self.table_ports.setVerticalHeaderLabels(self.header_vertical_ports)

            self.table_ports.removeRow(self.table_ports.rowCount() - 1)

        else:
            print("not enough element to delete!")

    def clear_all(self):
        self.table_ports.setRowCount(0)
        self.new_port()
