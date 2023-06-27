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
from utils.GUI import gui_pages, gui_width, gui_height


class Window(QtWidgets.QWidget):
    """This class defines the add control_port and cocontrol_port page of the GUI and asks to insert
    the control_ports and co-control_port realted to the new distribuited control_port-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self):
        QtWidgets.QWidget.__init__(self)

        self.setWindowTitle("Definition of Control control_port/s")
        self.setFixedWidth(gui_width)
        self.setFixedHeight(gui_height)
        # self.setGeometry(100, 100, 600, 300)

        layout = QGridLayout()

        # self.line_edit = QLineEdit()
        # layout.addWidget(self.line_edit)

        # create a QTableWidget control_ports
        self.table_control_ports = QTableWidget()
        self.table_control_ports.setRowCount(1)

        # adding header to the table
        header_horizontal_control_ports = [
            "Name",
            "Name Control",
            "Description Control",
            "Name Observation",
            "Description Observation",
            "Kind",
            "Region",
            "Mesh ID",
        ]
        self.table_control_ports.setColumnCount(len(header_horizontal_control_ports))

        self.header_vertical_control_ports = ["control_port"]
        self.table_control_ports.setHorizontalHeaderLabels(
            header_horizontal_control_ports
        )
        self.table_control_ports.setVerticalHeaderLabels(
            self.header_vertical_control_ports
        )

        # adjust size columns of horizontal header
        for i, _ in enumerate(header_horizontal_control_ports):
            self.table_control_ports.setColumnWidth(i, 155)

        self.button_add_control_port = QPushButton("Add control_port")
        self.button_add_control_port.clicked.connect(self.new_control_port)

        self.button_delete_control_port = QPushButton("Remove control_port")
        self.button_delete_control_port.clicked.connect(self.delete_control_port)

        self.button_clear_all = QPushButton("Clear All")
        self.button_clear_all.clicked.connect(self.clear_all)

        # layout_buttons_control_port = QHBoxLayout()

        # layout_buttons_control_port.addWidget(self.button_add_control_port)
        # layout_buttons_control_port.addWidget(self.button_delete_control_port)

        # cell_double = QTableWidget(layout_buttons_control_port)

        self.button_next = QPushButton("Next >")
        self.button_next.clicked.connect(self.next_page)

        self.button_prev = QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        layout.addWidget(self.table_control_ports, 1, 0, 1, 3)
        # layout.addWidget(cell_double, 1, 3)
        layout.addWidget(self.button_clear_all, 0, 1)
        layout.addWidget(self.button_add_control_port, 0, 2, Qt.AlignTop)
        layout.addWidget(self.button_delete_control_port, 0, 3, Qt.AlignTop)

        layout.addWidget(self.button_next, 4, 3)
        layout.addWidget(self.button_prev, 4, 2)

        # create navigation list
        self.comboBox = QComboBox()
        self.comboBox.addItems(gui_pages)
        self.comboBox.setCurrentText("add_control_port_page")

        # There is an alternate signal to send the text.
        self.comboBox.currentTextChanged.connect(self.text_changed)
        layout.addWidget(self.comboBox, 4, 1)

        self.setLayout(layout)

    def text_changed(self, page):  # s is a str
        self.comboBox.setCurrentText("add_control_port_page")
        self.switch_window.emit(page)
        self.hide()

    def next_page(self):
        """This funciont emit the signal to navigate to the next page."""
        self.switch_window.emit("add_fem_page")
        self.hide()

    def previous_page(self):
        """This funcion emits the signal to navigate to the prvious page."""
        self.switch_window.emit("add_parameter_page")
        self.hide()

    def new_control_port(self):
        """This function adds 2 rows in the table (1 for control_port, 1 for co-control_port)"""
        count = self.table_control_ports.rowCount()
        self.table_control_ports.insertRow(count)
        self.header_vertical_control_ports += ["control_port"]
        self.table_control_ports.setVerticalHeaderLabels(
            self.header_vertical_control_ports
        )

    def delete_control_port(self):
        """This function removes 2 rows in the table (1 for control_port, 1 for co-control_port)"""
        if len(self.header_vertical_control_ports) > 1:
            self.header_vertical_control_ports.pop()
            self.table_control_ports.setVerticalHeaderLabels(
                self.header_vertical_control_ports
            )

            self.table_control_ports.removeRow(self.table_control_ports.rowCount() - 1)

        else:
            print("not enough element to delete!")

    def clear_all(self):
        self.table_control_ports.setRowCount(0)
        self.new_control_port()
