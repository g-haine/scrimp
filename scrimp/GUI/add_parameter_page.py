from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import (
    QHBoxLayout,
    QPushButton,
    QLineEdit,
    QGridLayout,
    QTableWidget,
)
from PyQt5.QtCore import Qt


class Window(QtWidgets.QWidget):
    """This class defines the add parameter and coparameter page of the GUI and asks to insert
    the parameters and co-parameter realted to the new distribuited parameter-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self):
        QtWidgets.QWidget.__init__(self)

        self.setWindowTitle("Definition of Parameter/s")
        self.setFixedWidth(900)
        self.setFixedHeight(600)
        # self.setGeometry(100, 100, 600, 300)

        layout = QGridLayout()

        # self.line_edit = QLineEdit()
        # layout.addWidget(self.line_edit)

        # create a QTableWidget parameters
        self.table_parameters = QTableWidget()
        self.table_parameters.setRowCount(1)

        # adding header to the table
        header_horizontal_parameters = [
            "Name",
            "Description",
            "Kind",
            "Expression",
            "Name parameter",
        ]
        self.table_parameters.setColumnCount(len(header_horizontal_parameters))

        self.header_vertical_parameters = ["parameter"]
        self.table_parameters.setHorizontalHeaderLabels(header_horizontal_parameters)
        self.table_parameters.setVerticalHeaderLabels(self.header_vertical_parameters)

        # adjust size columns of horizontal header
        for i, _ in enumerate(header_horizontal_parameters):
            self.table_parameters.setColumnWidth(i, 115)

        self.button_add_parameter = QPushButton("Add parameter")
        self.button_add_parameter.clicked.connect(self.new_parameter)

        self.button_delete_parameter = QPushButton("Remove parameter")
        self.button_delete_parameter.clicked.connect(self.delete_parameter)

        self.button_clear_all = QPushButton("Clear All")
        self.button_clear_all.clicked.connect(self.clear_all)

        # layout_buttons_parameter = QHBoxLayout()

        # layout_buttons_parameter.addWidget(self.button_add_parameter)
        # layout_buttons_parameter.addWidget(self.button_delete_parameter)

        # cell_double = QTableWidget(layout_buttons_parameter)

        self.button_next = QPushButton("Next >")
        self.button_next.clicked.connect(self.next_page)

        self.button_prev = QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        layout.addWidget(self.table_parameters, 1, 0, 1, 3)
        # layout.addWidget(cell_double, 1, 3)
        layout.addWidget(self.button_clear_all, 0, 1)
        layout.addWidget(self.button_add_parameter, 0, 2, Qt.AlignTop)
        layout.addWidget(self.button_delete_parameter, 0, 3, Qt.AlignTop)

        layout.addWidget(self.button_next, 4, 3)
        layout.addWidget(self.button_prev, 4, 2)

        self.setLayout(layout)

    def next_page(self):
        """This funciont emit the signal to navigate to the next page."""
        self.switch_window.emit("add_control_port_page")
        self.hide()

    def previous_page(self):
        """This funcion emits the signal to navigate to the prvious page."""
        self.switch_window.emit("add_port_page")
        self.hide()

    def new_parameter(self):
        """This function adds 2 rows in the table (1 for parameter, 1 for co-parameter)"""
        count = self.table_parameters.rowCount()
        self.table_parameters.insertRow(count)
        self.header_vertical_parameters += ["parameter"]
        self.table_parameters.setVerticalHeaderLabels(self.header_vertical_parameters)

    def delete_parameter(self):
        """This function removes 2 rows in the table (1 for parameter, 1 for co-parameter)"""
        if len(self.header_vertical_parameters) > 1:
            self.header_vertical_parameters.pop()
            self.table_parameters.setVerticalHeaderLabels(
                self.header_vertical_parameters
            )

            self.table_parameters.removeRow(self.table_parameters.rowCount() - 1)

        else:
            print("not enough element to delete!")

    def clear_all(self):
        self.table_parameters.setRowCount(0)
        self.new_parameter()
