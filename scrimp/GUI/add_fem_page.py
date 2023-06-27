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
    """This class defines the add FEM and coFEM page of the GUI and asks to insert
    the FEMs and co-FEM realted to the new distribuited FEM-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self):
        QtWidgets.QWidget.__init__(self)

        self.setWindowTitle("Definition of FEM/s")
        self.setFixedWidth(gui_width)
        self.setFixedHeight(gui_height)
        # self.setGeometry(100, 100, 600, 300)

        layout = QGridLayout()

        # self.line_edit = QLineEdit()
        # layout.addWidget(self.line_edit)

        # create a QTableWidget FEMs
        self.table_FEMs = QTableWidget()
        # self.table_FEMs.setRowCount(1)

        # adding header to the table
        header_horizontal_FEMs = ["Name", "Order", "FEM"]
        self.table_FEMs.setColumnCount(len(header_horizontal_FEMs))

        self.header_vertical_FEMs = ["FEM"]
        self.table_FEMs.setHorizontalHeaderLabels(header_horizontal_FEMs)
        self.table_FEMs.setVerticalHeaderLabels(self.header_vertical_FEMs)

        # adjust size columns of horizontal header
        for i, _ in enumerate(header_horizontal_FEMs):
            self.table_FEMs.setColumnWidth(i, 150)

        self.button_add_FEM = QPushButton("Add FEM")
        self.button_add_FEM.clicked.connect(self.new_FEM)

        self.button_delete_FEM = QPushButton("Remove FEM")
        self.button_delete_FEM.clicked.connect(self.delete_FEM)

        self.button_clear_all = QPushButton("Clear All")
        self.button_clear_all.clicked.connect(self.clear_all)

        # layout_buttons_FEM = QHBoxLayout()

        # layout_buttons_FEM.addWidget(self.button_add_FEM)
        # layout_buttons_FEM.addWidget(self.button_delete_FEM)

        # cell_double = QTableWidget(layout_buttons_FEM)

        self.button_next = QPushButton("Next >")
        self.button_next.clicked.connect(self.next_page)

        self.button_prev = QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        layout.addWidget(self.table_FEMs, 1, 0, 1, 3)
        # layout.addWidget(cell_double, 1, 3)
        layout.addWidget(self.button_clear_all, 0, 1)
        layout.addWidget(self.button_add_FEM, 0, 2, Qt.AlignTop)
        layout.addWidget(self.button_delete_FEM, 0, 3, Qt.AlignTop)

        layout.addWidget(self.button_next, 4, 3)
        layout.addWidget(self.button_prev, 4, 2)

        # create navigation list
        self.comboBox = QComboBox()
        self.comboBox.addItems(gui_pages)
        self.comboBox.setCurrentText("add_fem_page")

        # There is an alternate signal to send the text.
        self.comboBox.currentTextChanged.connect(self.text_changed)
        layout.addWidget(self.comboBox, 4, 1)

        self.setLayout(layout)
        self.new_FEM()

    def text_changed(self, page):  # s is a str
        self.comboBox.setCurrentText("add_fem_page")
        self.switch_window.emit(page)
        self.hide()

    def next_page(self):
        """This funciont emit the signal to navigate to the next page."""
        self.switch_window.emit("set_hamiltonian_page")
        self.hide()

    def previous_page(self):
        """This funcion emits the signal to navigate to the prvious page."""
        self.switch_window.emit("add_control_port_page")
        self.hide()

    def new_FEM(self):
        """This function adds 2 rows in the table (1 for FEM, 1 for co-FEM)"""
        count = self.table_FEMs.rowCount()
        self.table_FEMs.insertRow(count)
        self.header_vertical_FEMs += ["FEM"]
        self.table_FEMs.setVerticalHeaderLabels(self.header_vertical_FEMs)
        # create table to add in cell of table
        fem_choice = QComboBox()
        fem_choice.addItems(["CG", "DG"])
        self.table_FEMs.setCellWidget(count, 2, fem_choice)

    def delete_FEM(self):
        """This function removes 2 rows in the table (1 for FEM, 1 for co-FEM)"""
        if len(self.header_vertical_FEMs) > 1:
            self.header_vertical_FEMs.pop()
            self.table_FEMs.setVerticalHeaderLabels(self.header_vertical_FEMs)

            self.table_FEMs.removeRow(self.table_FEMs.rowCount() - 1)

        else:
            print("not enough element to delete!")

    def clear_all(self):
        self.table_FEMs.setRowCount(0)
        self.new_FEM()
