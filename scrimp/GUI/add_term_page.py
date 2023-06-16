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
    """This class defines the add term and coterm page of the GUI and asks to insert
    the terms and co-term realted to the new distribuited term-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self):
        QtWidgets.QWidget.__init__(self)

        self.setWindowTitle("Definition of Term/s")
        self.setFixedWidth(600)
        self.setFixedHeight(600)
        # self.setGeometry(100, 100, 600, 300)

        layout = QGridLayout()

        # self.line_edit = QLineEdit()
        # layout.addWidget(self.line_edit)

        # create a QTableWidget terms
        self.table_terms = QTableWidget()
        self.table_terms.setRowCount(1)

        # adding header to the table
        header_horizontal_terms = ["Description", "Expression", "Regions", "Mesh ID"]
        self.table_terms.setColumnCount(len(header_horizontal_terms))

        self.header_vertical_terms = ["term"]
        self.table_terms.setHorizontalHeaderLabels(header_horizontal_terms)
        self.table_terms.setVerticalHeaderLabels(self.header_vertical_terms)

        # adjust size columns of horizontal header
        for i, _ in enumerate(header_horizontal_terms):
            self.table_terms.setColumnWidth(i, 130)

        self.button_add_term = QPushButton("Add term")
        self.button_add_term.clicked.connect(self.new_term)

        self.button_delete_term = QPushButton("Remove term")
        self.button_delete_term.clicked.connect(self.delete_term)

        self.button_clear_all = QPushButton("Clear All")
        self.button_clear_all.clicked.connect(self.clear_all)

        # layout_buttons_term = QHBoxLayout()

        # layout_buttons_term.addWidget(self.button_add_term)
        # layout_buttons_term.addWidget(self.button_delete_term)

        # cell_double = QTableWidget(layout_buttons_term)

        self.button_next = QPushButton("Next >")
        self.button_next.clicked.connect(self.next_page)

        self.button_prev = QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        layout.addWidget(self.table_terms, 1, 0, 1, 3)
        # layout.addWidget(cell_double, 1, 3)
        layout.addWidget(self.button_clear_all, 0, 1)
        layout.addWidget(self.button_add_term, 0, 2, Qt.AlignTop)
        layout.addWidget(self.button_delete_term, 0, 3, Qt.AlignTop)

        layout.addWidget(self.button_next, 4, 3)
        layout.addWidget(self.button_prev, 4, 2)

        self.setLayout(layout)

    def next_page(self):
        """This funciont emit the signal to navigate to the next page."""
        self.switch_window.emit("add_brick_page")
        self.hide()

    def previous_page(self):
        """This funcion emits the signal to navigate to the prvious page."""
        self.switch_window.emit("set_hamiltonian_page")
        self.hide()

    def new_term(self):
        """This function adds 2 rows in the table (1 for term, 1 for co-term)"""
        count = self.table_terms.rowCount()
        self.table_terms.insertRow(count)
        self.header_vertical_terms += ["term"]
        self.table_terms.setVerticalHeaderLabels(self.header_vertical_terms)

    def delete_term(self):
        """This function removes 2 rows in the table (1 for term, 1 for co-term)"""
        if len(self.header_vertical_terms) > 1:
            self.header_vertical_terms.pop()
            self.table_terms.setVerticalHeaderLabels(self.header_vertical_terms)

            self.table_terms.removeRow(self.table_terms.rowCount() - 1)

        else:
            print("not enough element to delete!")

    def clear_all(self):
        self.table_terms.setRowCount(0)
        self.new_term()
