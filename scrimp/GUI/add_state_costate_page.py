from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QPushButton, QLineEdit, QGridLayout, QTableWidget


class Window(QtWidgets.QWidget):
    """This class defines the add state and costate page of the GUI and asks to insert
    the states and co-state realted to the new distribuited port-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self):
        QtWidgets.QWidget.__init__(self)

        self.setWindowTitle("Definition of State/s and Co-state/s")
        self.setGeometry(100, 100, 600, 300)

        layout = QGridLayout()

        # self.line_edit = QLineEdit()
        # layout.addWidget(self.line_edit)

        # create a QTableWidget
        self.table = QTableWidget()
        self.table.setRowCount(2)
        self.table.setColumnCount(3)
        self.table.setGeometry(50, 100, 300, 300)

        # adding header to the table
        header_horizontal = ["Name", "Description", "Kind"]
        self.header_vertical = ["state", "co-state"]
        self.table.setHorizontalHeaderLabels(header_horizontal)
        self.table.setVerticalHeaderLabels(self.header_vertical)

        # adjust size columns of horizontal header
        for i, _ in enumerate(header_horizontal):
            self.table.setColumnWidth(i, 115)

        # adjust size columns of horizontal header
        self.table.setColumnWidth(i, 140)

        layout.addWidget(self.table, 1, 0, 1, 3)

        self.button_next = QPushButton("New")
        self.button_next.clicked.connect(self.new_rows)

        layout.addWidget(self.button_next, 0, 2)

        self.button_prev = QPushButton("Delete")
        self.button_prev.clicked.connect(self.delete_rows)

        layout.addWidget(self.button_prev, 0, 3)

        self.button_next = QPushButton("Next >")
        self.button_next.clicked.connect(self.next_page)

        layout.addWidget(self.button_next, 3, 3)

        self.button_prev = QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        layout.addWidget(self.button_prev, 3, 2)

        self.setLayout(layout)

    def next_page(self):
        """This funciont emit the signal to navigate to the next page."""
        self.switch_window.emit("add_port_page")
        self.hide()

    def previous_page(self):
        """This funcion emits the signal to navigate to the prvious page."""
        self.switch_window.emit("set_domain_page")
        self.hide()

    def new_rows(self):
        """This function adds 2 rows in the table (1 for state, 1 for co-state)"""
        count = self.table.rowCount()
        for _ in range(2):
            self.table.insertRow(count)
        self.header_vertical += ["state", "co-state"]
        self.table.setVerticalHeaderLabels(self.header_vertical)

    def delete_rows(self):
        """This function removes 2 rows in the table (1 for state, 1 for co-state)"""
        if len(self.header_vertical) > 2:
            self.header_vertical.pop()
            self.header_vertical.pop()
            self.table.setVerticalHeaderLabels(self.header_vertical)

            for _ in range(2):
                self.table.removeRow(self.table.rowCount() - 1)
        else:
            print("not enough element to delete!")
