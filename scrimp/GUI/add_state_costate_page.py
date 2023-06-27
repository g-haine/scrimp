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
    """This class defines the add state and costate page of the GUI and asks to insert
    the states and co-state realted to the new distribuited port-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self):
        QtWidgets.QWidget.__init__(self)

        self.setWindowTitle("Definition of State/s and Costate/s")
        self.setFixedWidth(gui_width)
        self.setFixedHeight(gui_height)
        # self.setGeometry(100, 100, 600, 300)

        layout = QGridLayout()

        # self.line_edit = QLineEdit()
        # layout.addWidget(self.line_edit)

        # create a QTableWidget States
        self.table_states = QTableWidget()
        self.table_states.setRowCount(1)
        self.table_states.setColumnCount(5)
        # self.table_states.setGeometry(50, 100, 300, 300)

        # adding header to the table
        header_horizontal_states = ["Name", "Description", "Kind", "Region", "Mesh ID"]
        self.header_vertical_states = ["state"]
        self.table_states.setHorizontalHeaderLabels(header_horizontal_states)
        self.table_states.setVerticalHeaderLabels(self.header_vertical_states)

        # adjust size columns of horizontal header
        for i, _ in enumerate(header_horizontal_states):
            self.table_states.setColumnWidth(i, 150)

        self.button_add_state = QPushButton("Add State/Costate")
        self.button_add_state.clicked.connect(self.new_state)

        self.button_delete_state = QPushButton("Remove State/Costate")
        self.button_delete_state.clicked.connect(self.delete_state)

        self.button_clear_all = QPushButton("Clear All")
        self.button_clear_all.clicked.connect(self.clear_all)

        # layout_buttons_state = QHBoxLayout()

        # layout_buttons_state.addWidget(self.button_add_state)
        # layout_buttons_state.addWidget(self.button_delete_state)

        # cell_double = QTableWidget(layout_buttons_state)

        # create a QTableWidget Co-States
        self.table_costates = QTableWidget()
        self.table_costates.setRowCount(1)
        self.table_costates.setColumnCount(5)
        # self.table_costates.setGeometry(50, 100, 300, 300)

        # adding header to the table
        header_horizontal_costates = [
            "Name",
            "Description",
            "State",
            "Substituted",
            "Mesh ID",
        ]
        self.header_vertical_costates = ["costate"]
        self.table_costates.setHorizontalHeaderLabels(header_horizontal_costates)
        self.table_costates.setVerticalHeaderLabels(self.header_vertical_costates)

        # adjust size columns of horizontal header
        for i, _ in enumerate(header_horizontal_costates):
            self.table_costates.setColumnWidth(i, 150)

        # self.button_add_costate = QPushButton("Add Costate")
        # self.button_add_costate.clicked.connect(self.new_costate)

        # self.button_delete_costate = QPushButton("Delete Costate")
        # self.button_delete_costate.clicked.connect(self.delete_costate)

        self.button_next = QPushButton("Next >")
        self.button_next.clicked.connect(self.next_page)

        self.button_prev = QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        layout.addWidget(self.table_states, 1, 0, 1, 3)
        # layout.addWidget(cell_double, 1, 3)
        layout.addWidget(self.button_clear_all, 0, 1)
        layout.addWidget(self.button_add_state, 0, 2, Qt.AlignTop)
        layout.addWidget(self.button_delete_state, 0, 3, Qt.AlignTop)
        layout.addWidget(self.table_costates, 3, 0, 1, 3)
        # layout.addWidget(self.button_add_costate, 2, 2)
        # layout.addWidget(self.button_delete_costate, 2, 3)
        layout.addWidget(self.button_next, 4, 3)
        layout.addWidget(self.button_prev, 4, 2)

        # create navigation list
        self.comboBox = QComboBox()
        self.comboBox.addItems(gui_pages)
        self.comboBox.setCurrentText("add_state_costate_page")

        # There is an alternate signal to send the text.
        self.comboBox.currentTextChanged.connect(self.text_changed)
        layout.addWidget(self.comboBox, 4, 1)

        self.setLayout(layout)

    def text_changed(self, page):  # s is a str
        self.comboBox.setCurrentText("add_state_costate_page")
        self.switch_window.emit(page)
        self.hide()

    def next_page(self):
        """This funciont emit the signal to navigate to the next page."""
        self.switch_window.emit("add_port_page")
        self.hide()

    def previous_page(self):
        """This funcion emits the signal to navigate to the prvious page."""
        self.switch_window.emit("set_domain_page")
        self.hide()

    def new_state(self):
        """This function adds 2 rows in the table (1 for state, 1 for co-state)"""
        count = self.table_states.rowCount()
        self.table_states.insertRow(count)
        self.header_vertical_states += ["state"]
        self.table_states.setVerticalHeaderLabels(self.header_vertical_states)
        self.new_costate()

    def delete_state(self):
        """This function removes 2 rows in the table (1 for state, 1 for co-state)"""
        if len(self.header_vertical_states) > 1:
            self.header_vertical_states.pop()
            self.table_states.setVerticalHeaderLabels(self.header_vertical_states)

            self.table_states.removeRow(self.table_states.rowCount() - 1)
            self.delete_costate()
        else:
            print("not enough element to delete!")

    def new_costate(self):
        """This function adds 2 rows in the table (1 for state, 1 for co-state)"""
        count = self.table_costates.rowCount()
        self.table_costates.insertRow(count)
        self.header_vertical_costates += ["costate"]
        self.table_costates.setVerticalHeaderLabels(self.header_vertical_costates)

    def delete_costate(self):
        """This function removes 2 rows in the table (1 for state, 1 for co-state)"""
        if len(self.header_vertical_costates) > 1:
            self.header_vertical_costates.pop()
            self.table_costates.setVerticalHeaderLabels(self.header_vertical_costates)

            self.table_costates.removeRow(self.table_costates.rowCount() - 1)
        else:
            print("not enough element to delete!")

    def clear_all(self):
        self.table_states.setRowCount(0)
        self.table_costates.setRowCount(0)
        self.new_state()
