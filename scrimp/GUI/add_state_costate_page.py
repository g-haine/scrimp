from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import (
    QHBoxLayout,
    QPushButton,
    QLineEdit,
    QGridLayout,
    QTableWidget,
    QTableWidgetItem,
    QComboBox,
)
from PyQt5.QtCore import Qt
from utils.GUI import gui_pages, gui_width, gui_height, Help


class Window(QtWidgets.QWidget):
    """This class defines the add state and costate page of the GUI and asks to insert
    the states and co-state realted to the new distribuited port-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self,session):
        QtWidgets.QWidget.__init__(self)
        self.session  = session

        self.setWindowTitle("Definition of State/s and Costate/s")
        self.setFixedWidth(gui_width)
        self.setFixedHeight(gui_height)
        # self.setGeometry(100, 100, 600, 300)

        self.layout = QGridLayout()

        # self.line_edit = QLineEdit()
        # self.layout.addWidget(self.line_edit)

        # create a QTableWidget States
        self.table_states = QTableWidget()
        # self.table_states.setRowCount(1)
        self.table_states.setColumnCount(5)
        # self.table_states.setGeometry(50, 100, 300, 300)

        # adding header to the table
        header_horizontal_states = [
            "Name", "Description", "Kind", "Region", "Mesh ID"]
        self.header_vertical_states = ["state"]
        self.table_states.setHorizontalHeaderLabels(header_horizontal_states)
        self.table_states.setVerticalHeaderLabels(self.header_vertical_states)

        # adjust size columns of horizontal header
        for i, _ in enumerate(header_horizontal_states):
            self.table_states.setColumnWidth(i, 150)

        self.button_add_state = QPushButton("Add State/Costate")
        self.button_add_state.clicked.connect(self.new_state)

        self.button_delete_state = QPushButton("Remove selected State/Costate")
        self.button_delete_state.clicked.connect(self.delete_state)

        self.button_clear_all = QPushButton("Clear All")
        self.button_clear_all.clicked.connect(self.clear_all)

        # layout_buttons_state = QHBoxLayout()

        # layout_buttons_state.addWidget(self.button_add_state)
        # layout_buttons_state.addWidget(self.button_delete_state)

        # cell_double = QTableWidget(layout_buttons_state)

        # create a QTableWidget Co-States
        self.table_costates = QTableWidget()
        # self.table_costates.setRowCount(1)
        self.table_costates.setColumnCount(4)
        # self.table_costates.setGeometry(50, 100, 300, 300)

        # adding header to the table
        header_horizontal_costates = [
            "Name",
            "Description",
            "State",
            "Substituted",
        ]
        self.header_vertical_costates = ["costate"]
        self.table_costates.setHorizontalHeaderLabels(
            header_horizontal_costates)
        self.table_costates.setVerticalHeaderLabels(
            self.header_vertical_costates)

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

        self.layout.addWidget(self.table_states, 1, 0, 1, 3)
        # self.layout.addWidget(cell_double, 1, 3)
        self.layout.addWidget(self.button_clear_all, 0, 1)
        self.layout.addWidget(self.button_add_state, 0, 2, Qt.AlignTop)
        self.layout.addWidget(self.button_delete_state, 0, 3, Qt.AlignTop)
        self.layout.addWidget(self.table_costates, 3, 0, 1, 3)
        # self.layout.addWidget(self.button_add_costate, 2, 2)
        # self.layout.addWidget(self.button_delete_costate, 2, 3)
        self.layout.addWidget(self.button_next, 4, 4)
        self.layout.addWidget(self.button_prev, 4, 3)

        # create navigation list
        self.comboBox = QComboBox()
        self.comboBox.addItems(gui_pages)
        self.comboBox.setCurrentText("add_state_costate_page")

        # There is an alternate signal to send the text.
        self.comboBox.currentTextChanged.connect(self.text_changed)
        self.layout.addWidget(self.comboBox, 4, 2)

        self.setLayout(self.layout)

        self.help = Help(self.layout, 3, 3)
        self.table_costates.cellClicked.connect(self.update_help_costate)
        self.table_states.cellClicked.connect(self.update_help_state)
        self.table_costates.cellClicked.connect(
            self.update_costate_table_by_state)

        self.new_state()

    def update_costate_table_by_state(self):
        row = self.table_costates.currentRow()
        self.table_costates.setItem(
            row, 2, self.table_states.item(row, 0).clone())

    def update_help_state(self):
        # item = self.table_states.currentItem()
        example = ""
        col = self.table_states.currentColumn()

        if col is not None:
            text = self.table_states.horizontalHeaderItem(col).text()
            print(f"col:{col},text:{text},selection:state")

            self.layout.itemAt(self.layout.count() - 1).widget().show()
            if col == 0:
                description = "Choose a name for the State"

            elif col == 1:
                description = "Choose a set of words that can describe your State"
                example = (
                    f"If your State is named 'T' you may describe it as 'Temperature'."
                )

            elif col == 2:
                description = """Choose what is the kind of your state.\nIt could be one of the following list:
                \n- scalar-field
                \n- vector-field
                \n- tensor-field"""
                example = "In 1D everything must be scalar-field."
            elif col == 3:
                description = "Choose which is the region that interest your state."

            elif col == 4:
                description = (
                    "Choose which is the ID for the mesh that interest your state."
                )
                example = "Default is 0."

            self.help.updateFields(text, description, example)

        else:
            self.help.clear()
            self.layout.itemAt(self.layout.count() - 1).widget().hide()

    def update_help_costate(self):
        # item = self.table_states.currentItem()
        example = ""
        col = self.table_costates.currentColumn()

        if col is not None:
            text = self.table_costates.horizontalHeaderItem(col).text()
            print(f"col:{col},text:{text},selection:costate")

            self.layout.itemAt(self.layout.count() - 1).widget().show()
            if col == 0:
                description = "Choose a name for the Costate"

            elif col == 1:
                description = "Choose a set of words that can describe your Costate"
                example = (
                    f"If your State is named 'T' you may describe it as 'Temperature'."
                )

            elif col == 2:
                description = "This is the State at which the Costate will be bounded."

            elif col == 3:
                description = "It is a boolean that defines whether to substitute the variable. Defaults to False"

            self.help.updateFields(text, description, example)

        else:
            self.help.clear()
            self.layout.itemAt(self.layout.count() - 1).widget().hide()

    def text_changed(self, page):  # s is a str
        self.comboBox.setCurrentText("add_state_costate_page")
        self.switch_window.emit(page)
        self.hide()

    def update_page(self):
        for row in range(self.table_states.rowCount()):
            comboBox = self.table_states.cellWidget(row,2)
            comboBox.clear()
            if "domain" in self.session.keys() and self.session["domain"] == "Segment":
                comboBox.addItems(
                    ["scalar-field"]
                )
            else:
                comboBox.addItems(
                    ["scalar-field", "vector-field", "tensor-field"]
                )

    def next_page(self):
        """This funciont emit the signal to navigate to the next page."""
        self.switch_window.emit("add_port_page")
        self.hide()

    def previous_page(self):
        """This funcion emits the signal to navigate to the prvious page."""
        self.switch_window.emit("set_domain_page")
        self.hide()

    def choice_clicked(self, text):
        def foo():
            print(text)
            description = ""
            example = ""

            if text == "Kind":
                description = """Choose what is the kind of your state.\nIt could be one of the following list:
                \n- scalar-field
                \n- vector-field
                \n- tensor-field"""
                example = "In 1D everything must be scalar-field."

            elif text == "Substituted":
                description = "It is a boolean that defines whether to substitute the variable. Defaults to False"
                example = ""
            self.help.updateFields(text, description, example)

        return foo

    def new_state(self):
        """This function adds 1 row per each table (1 for state, 1 for co-state)"""
        count = self.table_states.rowCount()
        self.table_states.insertRow(count)
        self.header_vertical_states += ["state"]
        self.table_states.setVerticalHeaderLabels(self.header_vertical_states)
        state_choice_kind = QComboBox()
        if "domain" in self.session.keys() and self.session["domain"] == "Segment":
            state_choice_kind.addItems(
                ["scalar-field"]
            )
        else:
            state_choice_kind.addItems(
                ["scalar-field", "vector-field", "tensor-field"]
            )

        state_choice_kind.textHighlighted.connect(self.choice_clicked("Kind"))
        self.table_states.setCellWidget(count, 2, state_choice_kind)
        # set defaults
        new_value = QTableWidgetItem("0")
        self.table_states.setItem(count, 4, new_value)
        self.new_costate()

        for i in range(self.table_states.columnCount()):
            if i not in [2, 4]:
                # others
                new_value = QTableWidgetItem("")
                self.table_states.setItem(count, i, new_value)

    def delete_state(self):
        """This function removes 2 rows in the table (1 for state, 1 for co-state)"""
        if len(self.header_vertical_states) > 1:
            self.header_vertical_states.pop()
            self.table_states.setVerticalHeaderLabels(
                self.header_vertical_states)

            self.table_states.removeRow(self.table_states.currentRow())
            self.delete_costate()
        else:
            print("not enough element to delete!")

    def new_costate(self):
        """This function adds 1 row in the table (1 for state, 1 for co-state)"""
        count = self.table_costates.rowCount()
        self.table_costates.insertRow(count)
        self.header_vertical_costates += ["costate"]
        self.table_costates.setVerticalHeaderLabels(
            self.header_vertical_costates)
        costate_choice_kind = QComboBox()
        costate_choice_kind.addItems(["False", "True"])

        costate_choice_kind.textHighlighted.connect(
            self.choice_clicked("Substituted"))
        self.table_costates.setCellWidget(count, 3, costate_choice_kind)

        for i in range(self.table_costates.columnCount()):
            if i not in [2, 3]:
                # others
                new_value = QTableWidgetItem("")
                self.table_costates.setItem(count, i, new_value)

    def delete_costate(self):
        """This function removes 2 rows in the table (1 for state, 1 for co-state)"""
        if len(self.header_vertical_costates) > 1:
            self.header_vertical_costates.pop()
            self.table_costates.setVerticalHeaderLabels(
                self.header_vertical_costates)

            self.table_costates.removeRow(self.table_costates.currentRow())
        else:
            print("not enough element to delete!")

    def clear_all(self):
        self.table_states.setRowCount(0)
        self.table_costates.setRowCount(0)
        self.new_state()
