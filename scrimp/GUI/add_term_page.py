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
    """This class defines the add term and coterm page of the GUI and asks to insert
    the terms and co-term realted to the new distribuited term-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self,session):
        QtWidgets.QWidget.__init__(self)
        self.session  = session

        self.setWindowTitle("Definition of Term/s")
        self.setFixedWidth(gui_width)
        self.setFixedHeight(gui_height)
        # self.setGeometry(100, 100, 600, 300)

        self.layout = QGridLayout()

        # self.line_edit = QLineEdit()
        # layout.addWidget(self.line_edit)

        # create a QTableWidget terms
        self.table_terms = QTableWidget()
        # self.table_terms.setRowCount(1)

        # adding header to the table
        header_horizontal_terms = ["Description",
                                   "Expression", "Regions", "Mesh ID"]
        self.table_terms.setColumnCount(len(header_horizontal_terms))

        self.header_vertical_terms = ["term"]
        self.table_terms.setHorizontalHeaderLabels(header_horizontal_terms)
        self.table_terms.setVerticalHeaderLabels(self.header_vertical_terms)

        # adjust size columns of horizontal header
        for i, _ in enumerate(header_horizontal_terms):
            self.table_terms.setColumnWidth(i, 150)

        self.button_add_term = QPushButton("Add term")
        self.button_add_term.clicked.connect(self.new_term)

        self.button_delete_term = QPushButton("Remove selected term")
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

        self.layout.addWidget(self.table_terms, 1, 0, 1, 3)
        # layout.addWidget(cell_double, 1, 3)
        self.layout.addWidget(self.button_clear_all, 0, 1)
        self.layout.addWidget(self.button_add_term, 0, 2, Qt.AlignTop)
        self.layout.addWidget(self.button_delete_term, 0, 3, Qt.AlignTop)

        self.layout.addWidget(self.button_next, 4, 3)
        self.layout.addWidget(self.button_prev, 4, 2)

        # create navigation list
        self.comboBox = QComboBox()
        self.comboBox.addItems(gui_pages)
        self.comboBox.setCurrentText("add_term_page")

        # There is an alternate signal to send the text.
        self.comboBox.currentTextChanged.connect(self.text_changed)
        self.layout.addWidget(self.comboBox, 4, 1)

        self.setLayout(self.layout)

        self.help = Help(self.layout, 3, 3)
        self.table_terms.cellClicked.connect(self.update_help)

        # self.new_term()

    def update_help(self):
        example = ""
        col = self.table_terms.currentColumn()

        if col is not None:
            text = self.table_terms.horizontalHeaderItem(col).text()
            print(f"col:{col},text:{text},selection:Term")

            self.layout.itemAt(self.layout.count() - 1).widget().show()
            if col == 0:
                description = "Choose a name or a description for the Term"
                example = "Kinetic energy"

            elif col == 1:
                description = "Choose the formula, using the Model variables, to define the term."
                example = " Parameters are allowed: 0.5*q.T.q"

            elif col == 2:
                description = "The region IDs of the mesh where the expression has to be evaluated. If more than one use a coma ',' to separate them with no spaces in between."

            elif col == 3:
                description = "Tthe mesh ID of the mesh where the regions belong to."
                example = "Default is 0.<br>If multiple: 0,1,2"

            self.help.updateFields(text, description, example)

        else:
            self.help.clear()
            self.layout.itemAt(self.layout.count() - 1).widget().hide()

    def text_changed(self, page):  # s is a str
        self.comboBox.setCurrentText("add_term_page")
        self.switch_window.emit(page)
        self.hide()

    def update_page(self):
        pass

    def next_page(self):
        """This funciont emit the signal to navigate to the next page."""
        self.switch_window.emit("add_brick_page")
        self.hide()

    def previous_page(self):
        """This funcion emits the signal to navigate to the prvious page."""
        self.switch_window.emit("set_hamiltonian_page")
        self.hide()

    def new_term(self):
        """This function adds 1 row in the table for term"""
        count = self.table_terms.rowCount()
        self.table_terms.insertRow(count)
        self.header_vertical_terms += ["term"]
        self.table_terms.setVerticalHeaderLabels(self.header_vertical_terms)
        # set defaults
        # mesh_id
        new_value = QTableWidgetItem("0")
        self.table_terms.setItem(count, 3, new_value)

        for i in range(self.table_terms.columnCount()):
            if i not in [3]:
                # others
                new_value = QTableWidgetItem("")
                self.table_terms.setItem(count, i, new_value)

    def delete_term(self):
        """This function removes 2 rows in the table (1 for term, 1 for co-term)"""
        if len(self.header_vertical_terms) > 1:
            self.header_vertical_terms.pop()
            self.table_terms.setVerticalHeaderLabels(
                self.header_vertical_terms)

            self.table_terms.removeRow(self.table_terms.currentRow())

        else:
            print("not enough element to delete!")

    def clear_all(self):
        self.table_terms.setRowCount(0)
        self.new_term()
