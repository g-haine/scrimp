from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import (
    QHBoxLayout,
    QPushButton,
    QLineEdit,
    QGridLayout,
    QTableWidget,
    QComboBox,
    QTableWidgetItem,
)
from PyQt5.QtCore import Qt
from utils.GUI import gui_pages, gui_width, gui_height, Help, check_black_listed_words


class Window(QtWidgets.QWidget):
    """This class defines the add FEM and coFEM page of the GUI and asks to insert
    the FEMs and co-FEM realted to the new distribuited FEM-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self, session):
        QtWidgets.QWidget.__init__(self)
        self.session = session

        self.setWindowTitle("Definition of FEM/s")
        self.setFixedWidth(gui_width)
        self.setFixedHeight(gui_height)

        self.layout = QGridLayout()

        # create a QTableWidget FEMs
        self.table_FEMs = QTableWidget()

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

        self.button_delete_FEM = QPushButton("Remove selected FEM")
        self.button_delete_FEM.clicked.connect(self.delete_FEM)

        self.button_clear_all = QPushButton("Clear All")
        self.button_clear_all.clicked.connect(self.clear_all)

        self.button_next = QPushButton("Next >")
        self.button_next.clicked.connect(self.next_page)

        self.button_prev = QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        self.layout.addWidget(self.table_FEMs, 1, 0, 1, 3)
        self.layout.addWidget(self.button_clear_all, 0, 1)
        self.layout.addWidget(self.button_add_FEM, 0, 2, Qt.AlignTop)
        self.layout.addWidget(self.button_delete_FEM, 0, 3, Qt.AlignTop)

        self.layout.addWidget(self.button_next, 4, 3)
        self.layout.addWidget(self.button_prev, 4, 2)

        # create navigation list
        self.comboBox = QComboBox()
        self.comboBox.addItems(gui_pages)
        self.comboBox.setCurrentText("add_fem_page")

        # There is an alternate signal to send the text.
        self.comboBox.currentTextChanged.connect(self.text_changed)
        self.layout.addWidget(self.comboBox, 4, 1)

        self.setLayout(self.layout)

        self.help = Help(self.layout, 3, 3)
        self.table_FEMs.cellClicked.connect(self.update_help)

    def update_help(self):
        """This function updates the Help object through its update_fields method.
        A text, a description and an example are prepared to be passed to the abovementioned method.
        """
        example = ""
        col = self.table_FEMs.currentColumn()

        if col is not None:
            text = self.table_FEMs.horizontalHeaderItem(col).text()
            print(f"col:{col},text:{text},selection:FEM")

            self.layout.itemAt(self.layout.count() - 1).widget().show()
            if col == 0:
                description = "Choose a name for the FEM."
                example = """This has to be:\n
                the name of a State\n
                the name of a Port\n
                the name of a Control Port"""

            elif col == 1:
                description = "Choose the order of your FEM."
                example = "it must be an int"
            elif col == 2:
                description = "Select the type of FEM."

            self.help.updateFields(text, description, example)

        else:
            self.help.clear()
            self.layout.itemAt(self.layout.count() - 1).widget().hide()

    def text_changed(self, page):  # s is a str
        """This function allows the navigation trhough the navigation list.
        After checking the presence of black listed words, the function hides the current page for showing the selected one.

        Args:
            page (str): the name of the page.
        """
        self.comboBox.setCurrentText("add_fem_page")
        if not check_black_listed_words(self, self.table_FEMs, "FEMs"):
            self.switch_window.emit(page)
            self.hide()

    def update_page(self):
        """This function manages the update of the current page."""
        pass

    def next_page(self):
        """This function emits the signal to navigate to the next page."""
        if not check_black_listed_words(self, self.table_FEMs, "FEMs"):
            self.switch_window.emit("set_hamiltonian_page")
            self.hide()

    def previous_page(self):
        """This funcion emits the signal to navigate to the prvious page."""
        if not check_black_listed_words(self, self.table_FEMs, "FEMs"):
            self.switch_window.emit("add_control_port_page")
            self.hide()

    def fem_choice_clicked(self):
        description = "Select the type of FEM."
        text = "FEM"
        example = """\n
        CG stands for Continuous Galerkin\n
        DG stands for Discontinuous Galerkin"""
        self.help.updateFields(text, description, example)

    def new_FEM(self):
        """This function adds 1 row in the table for FEM"""
        count = self.table_FEMs.rowCount()
        self.table_FEMs.insertRow(count)
        self.header_vertical_FEMs += ["FEM"]
        self.table_FEMs.setVerticalHeaderLabels(self.header_vertical_FEMs)

        # create table to add in cell of table
        fem_choice = QComboBox()
        fem_choice.addItems(["CG", "DG"])
        fem_choice.textHighlighted.connect(self.fem_choice_clicked)
        self.table_FEMs.setCellWidget(count, 2, fem_choice)

        for i in range(self.table_FEMs.columnCount()):
            if i != 2:
                # others
                new_value = QTableWidgetItem("")
                self.table_FEMs.setItem(count, i, new_value)

    def delete_FEM(self):
        """This function removes 1 row from the table"""
        if len(self.header_vertical_FEMs) > 1:
            self.header_vertical_FEMs.pop()
            self.table_FEMs.setVerticalHeaderLabels(self.header_vertical_FEMs)

            self.table_FEMs.removeRow(self.table_FEMs.currentRow())

        else:
            print("not enough element to delete!")

    def clear_all(self):
        """This function removes all the rows from the table."""
        self.table_FEMs.setRowCount(0)
        self.new_FEM()
