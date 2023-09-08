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
    """This class defines the add initial_value and coinitial_value page of the GUI and asks to insert
    the initial_values and co-initial_value realted to the new distribuited initial_value-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self):
        QtWidgets.QWidget.__init__(self)

        self.setWindowTitle("Definition of Initial value/s")
        self.setFixedWidth(gui_width)
        self.setFixedHeight(gui_height)
        # self.setGeometry(100, 100, 600, 300)

        self.layout = QGridLayout()

        # self.line_edit = QLineEdit()
        # layout.addWidget(self.line_edit)

        # create a QTableWidget initial_values
        self.table_initial_values = QTableWidget()
        # self.table_initial_values.setRowCount(1)

        # adding header to the table
        header_horizontal_initial_values = ["Variable Name", "Initial Value"]
        self.table_initial_values.setColumnCount(
            len(header_horizontal_initial_values))

        self.header_vertical_initial_values = ["initial value"]
        self.table_initial_values.setHorizontalHeaderLabels(
            header_horizontal_initial_values
        )
        self.table_initial_values.setVerticalHeaderLabels(
            self.header_vertical_initial_values
        )

        # adjust size columns of horizontal header
        for i, _ in enumerate(header_horizontal_initial_values):
            self.table_initial_values.setColumnWidth(i, 150)

        self.button_add_initial_value = QPushButton("Add initial_value")
        self.button_add_initial_value.clicked.connect(self.new_initial_value)

        self.button_delete_initial_value = QPushButton("Remove initial_value")
        self.button_delete_initial_value.clicked.connect(
            self.delete_initial_value)

        self.button_clear_all = QPushButton("Clear All")
        self.button_clear_all.clicked.connect(self.clear_all)

        # layout_buttons_initial_value = QHBoxLayout()

        # layout_buttons_initial_value.addWidget(self.button_add_initial_value)
        # layout_buttons_initial_value.addWidget(self.button_delete_initial_value)

        # cell_double = QTableWidget(layout_buttons_initial_value)

        self.button_next = QPushButton("Next >")
        self.button_next.clicked.connect(self.next_page)

        self.button_prev = QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        self.layout.addWidget(self.table_initial_values, 1, 0, 1, 3)
        # layout.addWidget(cell_double, 1, 3)
        self.layout.addWidget(self.button_clear_all, 0, 1)
        self.layout.addWidget(self.button_add_initial_value, 0, 2, Qt.AlignTop)
        self.layout.addWidget(
            self.button_delete_initial_value, 0, 3, Qt.AlignTop)
        self.layout.addWidget(self.button_prev, 4, 2)
        self.layout.addWidget(self.button_next, 4, 3)

        # create navigation list
        self.comboBox = QComboBox()
        self.comboBox.addItems(gui_pages)
        self.comboBox.setCurrentText("add_initial_value_page")

        # There is an alternate signal to send the text.
        self.comboBox.currentTextChanged.connect(self.text_changed)
        self.layout.addWidget(self.comboBox, 4, 1)

        self.setLayout(self.layout)

        self.help = Help(self.layout, 3, 3)
        self.table_initial_values.cellClicked.connect(self.update_help)

    def update_help(self):
        example = ""
        col = self.table_initial_values.currentColumn()

        if col is not None:
            text = self.table_initial_values.horizontalHeaderItem(col).text()
            print(f"col:{col},text:{text},selection:Initial Value")

            self.layout.itemAt(self.layout.count() - 1).widget().show()
            if col == 0:
                description = "Choose the name of the variable for which the initial value has to be set"
                example = "It could be the name of a State."

            elif col == 1:
                description = "Choose the expression of the function to use for the indicated variable"
                example = "np.exp(-50*((x-1)*(x-1)+(y-0.5)*(y-0.5))**2)"

            self.help.updateFields(text, description, example)

        else:
            self.help.clear()
            self.layout.itemAt(self.layout.count() - 1).widget().hide()

    def text_changed(self, page):  # s is a str
        self.comboBox.setCurrentText("add_initial_value_page")
        self.switch_window.emit(page)
        self.hide()

    def next_page(self):
        """This funciont emit the signal to navigate to the next page."""
        self.switch_window.emit("set_time_scheme_page")
        self.hide()

    def previous_page(self):
        """This funcion emits the signal to navigate to the prvious page."""
        self.switch_window.emit("add_expression_page")
        self.hide()

    def new_initial_value(self):
        """This function adds 2 rows in the table (1 for initial_value, 1 for co-initial_value)"""
        count = self.table_initial_values.rowCount()
        self.table_initial_values.insertRow(count)
        self.header_vertical_initial_values += ["initial_value"]
        self.table_initial_values.setVerticalHeaderLabels(
            self.header_vertical_initial_values
        )

        for i in range(self.table_initial_values.columnCount()):
            new_value = QTableWidgetItem("")
            self.table_initial_values.setItem(count, i, new_value)

    def delete_initial_value(self):
        """This function removes 2 rows in the table (1 for initial_value, 1 for co-initial_value)"""
        if len(self.header_vertical_initial_values) > 1:
            self.header_vertical_initial_values.pop()
            self.table_initial_values.setVerticalHeaderLabels(
                self.header_vertical_initial_values
            )

            self.table_initial_values.removeRow(
                self.table_initial_values.rowCount() - 1
            )

        else:
            print("not enough element to delete!")

    def clear_all(self):
        self.table_initial_values.setRowCount(0)
        self.new_initial_value()
