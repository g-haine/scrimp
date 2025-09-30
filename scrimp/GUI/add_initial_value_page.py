from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import (
    QPushButton,
    QGridLayout,
    QTableWidget,
    QTableWidgetItem,
    QComboBox,
)
from PyQt5.QtCore import Qt
from utils.GUI import gui_pages, gui_width, gui_height, Help, check_black_listed_words


class Window(QtWidgets.QWidget):
    """This class defines the add initial_value and coinitial_value page of the GUI and asks to insert
    the initial_values and co-initial_value realted to the new distribuited initial_value-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self, session):
        QtWidgets.QWidget.__init__(self)
        self.session = session

        self.setWindowTitle("Definition of Initial value/s")
        self.setFixedWidth(gui_width)
        self.setFixedHeight(gui_height)

        self.layout = QGridLayout()

        # create a QTableWidget initial_values
        self.table_initial_values = QTableWidget()

        # adding header to the table
        header_horizontal_initial_values = ["Variable Name", "Initial Value"]
        self.table_initial_values.setColumnCount(len(header_horizontal_initial_values))

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

        self.button_delete_initial_value = QPushButton("Remove selected initial_value")
        self.button_delete_initial_value.clicked.connect(self.delete_initial_value)

        self.button_clear_all = QPushButton("Clear All")
        self.button_clear_all.clicked.connect(self.clear_all)

        self.button_next = QPushButton("Next >")
        self.button_next.clicked.connect(self.next_page)

        self.button_prev = QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        self.layout.addWidget(self.table_initial_values, 1, 0, 1, 3)
        self.layout.addWidget(self.button_clear_all, 0, 1)
        self.layout.addWidget(self.button_add_initial_value, 0, 2, Qt.AlignTop)
        self.layout.addWidget(self.button_delete_initial_value, 0, 3, Qt.AlignTop)
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
        """This function updates the Help object through its update_fields method.
        A text, a description and an example are prepared to be passed to the abovementioned method.
        """
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
        """This function allows the navigation trhough the navigation list.
        After checking the presence of black listed words, the function hides the current page for showing the selected one.

        Args:
            page (str): the name of the page.
        """
        self.comboBox.setCurrentText("add_initial_value_page")
        if not check_black_listed_words(
            self, self.table_initial_values, "Initial Values"
        ):
            self.switch_window.emit(page)
            self.hide()

    def update_page(self):
        """This function manages the update of the current page."""
        if (
            "read_from_file" in self.session.keys()
            and not self.session["add_initial_value_page"]["loaded_from_file"]
        ):
            self.load_session_from_file()

    def load_session_from_file(self):
        if (
            "initial_values"
            in self.session["read_from_file"]["dict"]["add_initial_value_page"].keys()
        ):

            self.table_initial_values.setRowCount(0)
            row = 0
            for initial_value in self.session["read_from_file"]["dict"][
                "add_initial_value_page"
            ]["initial_values"]:
                self.new_initial_value()
                for col, param in enumerate(initial_value):
                    self.table_initial_values.setItem(row, col, QTableWidgetItem(param))

                row += 1
        self.session["add_initial_value_page"]["loaded_from_file"] = True

    def next_page(self):
        """This function emits the signal to navigate to the next page."""
        if not check_black_listed_words(
            self, self.table_initial_values, "Initial Values"
        ):
            self.switch_window.emit("set_time_scheme_page")
            self.hide()

    def previous_page(self):
        """This funcion emits the signal to navigate to the prvious page."""
        if not check_black_listed_words(
            self, self.table_initial_values, "Initial Values"
        ):
            self.switch_window.emit("add_expression_page")
            self.hide()

    def new_initial_value(self):
        """This function adds 1 row in the table for initial_value"""
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
        """This function removes 1 row from the table"""
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
        """This function removes all the rows from the table."""
        self.table_initial_values.setRowCount(0)
        self.new_initial_value()
