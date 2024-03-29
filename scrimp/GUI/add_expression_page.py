from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import (
    QPushButton,
    QGridLayout,
    QTableWidget,
    QComboBox,
    QTableWidgetItem,
)
from PyQt5.QtCore import Qt
from utils.GUI import gui_pages, gui_width, gui_height, Help, check_black_listed_words


class Window(QtWidgets.QWidget):
    """This class defines the add expression and coexpression page of the GUI and asks to insert
    the expressions and co-expression realted to the new distribuited expression-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self, session):
        QtWidgets.QWidget.__init__(self)
        self.session = session

        self.setWindowTitle("Definition of Expression/s")
        self.setFixedWidth(gui_width)
        self.setFixedHeight(gui_height)

        self.layout = QGridLayout()

        # create a QTableWidget expressions
        self.table_expressions = QTableWidget()

        # adding header to the table
        header_horizontal_expressions = ["Control Port", "Expression"]
        self.table_expressions.setColumnCount(len(header_horizontal_expressions))

        self.header_vertical_expressions = ["expression"]
        self.table_expressions.setHorizontalHeaderLabels(header_horizontal_expressions)
        self.table_expressions.setVerticalHeaderLabels(self.header_vertical_expressions)

        # adjust size columns of horizontal header
        for i, _ in enumerate(header_horizontal_expressions):
            self.table_expressions.setColumnWidth(i, 150)

        self.button_add_expression = QPushButton("Add expression")
        self.button_add_expression.clicked.connect(self.new_expression)

        self.button_delete_expression = QPushButton("Remove selected expression")
        self.button_delete_expression.clicked.connect(self.delete_expression)

        self.button_clear_all = QPushButton("Clear All")
        self.button_clear_all.clicked.connect(self.clear_all)

        self.button_next = QPushButton("Next >")
        self.button_next.clicked.connect(self.next_page)

        self.button_prev = QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        self.layout.addWidget(self.table_expressions, 1, 0, 1, 3)
        self.layout.addWidget(self.button_clear_all, 0, 1)
        self.layout.addWidget(self.button_add_expression, 0, 2, Qt.AlignTop)
        self.layout.addWidget(self.button_delete_expression, 0, 3, Qt.AlignTop)
        self.layout.addWidget(self.button_prev, 4, 2)
        self.layout.addWidget(self.button_next, 4, 3)

        # create navigation list
        self.comboBox = QComboBox()
        self.comboBox.addItems(gui_pages)
        self.comboBox.setCurrentText("add_expression_page")

        # There is an alternate signal to send the text.
        self.comboBox.currentTextChanged.connect(self.text_changed)
        self.layout.addWidget(self.comboBox, 4, 1)

        self.setLayout(self.layout)

        self.help = Help(self.layout, 3, 3)
        self.table_expressions.cellClicked.connect(self.update_help)

    def update_help(self):
        """This function updates the Help object through its update_fields method.
        A text, a description and an example are prepared to be passed to the abovementioned method.
        """
        example = ""
        col = self.table_expressions.currentColumn()

        if col is not None:
            text = self.table_expressions.horizontalHeaderItem(col).text()
            print(f"col:{col},text:{text},selection:Expression")

            self.layout.itemAt(self.layout.count() - 1).widget().show()

            if col == 0:
                description = "Insert the name of the Control Port at which the expression is tied to"

            elif col == 1:
                description = "Choose a source term expression for the Control Port"
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
        self.comboBox.setCurrentText("add_expression_page")
        if not check_black_listed_words(self, self.table_expressions, "Expression"):
            self.switch_window.emit(page)
            self.hide()

    def update_page(self):
        """This function manages the update of the current page."""
        pass

    def next_page(self):
        """This function emits the signal to navigate to the next page."""
        if not check_black_listed_words(self, self.table_expressions, "Expression"):
            self.switch_window.emit("add_initial_value_page")
            self.hide()

    def previous_page(self):
        """This funcion emits the signal to navigate to the prvious page."""
        if not check_black_listed_words(self, self.table_expressions, "Expression"):
            self.switch_window.emit("add_brick_page")
            self.hide()

    def new_expression(self):
        """This function adds 1 row in the table for expression"""
        count = self.table_expressions.rowCount()
        self.table_expressions.insertRow(count)
        self.header_vertical_expressions += ["expression"]
        self.table_expressions.setVerticalHeaderLabels(self.header_vertical_expressions)

        for i in range(self.table_expressions.columnCount()):
            new_value = QTableWidgetItem("")
            self.table_expressions.setItem(count, i, new_value)

    def delete_expression(self):
        """This function removes 1 row from the table"""
        if len(self.header_vertical_expressions) > 1:
            self.header_vertical_expressions.pop()
            self.table_expressions.setVerticalHeaderLabels(
                self.header_vertical_expressions
            )

            self.table_expressions.removeRow(self.table_expressions.currentRow())

        else:
            print("not enough element to delete!")

    def clear_all(self):
        """This function removes all the rows from the table."""
        self.table_expressions.setRowCount(0)
        self.new_expression()
