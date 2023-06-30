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
from utils.GUI import gui_pages, gui_width, gui_height, Help


class Window(QtWidgets.QWidget):
    """This class defines the add expression and coexpression page of the GUI and asks to insert
    the expressions and co-expression realted to the new distribuited expression-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self):
        QtWidgets.QWidget.__init__(self)

        self.setWindowTitle("Definition of Expression/s")
        self.setFixedWidth(gui_width)
        self.setFixedHeight(gui_height)
        # self.setGeometry(100, 100, 600, 300)

        self.layout = QGridLayout()

        # self.line_edit = QLineEdit()
        # layout.addWidget(self.line_edit)

        # create a QTableWidget expressions
        self.table_expressions = QTableWidget()
        self.table_expressions.setRowCount(1)

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

        self.button_delete_expression = QPushButton("Remove expression")
        self.button_delete_expression.clicked.connect(self.delete_expression)

        self.button_clear_all = QPushButton("Clear All")
        self.button_clear_all.clicked.connect(self.clear_all)

        # layout_buttons_expression = QHBoxLayout()

        # layout_buttons_expression.addWidget(self.button_add_expression)
        # layout_buttons_expression.addWidget(self.button_delete_expression)

        # cell_double = QTableWidget(layout_buttons_expression)

        self.button_next = QPushButton("Next >")
        self.button_next.clicked.connect(self.next_page)

        self.button_prev = QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        self.layout.addWidget(self.table_expressions, 1, 0, 1, 3)
        # layout.addWidget(cell_double, 1, 3)
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
        example = ""
        col = self.table_expressions.currentColumn()

        if col is not None:
            text = self.table_expressions.horizontalHeaderItem(col).text()
            print(f"col:{col},text:{text},selection:Expression")

            self.layout.itemAt(self.layout.count() - 1).widget().show()

            if col == 0:
                description = "Insert the name of the Control Port at wich the expression has to be bind"

            elif col == 1:
                description = "Choose a source term expression for the control port"
            self.help.updateFields(text, description, example)

        else:
            self.help.clear()
            self.layout.itemAt(self.layout.count() - 1).widget().hide()

    def text_changed(self, page):  # s is a str
        self.comboBox.setCurrentText("add_expression_page")
        self.switch_window.emit(page)
        self.hide()

    def next_page(self):
        """This funciont emit the signal to navigate to the next page."""
        self.switch_window.emit("add_initial_value_page")
        self.hide()

    def previous_page(self):
        """This funcion emits the signal to navigate to the prvious page."""
        self.switch_window.emit("add_brick_page")
        self.hide()

    def new_expression(self):
        """This function adds 2 rows in the table (1 for expression, 1 for co-expression)"""
        count = self.table_expressions.rowCount()
        self.table_expressions.insertRow(count)
        self.header_vertical_expressions += ["expression"]
        self.table_expressions.setVerticalHeaderLabels(self.header_vertical_expressions)

    def delete_expression(self):
        """This function removes 2 rows in the table (1 for expression, 1 for co-expression)"""
        if len(self.header_vertical_expressions) > 1:
            self.header_vertical_expressions.pop()
            self.table_expressions.setVerticalHeaderLabels(
                self.header_vertical_expressions
            )

            self.table_expressions.removeRow(self.table_expressions.rowCount() - 1)

        else:
            print("not enough element to delete!")

    def clear_all(self):
        self.table_expressions.setRowCount(0)
        self.new_expression()