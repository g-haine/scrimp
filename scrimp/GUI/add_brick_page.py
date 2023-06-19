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
from utils.GUI import gui_pages


class Window(QtWidgets.QWidget):
    """This class defines the add brick and cobrick page of the GUI and asks to insert
    the bricks and co-brick realted to the new distribuited brick-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self):
        QtWidgets.QWidget.__init__(self)

        self.setWindowTitle("Definition of Control brick/s")
        self.setFixedWidth(1700)
        self.setFixedHeight(600)
        # self.setGeometry(100, 100, 600, 300)

        layout = QGridLayout()

        # self.line_edit = QLineEdit()
        # layout.addWidget(self.line_edit)

        # create a QTableWidget bricks
        self.table_bricks = QTableWidget()
        self.table_bricks.setRowCount(1)

        # adding header to the table
        header_horizontal_bricks = [
            "Name",
            "Form",
            "Regions",
            "Linear",
            "dt",
            "Position",
            "Mesh ID",
        ]
        self.table_bricks.setColumnCount(len(header_horizontal_bricks))

        self.header_vertical_bricks = ["brick"]
        self.table_bricks.setHorizontalHeaderLabels(header_horizontal_bricks)
        self.table_bricks.setVerticalHeaderLabels(self.header_vertical_bricks)

        # adjust size columns of horizontal header
        for i, _ in enumerate(header_horizontal_bricks):
            self.table_bricks.setColumnWidth(i, 150)

        self.button_add_brick = QPushButton("Add brick")
        self.button_add_brick.clicked.connect(self.new_brick)

        self.button_delete_brick = QPushButton("Remove brick")
        self.button_delete_brick.clicked.connect(self.delete_brick)

        self.button_clear_all = QPushButton("Clear All")
        self.button_clear_all.clicked.connect(self.clear_all)

        # layout_buttons_brick = QHBoxLayout()

        # layout_buttons_brick.addWidget(self.button_add_brick)
        # layout_buttons_brick.addWidget(self.button_delete_brick)

        # cell_double = QTableWidget(layout_buttons_brick)

        self.button_next = QPushButton("Next >")
        self.button_next.clicked.connect(self.next_page)

        self.button_prev = QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        layout.addWidget(self.table_bricks, 1, 0, 1, 3)
        # layout.addWidget(cell_double, 1, 3)
        layout.addWidget(self.button_clear_all, 0, 1)
        layout.addWidget(self.button_add_brick, 0, 2, Qt.AlignTop)
        layout.addWidget(self.button_delete_brick, 0, 3, Qt.AlignTop)

        layout.addWidget(self.button_next, 4, 3)
        layout.addWidget(self.button_prev, 4, 2)

        self.comboBox = QComboBox()
        self.comboBox.addItems(gui_pages)
        self.comboBox.setCurrentText("add_brick_page")

        # There is an alternate signal to send the text.
        self.comboBox.currentTextChanged.connect(self.text_changed)
        layout.addWidget(self.comboBox, 4, 1)

        self.setLayout(layout)

    def text_changed(self, page):  # s is a str
        self.comboBox.setCurrentText("add_brick_page")
        self.switch_window.emit(page)
        self.hide()

    def next_page(self):
        """This funciont emit the signal to navigate to the next page."""
        self.switch_window.emit("add_expression_page")
        self.hide()

    def previous_page(self):
        """This funcion emits the signal to navigate to the prvious page."""
        self.switch_window.emit("add_term_page")
        self.hide()

    def new_brick(self):
        """This function adds 2 rows in the table (1 for brick, 1 for co-brick)"""
        count = self.table_bricks.rowCount()
        self.table_bricks.insertRow(count)
        self.header_vertical_bricks += ["brick"]
        self.table_bricks.setVerticalHeaderLabels(self.header_vertical_bricks)

    def delete_brick(self):
        """This function removes 2 rows in the table (1 for brick, 1 for co-brick)"""
        if len(self.header_vertical_bricks) > 1:
            self.header_vertical_bricks.pop()
            self.table_bricks.setVerticalHeaderLabels(self.header_vertical_bricks)

            self.table_bricks.removeRow(self.table_bricks.rowCount() - 1)

        else:
            print("not enough element to delete!")

    def clear_all(self):
        self.table_bricks.setRowCount(0)
        self.new_brick()
