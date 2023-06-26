from PyQt5 import QtCore, QtWidgets, QtGui
from PyQt5.QtCore import Qt
from utils.GUI import gui_pages

from PyQt5.QtWidgets import (
    QListWidget,
    QListWidgetItem,
    QAbstractItemView,
    QPushButton,
    QLabel,
    QLineEdit,
    QGridLayout,
    QTableWidget,
    QComboBox,
    QTableWidgetItem,
)


class Window(QtWidgets.QWidget):
    """This class defines the set domain page of the GUI and asks to define domain of the new distribuited port-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self):
        QtWidgets.QWidget.__init__(self)

        self.setWindowTitle("Define the domain for the dpHs")
        self.setFixedWidth(1700)
        self.setFixedHeight(600)

        self.layout = QGridLayout()

        label_type_domain = QLabel(
            '<font size="4">'
            + "<p>Select a domain for your dpHs from the built-in<\p>"
            + "<p> list or define your own one:<\p>"
            + "</font>"
        )

        # creating a QListWidget
        self.list_widget = QListWidget(self)

        # list widget items
        item_separator_1D = QListWidgetItem("1D:")
        item_separator_2D = QListWidgetItem("2D:")
        item1 = QListWidgetItem("Rectangle")
        item2 = QListWidgetItem("Disk")
        item3 = QListWidgetItem("Annulus")
        item4 = QListWidgetItem("Segment")
        item5 = QListWidgetItem("Other")

        item_separator_1D.setFlags(Qt.NoItemFlags)
        item_separator_2D.setFlags(Qt.NoItemFlags)
        item_separator_1D.setForeground(QtGui.QBrush(QtGui.QColor(255, 255, 255)))
        item_separator_2D.setForeground(QtGui.QBrush(QtGui.QColor(255, 255, 255)))
        item_separator_1D.setBackground(QtGui.QColor(102, 178, 255))
        item_separator_2D.setBackground(QtGui.QColor(102, 178, 255))

        # adding items to the list widget
        self.list_widget.addItem(item_separator_1D)
        self.list_widget.addItem(item4)
        self.list_widget.addItem(item_separator_2D)
        self.list_widget.addItem(item1)
        self.list_widget.addItem(item2)
        self.list_widget.addItem(item3)
        self.list_widget.addItem(item5)

        # setting selection mode property
        self.list_widget.setSelectionMode(QAbstractItemView.SingleSelection)
        self.list_widget.itemClicked.connect(self.update_table)

        label_parameter = QLabel('<font size="4"> Parameters:</font>')

        self.button_next = QPushButton("Next >")
        self.button_next.clicked.connect(self.next_page)

        self.button_prev = QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        # create a QTableWidget
        self.table = QTableWidget()
        self.table.setRowCount(1)
        self.table.setColumnCount(2)

        # adding header to the table
        header_horizontal = ["Name", "Value"]

        self.table.setHorizontalHeaderLabels(header_horizontal)

        for i, _ in enumerate(header_horizontal):
            self.table.setColumnWidth(i, 150)

        self.button_add = QPushButton("Add")
        self.button_add.clicked.connect(self.new_rows)

        self.layout.addWidget(label_type_domain, 1, 0)
        self.layout.addWidget(self.list_widget, 2, 1)
        self.layout.addWidget(label_parameter, 3, 0)
        self.layout.addWidget(self.button_add, 4, 1, Qt.AlignTop)
        self.layout.addWidget(self.table, 4, 0)
        self.layout.addWidget(self.button_prev, 5, 3)
        self.layout.addWidget(self.button_next, 5, 4)

        # create navigation list
        self.comboBox = QComboBox()
        self.comboBox.addItems(gui_pages)
        self.comboBox.setCurrentText("set_domain_page")

        # There is an alternate signal to send the text.
        self.comboBox.currentTextChanged.connect(self.text_changed)
        self.layout.addWidget(self.comboBox, 5, 2)

        self.layout.itemAt(3).widget().hide()
        self.setLayout(self.layout)

    def text_changed(self, page):  # s is a str
        self.comboBox.setCurrentText("set_domain_page")
        self.switch_window.emit(page)
        self.hide()

    def update_table(self):
        selection = self.list_widget.currentItem().text()

        if selection == "Other":
            self.layout.itemAt(3).widget().show()
        else:
            self.layout.itemAt(3).widget().hide()

        if selection == "Rectangle":
            # remove all the rows
            self.table.setRowCount(0)
            # # add 3 rows
            for _ in range(3):
                self.table.insertRow(self.table.rowCount())

            self.table.setItem(0, 0, QTableWidgetItem("L"))
            self.table.setItem(1, 0, QTableWidgetItem("l"))
            self.table.setItem(2, 0, QTableWidgetItem("h"))

        elif selection == "Disk":
            # remove all the rows
            self.table.setRowCount(0)
            # # add 3 rows
            for _ in range(2):
                self.table.insertRow(self.table.rowCount())

            self.table.setItem(0, 0, QTableWidgetItem("R"))
            self.table.setItem(1, 0, QTableWidgetItem("h"))

        elif selection == "Annulus":
            # remove all the rows
            self.table.setRowCount(0)
            # # add 3 rows
            for _ in range(3):
                self.table.insertRow(self.table.rowCount())

            self.table.setItem(0, 0, QTableWidgetItem("R"))
            self.table.setItem(1, 0, QTableWidgetItem("r"))
            self.table.setItem(2, 0, QTableWidgetItem("h"))

        elif selection == "Segment":
            # remove all the rows
            self.table.setRowCount(0)
            # # add 3 rows
            for _ in range(2):
                self.table.insertRow(self.table.rowCount())

            self.table.setItem(0, 0, QTableWidgetItem("L"))
            self.table.setItem(1, 0, QTableWidgetItem("h"))

        else:
            # remove all the rows
            self.table.setRowCount(0)

    def new_rows(self):
        """This function adds 2 rows in the table (1 for state, 1 for co-state)"""
        count = self.table.rowCount()
        self.table.insertRow(count)

    def next_page(self):
        """This funciont emit the signal to navigate to the next page."""
        self.switch_window.emit("add_state_costate_page")
        self.hide()

    def previous_page(self):
        """This funciont emit the signal to navigate to the previous page."""
        self.switch_window.emit("create_dphs_page")
        self.hide()
