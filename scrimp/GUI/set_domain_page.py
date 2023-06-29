from PyQt5 import QtCore, QtWidgets, QtGui
from PyQt5.QtCore import Qt
from utils.GUI import gui_pages, gui_width, gui_height, Help

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
        self.setFixedWidth(gui_width)
        self.setFixedHeight(gui_height)

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
        self.table.setColumnCount(2)

        # adding header to the table
        header_horizontal = ["Name", "Value"]

        self.table.setHorizontalHeaderLabels(header_horizontal)

        for i, _ in enumerate(header_horizontal):
            self.table.setColumnWidth(i, 150)

        self.button_add = QPushButton("Add")
        self.button_add.clicked.connect(self.new_rows)
        self.button_add.resize(100, 50)

        self.button_delete = QPushButton("Remove")
        self.button_delete.clicked.connect(self.delete)

        # create navigation list
        self.comboBox = QComboBox()
        self.comboBox.addItems(gui_pages)
        self.comboBox.setCurrentText("set_domain_page")

        # There is an alternate signal to send the text.
        self.comboBox.currentTextChanged.connect(self.text_changed)

        self.layout.addWidget(label_type_domain, 1, 0)
        self.layout.addWidget(self.list_widget, 2, 0)
        self.layout.addWidget(label_parameter, 3, 0)
        self.layout.addWidget(self.button_add, 4, 1, Qt.AlignTop)
        self.layout.addWidget(self.button_delete, 4, 2, Qt.AlignTop)
        self.layout.addWidget(self.table, 4, 0)
        self.layout.addWidget(self.comboBox, 5, 3)
        self.layout.addWidget(self.button_prev, 5, 4)
        self.layout.addWidget(self.button_next, 5, 5)

        self.layout.itemAt(3).widget().hide()
        self.layout.itemAt(4).widget().hide()

        self.setLayout(self.layout)

        self.help = Help(self.layout, 4, 2)
        self.table.itemClicked.connect(self.update_help)

    def update_help(self):
        item = self.table.currentItem()
        if item is not None:
            text = item.text()
            col = item.column()
            selection = self.list_widget.currentItem().text()
            print(f"col:{col},text:{text},selection:{selection}")

            if selection != "Other":
                description = None
                example = None

                if col == 0:
                    self.layout.itemAt(self.layout.count() - 1).widget().show()

                    if text == "L":
                        if selection == "Segment":
                            description = "Lenght of the segment"
                            example = "in a segment AB, L=15, implies a distance of 15 meters between the vertex AB."

                        elif selection == "Rectangle":
                            description = "Lenght of the longer base of the rectangle"
                            example = "in a rectangle ABCD, if AB is the longer side, L is its lenght."

                    if text == "l":
                        description = "Lenght of the shorter base of the rectangle"
                        example = "in a rectangle ABCD, if BC is the shorter side, l is its lenght."

                    elif text == "h":
                        description = "The step of integration"
                        example = ""

                    elif text == "R":
                        if selection == "Disk":
                            description = "Rayon of the disk"
                            example = ""

                        elif selection == "Annulus":
                            description = "Lenght of the longer rayon of the annulus"
                            example = (
                                "in an annulus, R is the rayon of the outer circle."
                            )

                    elif text == "r":
                        description = "Lenght of the shorter rayon of the annulus"
                        example = "in an annulus, R is the rayon of the inner circle."

                    self.help.updateFields(text, description, example)

            else:
                self.help.clear()
                self.layout.itemAt(self.layout.count() - 1).widget().hide()

        else:
            self.help.clear()
            self.layout.itemAt(self.layout.count() - 1).widget().hide()

    def text_changed(self, page):  # s is a str
        self.comboBox.setCurrentText("set_domain_page")
        self.switch_window.emit(page)
        self.hide()

    def delete(self):
        """This function removes 2 rows in the table (1 for state, 1 for co-state)"""
        if self.table.rowCount() >= 1:
            self.table.removeRow(self.table.rowCount() - 1)
        else:
            print("not enough element to delete!")

    def update_table(self):
        selection = self.list_widget.currentItem().text()
        self.help.clear()
        # self.update_help()

        if selection == "Other":
            self.layout.itemAt(3).widget().show()
            self.layout.itemAt(4).widget().show()
        else:
            self.layout.itemAt(3).widget().hide()
            self.layout.itemAt(4).widget().hide()

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
