from PyQt5 import QtCore, QtWidgets, QtGui
from PyQt5.QtCore import Qt
from utils.GUI import gui_pages, gui_width, gui_height, Help, check_black_listed_words

from PyQt5.QtWidgets import (
    QListWidget,
    QListWidgetItem,
    QAbstractItemView,
    QPushButton,
    QLabel,
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

    def __init__(self, session):
        QtWidgets.QWidget.__init__(self)
        self.session = session

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

        self.button_delete = QPushButton("Remove selected")
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
        """This function updates the Help object through its update_fields method.
        A text, a description and an example are prepared to be passed to the abovementioned method.
        """
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
                        description = "The mesh size parameter"
                        example = ""

                    elif text == "R":
                        if selection == "Disk":
                            description = "radius of the disk"
                            example = ""

                        elif selection == "Annulus":
                            description = "Lenght of the longer radius of the annulus"
                            example = (
                                "in an annulus, R is the radius of the outer circle."
                            )

                    elif text == "r":
                        description = "Lenght of the shorter radius of the annulus"
                        example = "in an annulus, R is the radius of the inner circle."

                    self.help.updateFields(text, description, example)

            else:
                self.help.clear()
                self.layout.itemAt(self.layout.count() - 1).widget().hide()

        else:
            self.help.clear()
            self.layout.itemAt(self.layout.count() - 1).widget().hide()

    def text_changed(self, page):  # s is a str
        """This function allows the navigation trhough the navigation list.
        After checking the presence of black listed words, the function hides the current page for showing the selected one.

        Args:
            page (str): the name of the page.
        """
        self.comboBox.setCurrentText("set_domain_page")
        if not check_black_listed_words(self, self.table, "Domain"):
            self.switch_window.emit(page)
            self.hide()

    def delete(self):
        """This function removes 1 row from the table"""
        if self.table.rowCount() >= 1:
            self.table.removeRow(self.table.currentRow())
        else:
            print("not enough element to delete!")

    def update_table(self):
        """This function updates the fields of the domain table."""
        selection = self.list_widget.currentItem().text()
        self.session["domain"] = selection

        self.help.clear()

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
        """This function adds 1 row from the table"""
        count = self.table.rowCount()
        self.table.insertRow(count)

    def update_page(self):
        """This function manages the update of the current page."""
        if (
            "read_from_file" in self.session.keys()
            and not self.session["set_domain_page"]["loaded_from_file"]
        ):
            self.load_session_from_file()

    def load_session_from_file(self):
        # set domain page
        domain_to_row = {
            "Segment": 1,
            "Rectangle": 3,
            "Disk": 4,
            "Annulus": 5,
            "Other": 6,
        }

        row = domain_to_row[
            self.session["read_from_file"]["dict"]["set_domain_page"]["list_widget"]
        ]
        self.list_widget.setCurrentRow(row)

        self.update_table()
        for i, param in enumerate(
            self.session["read_from_file"]["dict"]["set_domain_page"]["parameters"]
        ):
            self.table.setItem(i, 1, QTableWidgetItem(param))

        self.session["set_domain_page"]["loaded_from_file"] = True

    def next_page(self):
        """This function emits the signal to navigate to the next page."""
        if not check_black_listed_words(self, self.table, "Domain"):
            self.switch_window.emit("add_state_costate_page")
            self.hide()

    def previous_page(self):
        """This function emits the signal to navigate to the previous page."""
        if not check_black_listed_words(self, self.table, "Domain"):
            self.switch_window.emit("create_dphs_page")
            self.hide()
