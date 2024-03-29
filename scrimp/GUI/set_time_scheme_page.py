from PyQt5 import QtCore, QtWidgets
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
    QTableWidgetItem,
    QComboBox,
    QTableWidgetItem,
    QCheckBox,
)


class Window(QtWidgets.QWidget):
    """This class defines the set domain page of the GUI and asks to define domain of the new distribuited port-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def check_state(self):
        """This function checks wether or not the checkbox has been checked and update the page accordingly."""
        print("clicked")
        if self.checkBox_answer.isChecked():
            self.set_self_layout()
        else:
            self.clear_self_layout()

    def clear_self_layout(self):
        """This function clears the layout object."""
        for i in range(2, 7):
            self.layout.itemAt(i).widget().hide()

    def set_self_layout(self):
        """This function presets the layout object."""
        for i in range(2, 7):
            self.layout.itemAt(i).widget().show()

    def __init__(self, session):
        QtWidgets.QWidget.__init__(self)
        self.session = session

        self.setWindowTitle("Define the Time scheme")
        self.setFixedWidth(gui_width)
        self.setFixedHeight(gui_height)

        self.layout = QGridLayout()

        self.label_question = QLabel("Check the box if you need a time scheme:")

        self.checkBox_answer = QCheckBox()
        self.checkBox_answer.toggled.connect(self.check_state)

        self.label_type_domain = QLabel(
            '<font size="4">'
            + "<p>Select the type for your time scheme:<\p>"
            + "</font>"
        )

        # creating a QListWidget
        self.list_widget = QListWidget(self)

        # list widget items
        item1 = QListWidgetItem("First")
        item2 = QListWidgetItem("Second")
        item3 = QListWidgetItem("Third")
        item4 = QListWidgetItem("Fourth")
        item5 = QListWidgetItem("Other")

        # adding items to the list widget
        self.list_widget.addItem(item1)
        self.list_widget.addItem(item2)
        self.list_widget.addItem(item3)
        self.list_widget.addItem(item4)
        self.list_widget.addItem(item5)

        # setting selection mode property
        self.list_widget.setSelectionMode(QAbstractItemView.SingleSelection)
        self.list_widget.itemClicked.connect(self.update_table)

        self.label_parameter = QLabel('<font size="4"> Parameters:</font>')

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

        self.layout.addWidget(self.label_question, 0, 0)
        self.layout.addWidget(self.checkBox_answer, 0, 1)
        self.layout.addWidget(self.label_type_domain, 1, 0)
        self.layout.addWidget(self.list_widget, 2, 1)
        self.layout.addWidget(self.label_parameter, 3, 0)
        self.layout.addWidget(self.button_add, 4, 1, Qt.AlignTop)
        self.layout.addWidget(self.table, 4, 0)
        self.layout.addWidget(self.button_prev, 5, 2)
        self.layout.addWidget(self.button_next, 5, 3)

        # create navigation list
        self.comboBox = QComboBox()
        self.comboBox.addItems(gui_pages)
        self.comboBox.setCurrentText("set_time_scheme_page")

        # There is an alternate signal to send the text.
        self.comboBox.currentTextChanged.connect(self.text_changed)
        self.layout.addWidget(self.comboBox, 5, 1)

        self.clear_self_layout()
        self.setLayout(self.layout)

        self.help = Help(self.layout, 4, 1)
        self.table.cellClicked.connect(self.update_help)

    def update_help(self):
        """This function updates the Help object through its update_fields method.
        A text, a description and an example are prepared to be passed to the abovementioned method.
        """

        item = self.table.currentItem()
        if item is not None:
            text = item.text()
            col = item.column()
            row = item.row()
            selection = self.list_widget.currentItem().text()
            print(f"col:{col},text:{text},selection:{selection}")
            if col == 0:
                if selection == "First":
                    description = None
                    example = None

                    if row == 0:
                        description = "This is ..."

                    elif row == 1:
                        description = "This is ..."

                    elif row == 2:
                        description = "This is ..."

                    elif row == 3:
                        description = "This is ..."

                    elif row == 4:
                        description = "This is ..."

                    elif row == 5:
                        description = "This is ..."

                    elif row == 6:
                        description = "This is ..."

                    elif row == 7:
                        description = "This is ..."

                    elif row == 8:
                        description = "This is ..."

                    elif row == 9:
                        description = "This is ..."

                    self.help.updateFields(text, description, example)

                elif selection == "Second":
                    description = None
                    example = None

                    if row == 0:
                        description = "This is ..."

                    elif row == 1:
                        description = "This is ..."

                    elif row == 2:
                        description = "This is ..."

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

        self.comboBox.setCurrentText("set_time_scheme_page")
        if not check_black_listed_words(self, self.table, "Set Time scheme tabke"):
            self.switch_window.emit(page)
            self.hide()

    def update_table(self):
        """This function updates the field of the selected time scheme."""
        selection = self.list_widget.currentItem().text()
        self.help.clear()

        if selection == "Other":
            self.layout.itemAt(5).widget().show()
        else:
            self.layout.itemAt(5).widget().hide()

        if selection == "First":
            # remove all the rows
            self.table.setRowCount(0)
            # # add 3 rows
            for _ in range(11):
                self.table.insertRow(self.table.rowCount())

            self.table.setItem(0, 0, QTableWidgetItem("ts_equation_type"))
            self.table.setItem(1, 0, QTableWidgetItem("ts_type"))
            self.table.setItem(2, 0, QTableWidgetItem("ksp_type"))
            self.table.setItem(3, 0, QTableWidgetItem("pc_type"))
            self.table.setItem(4, 0, QTableWidgetItem("pc_factor_mat_solver_type"))
            self.table.setItem(5, 0, QTableWidgetItem("t_0"))
            self.table.setItem(6, 0, QTableWidgetItem("t_f"))
            self.table.setItem(7, 0, QTableWidgetItem("dt"))
            self.table.setItem(8, 0, QTableWidgetItem("dt_save"))
            self.table.setItem(9, 0, QTableWidgetItem("ts_adapt_dt_max"))
            self.table.setItem(10, 0, QTableWidgetItem("init_step"))

        elif selection == "Second":
            # remove all the rows
            self.table.setRowCount(0)
            # # add 3 rows
            for _ in range(3):
                self.table.insertRow(self.table.rowCount())

            self.table.setItem(0, 0, QTableWidgetItem("t_f"))
            self.table.setItem(1, 0, QTableWidgetItem("dt"))
            self.table.setItem(2, 0, QTableWidgetItem("pc_type"))

        elif selection == "Third":
            # remove all the rows
            self.table.setRowCount(0)
            # # add 3 rows
            for _ in range(3):
                self.table.insertRow(self.table.rowCount())

            self.table.setItem(0, 0, QTableWidgetItem("R"))
            self.table.setItem(1, 0, QTableWidgetItem("r"))
            self.table.setItem(2, 0, QTableWidgetItem("h"))

        elif selection == "Fourth":
            # remove all the rows
            self.table.setRowCount(0)
            # # add 3 rows
            for _ in range(2):
                self.table.insertRow(self.table.rowCount())

            self.table.setItem(0, 0, QTableWidgetItem("L"))
            self.table.setItem(1, 0, QTableWidgetItem("h"))

        else:
            # remove all the rows
            self.table.setRowCount(1)

        rows = self.table.rowCount()
        for row in range(rows):
            new_value = QTableWidgetItem("")
            self.table.setItem(row, 1, new_value)

    def new_rows(self):
        """This function adds 1 row from the table"""
        count = self.table.rowCount()
        self.table.insertRow(count)

    def update_page(self):
        """This function manages the update of the current page."""
        pass

    def next_page(self):
        """This function emits the signal to navigate to the next page."""
        if not check_black_listed_words(self, self.table, "Set Time scheme tabke"):
            self.switch_window.emit("generate_code_page")
            self.hide()

    def previous_page(self):
        """This function emits the signal to navigate to the previous page."""
        if not check_black_listed_words(self, self.table, "Set Time scheme tabke"):
            self.switch_window.emit("add_initial_value_page")
            self.hide()
