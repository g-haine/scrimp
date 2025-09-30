from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import (
    QPushButton,
    QGridLayout,
    QTableWidget,
    QComboBox,
    QLabel,
    QTableWidgetItem,
)
from PyQt5.QtCore import Qt
from utils.GUI import gui_pages, gui_width, gui_height, Help, check_black_listed_words


class Window(QtWidgets.QWidget):
    """This class defines the add parameter and coparameter page of the GUI and asks to insert
    the parameters and co-parameter realted to the new distribuited parameter-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self, session):
        QtWidgets.QWidget.__init__(self)
        self.session = session

        self.setWindowTitle("Definition of Parameter/s")
        self.setFixedWidth(gui_width)
        self.setFixedHeight(gui_height)

        self.layout = QGridLayout()

        # create a QTableWidget parameters
        self.table_parameters = QTableWidget()

        # adding header to the table
        header_horizontal_parameters = [
            "Name",
            "Description",
            "Kind",
            "Expression",
            "Name Port",
        ]
        self.table_parameters.setColumnCount(len(header_horizontal_parameters))

        self.header_vertical_parameters = ["parameter"]
        self.table_parameters.setHorizontalHeaderLabels(header_horizontal_parameters)
        self.table_parameters.setVerticalHeaderLabels(self.header_vertical_parameters)

        # adjust size columns of horizontal header
        for i, _ in enumerate(header_horizontal_parameters):
            self.table_parameters.setColumnWidth(i, 150)

        self.button_add_parameter = QPushButton("Add parameter")
        self.button_add_parameter.clicked.connect(self.new_parameter)

        self.button_delete_parameter = QPushButton("Remove selected parameter")
        self.button_delete_parameter.clicked.connect(self.delete_parameter)

        self.button_clear_all = QPushButton("Clear All")
        self.button_clear_all.clicked.connect(self.clear_all)

        self.button_next = QPushButton("Next >")
        self.button_next.clicked.connect(self.next_page)

        self.button_prev = QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        self.layout.addWidget(self.table_parameters, 1, 0, 1, 3)
        self.layout.addWidget(self.button_clear_all, 0, 1)
        self.layout.addWidget(self.button_add_parameter, 0, 2, Qt.AlignTop)
        self.layout.addWidget(self.button_delete_parameter, 0, 3, Qt.AlignTop)

        self.layout.addWidget(self.button_next, 4, 3)
        self.layout.addWidget(self.button_prev, 4, 2)

        # resume session
        self.label_session = QLabel()
        self.layout.addWidget(self.label_session, 3, 1)

        # create navigation list
        self.comboBox = QComboBox()
        self.comboBox.addItems(gui_pages)
        self.comboBox.setCurrentText("add_parameter_page")

        # There is an alternate signal to send the text.
        self.comboBox.currentTextChanged.connect(self.text_changed)
        self.layout.addWidget(self.comboBox, 4, 1)

        self.setLayout(self.layout)

        self.help = Help(self.layout, 3, 3)
        self.table_parameters.cellClicked.connect(self.update_help)

    def update_help(self):
        """This function updates the Help object through its update_fields method.
        A text, a description and an example are prepared to be passed to the abovementioned method.
        """
        example = ""
        col = self.table_parameters.currentColumn()

        if col is not None:
            text = self.table_parameters.horizontalHeaderItem(col).text()
            print(f"col:{col},text:{text},selection:Parameter ")

            self.layout.itemAt(self.layout.count() - 1).widget().show()
            if col == 0:
                description = "Choose a name for the Parameter "

            elif col == 1:
                description = "Choose a set of words to describe the Parameter."

            elif col == 2:
                description = "Choose what is the kind of your Parameter."
                example = """It could be one of the following list:
                \n- scalar-field
                \n- vector-field
                \n- tensor-field"""

            elif col == 3:
                description = "Choose an expression for your Parameter ."

            elif col == 4:
                description = "This is the port at which the parameter is tied to."

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
        self.comboBox.setCurrentText("add_parameter_page")
        if not check_black_listed_words(self, self.table_parameters, "Parameters"):
            self.switch_window.emit(page)
            self.hide()

    def update_page(self, gui):
        """This function manages the update of the current page.

        Args:
            gui (obj): the GUI object

        """
        for row in range(self.table_parameters.rowCount()):
            comboBox = self.table_parameters.cellWidget(row, 2)
            comboBox.clear()
            if "domain" in self.session.keys() and self.session["domain"] == "Segment":
                comboBox.addItems(["scalar-field"])
            else:
                comboBox.addItems(["scalar-field", "vector-field", "tensor-field"])

        table_ports = gui.add_port_page.table_ports

        s = ""

        def itemToString(s, table, type_table):
            cols = [0]
            rows = table.rowCount()
            s += f"{type_table}: "
            if "ontrol" in type_table:
                cols = [1, 3]
            for row in range(rows):
                for col in cols:
                    item = table.item(row, col)
                    if item is not None:
                        s += item.text()
                    if row < rows - 1:
                        s += ", "
                    else:
                        s += "\n"
            return s

        s = itemToString(s, table_ports, "Ports")

        self.label_session.setText(s)

        if (
            "read_from_file" in self.session.keys()
            and not self.session["add_parameter_page"]["loaded_from_file"]
        ):
            self.load_session_from_file()

    def load_session_from_file(self):
        if (
            "parameters"
            in self.session["read_from_file"]["dict"]["add_parameter_page"].keys()
        ):
            type_to_index = {
                "scalar-field": 0,
                "vector-field": 1,
                "tensor-field": 2,
            }

            self.table_parameters.setRowCount(0)
            row = 0
            for parameter in self.session["read_from_file"]["dict"][
                "add_parameter_page"
            ]["parameters"]:
                self.new_parameter()
                for col, param in enumerate(parameter):
                    if col != 2:
                        self.table_parameters.setItem(row, col, QTableWidgetItem(param))
                    else:
                        index = type_to_index[param]
                        self.table_parameters.cellWidget(row, col).setCurrentIndex(
                            index
                        )
                row += 1
        self.session["add_parameter_page"]["loaded_from_file"] = True

    def next_page(self):
        """This function emits the signal to navigate to the next page."""
        if not check_black_listed_words(self, self.table_parameters, "Parameters"):
            self.switch_window.emit("add_control_port_page")
            self.hide()

    def previous_page(self):
        """This funcion emits the signal to navigate to the prvious page."""
        if not check_black_listed_words(self, self.table_parameters, "Parameters"):
            self.switch_window.emit("add_port_page")
            self.hide()

    def choice_clicked(self, text):
        """This function is responsible of the Help object updates.

        Args:
            text (str): the name of the selected column.
        """

        def make_update():
            print(text)
            description = ""
            example = ""

            if text == "Kind":
                description = "Choose what is the kind of your Parameter."
                example = """It could be one of the following list:
                \n- scalar-field
                \n- vector-field
                \n- tensor-field"""

            self.help.updateFields(text, description, example)

        return make_update

    def new_parameter(self):
        """This function adds 1 row in the table for parameter"""
        count = self.table_parameters.rowCount()
        self.table_parameters.insertRow(count)
        self.header_vertical_parameters += ["parameter"]
        self.table_parameters.setVerticalHeaderLabels(self.header_vertical_parameters)
        parameter_choice_kind = QComboBox()
        if "domain" in self.session.keys() and self.session["domain"] == "Segment":
            parameter_choice_kind.addItems(["scalar-field"])
        else:
            parameter_choice_kind.addItems(
                ["scalar-field", "vector-field", "tensor-field"]
            )
        parameter_choice_kind.textHighlighted.connect(self.choice_clicked("Kind"))
        self.table_parameters.setCellWidget(count, 2, parameter_choice_kind)

    def delete_parameter(self):
        """This function removes 1 row from the table"""
        if len(self.header_vertical_parameters) > 1:
            self.header_vertical_parameters.pop()
            self.table_parameters.setVerticalHeaderLabels(
                self.header_vertical_parameters
            )

            self.table_parameters.removeRow(self.table_parameters.currentRow())

        else:
            print("not enough element to delete!")

    def clear_all(self):
        """This function removes all the rows from the table."""
        self.table_parameters.setRowCount(0)
        self.new_parameter()
