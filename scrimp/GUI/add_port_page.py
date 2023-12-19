from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import (
    QPushButton,
    QGridLayout,
    QTableWidget,
    QTableWidgetItem,
    QComboBox,
)
from PyQt5.QtCore import Qt
from utils.GUI import (
    gui_pages,
    gui_width,
    gui_height,
    Help,
    check_black_listed_words,
    update_list_variables,
)


class Window(QtWidgets.QWidget):
    """This class defines the add port and coport page of the GUI and asks to insert
    the ports and co-port realted to the new distribuited port-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self, session):
        QtWidgets.QWidget.__init__(self)
        self.session = session

        self.setWindowTitle("Definition of Port/s")
        self.setFixedWidth(gui_width)
        self.setFixedHeight(gui_height)

        self.layout = QGridLayout()

        # create a QTableWidget ports
        self.table_ports = QTableWidget()

        # adding header to the table
        header_horizontal_ports = [
            "Name",
            "Flow",
            "Effort",
            "Kind",
            "Mesh ID",
            "Algebraic",
            "Substituted",
            "Region",
        ]
        self.table_ports.setColumnCount(len(header_horizontal_ports))

        self.header_vertical_ports = ["port"]
        self.table_ports.setHorizontalHeaderLabels(header_horizontal_ports)
        self.table_ports.setVerticalHeaderLabels(self.header_vertical_ports)

        # adjust size columns of horizontal header
        for i, _ in enumerate(header_horizontal_ports):
            self.table_ports.setColumnWidth(i, 150)

        self.button_add_port = QPushButton("Add port")
        self.button_add_port.clicked.connect(self.new_port)

        self.button_delete_port = QPushButton("Remove selected port")
        self.button_delete_port.clicked.connect(self.delete_port)

        self.button_clear_all = QPushButton("Clear All")
        self.button_clear_all.clicked.connect(self.clear_all)

        self.button_next = QPushButton("Next >")
        self.button_next.clicked.connect(self.next_page)

        self.button_prev = QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        self.layout.addWidget(self.table_ports, 1, 0, 1, 3)
        self.layout.addWidget(self.button_clear_all, 0, 1)
        self.layout.addWidget(self.button_add_port, 0, 2, Qt.AlignTop)
        self.layout.addWidget(self.button_delete_port, 0, 3, Qt.AlignTop)

        self.layout.addWidget(self.button_next, 4, 3)
        self.layout.addWidget(self.button_prev, 4, 2)

        # create navigation list
        self.comboBox = QComboBox()
        self.comboBox.addItems(gui_pages)
        self.comboBox.setCurrentText("add_port_page")

        # There is an alternate signal to send the text.
        self.comboBox.currentTextChanged.connect(self.text_changed)
        self.layout.addWidget(self.comboBox, 4, 1)

        self.setLayout(self.layout)

        self.help = Help(self.layout, 3, 3)
        self.table_ports.cellClicked.connect(self.update_help)

    def update_help(self):
        """This function updates the Help object through its update_fields method.
        A text, a description and an example are prepared to be passed to the abovementioned method.
        """
        example = ""
        col = self.table_ports.currentColumn()

        if col is not None:
            text = self.table_ports.horizontalHeaderItem(col).text()
            print(f"col:{col},text:{text},selection:Port")

            self.layout.itemAt(self.layout.count() - 1).widget().show()
            if col == 0:
                description = "Choose a name for the Port"

            elif col == 1:
                description = "Choose the name of the Flow variable."

            elif col == 2:
                description = "Choose the name of the Effort variable."

            elif col == 3:
                description = """Choose what is the kind of your state.\nIt could be one of the following list:
                \n- scalar-field
                \n- vector-field
                \n- tensor-field"""
                example = "In 1D everything must be scalar-field."

            elif col == 4:
                description = (
                    "Choose which is the ID for the mesh that interest your Port."
                )
                example = "Default is 0."

            elif col == 5:
                description = "if `False`, the flow variable will be derivated in time at resolution"
                example = "Default is True."

            elif col == 6:
                description = "if `True`, the constitutive relation is substituted and there is only a getfem variable for the effort"
                example = "Default is False."

            elif col == 7:
                description = "the integer identifying the region in the mesh where the port belong, useful for boundary ports"
                example = "Default is None."

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
        self.comboBox.setCurrentText("add_port_page")
        if not check_black_listed_words(self, self.table_ports, "Ports"):
            update_list_variables(self.session["variables"], self.table_ports, "port")
            self.switch_window.emit(page)
            self.hide()

    def update_page(self):
        for row in range(self.table_ports.rowCount()):
            comboBox = self.table_ports.cellWidget(row, 3)
            comboBox.clear()
            if "domain" in self.session.keys() and self.session["domain"] == "Segment":
                comboBox.addItems(["scalar-field"])
            else:
                comboBox.addItems(["scalar-field", "vector-field", "tensor-field"])

    def next_page(self):
        """This function emits the signal to navigate to the next page."""
        if not check_black_listed_words(self, self.table_ports, "Ports"):
            update_list_variables(self.session["variables"], self.table_ports, "port")
            self.switch_window.emit("add_parameter_page")
            self.hide()

    def previous_page(self):
        """This funcion emits the signal to navigate to the prvious page."""
        if not check_black_listed_words(self, self.table_ports, "Ports"):
            update_list_variables(self.session["variables"], self.table_ports, "port")
            self.switch_window.emit("add_state_costate_page")
            self.hide()

    def choice_clicked(self, text):
        def foo():
            print(text)
            description = ""
            example = ""

            if text == "Algebraic":
                description = "if `False`, the flow variable will be derivated in time at resolution"
                example = "Default is True."

            elif text == "Kind":
                description = "Choose what is the kind of your Control Port."
                example = """It could be one of the following list:
                \n- scalar-field
                \n- vector-field
                \n- tensor-field"""

            elif text == "Substituted":
                description = "It is a boolean that defines whether to substitute the variable. Defaults to False"
                example = "Default is False"

            self.help.updateFields(text, description, example)

        return foo

    def new_port(self):
        """This function adds 1 row in the table for port"""
        count = self.table_ports.rowCount()
        self.table_ports.insertRow(count)
        self.header_vertical_ports += ["port"]
        self.table_ports.setVerticalHeaderLabels(self.header_vertical_ports)

        # create table to add in cell of table
        port_choice_algebraic = QComboBox()
        port_choice_algebraic.addItems(["True", "False"])
        port_choice_algebraic.textHighlighted.connect(self.choice_clicked("Algebraic"))
        self.table_ports.setCellWidget(count, 5, port_choice_algebraic)

        port_choice_kind = QComboBox()
        if "domain" in self.session.keys() and self.session["domain"] == "Segment":
            port_choice_kind.addItems(["scalar-field"])
        else:
            port_choice_kind.addItems(["scalar-field", "vector-field", "tensor-field"])
        port_choice_kind.textHighlighted.connect(self.choice_clicked("Kind"))
        self.table_ports.setCellWidget(count, 3, port_choice_kind)

        port_choice_substituted = QComboBox()
        port_choice_substituted.addItems(["False", "True"])
        port_choice_substituted.textHighlighted.connect(
            self.choice_clicked("Substituted")
        )
        self.table_ports.setCellWidget(count, 6, port_choice_substituted)

        # set defaults
        # mesh_id
        new_value = QTableWidgetItem("0")
        self.table_ports.setItem(count, 4, new_value)

        for i in range(self.table_ports.columnCount()):
            if i not in [3, 4, 5, 6]:
                # others
                new_value = QTableWidgetItem("")
                self.table_ports.setItem(count, i, new_value)

    def delete_port(self):
        """This function removes 2 rows in the table (1 for port, 1 for co-port)"""
        if len(self.header_vertical_ports) > 1:
            self.header_vertical_ports.pop()
            self.table_ports.setVerticalHeaderLabels(self.header_vertical_ports)

            self.table_ports.removeRow(self.table_ports.currentRow())

        else:
            print("not enough element to delete!")

    def clear_all(self):
        self.table_ports.setRowCount(0)
        self.new_port()
