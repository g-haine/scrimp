from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import (
    QHBoxLayout,
    QPushButton,
    QLineEdit,
    QGridLayout,
    QTableWidget,
    QTableWidgetItem,
    QComboBox,
)
from PyQt5.QtCore import Qt
from utils.GUI import gui_pages, gui_width, gui_height, Help


class Window(QtWidgets.QWidget):
    """This class defines the add control_port and cocontrol_port page of the GUI and asks to insert
    the control_ports and co-control_port realted to the new distribuited control_port-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self,session):
        QtWidgets.QWidget.__init__(self)
        self.session  = session

        self.setWindowTitle("Definition of Control control_port/s")
        self.setFixedWidth(gui_width)
        self.setFixedHeight(gui_height)
        # self.setGeometry(100, 100, 600, 300)

        self.layout = QGridLayout()

        # self.line_edit = QLineEdit()
        # layout.addWidget(self.line_edit)

        # create a QTableWidget control_ports
        self.table_control_ports = QTableWidget()
        # self.table_control_ports.setRowCount(1)

        # adding header to the table
        header_horizontal_control_ports = [
            "Name",
            "Name Control",
            "Description Control",
            "Name Observation",
            "Description Observation",
            "Kind",
            "Region",
            "Position",
            "Mesh ID",
        ]
        self.table_control_ports.setColumnCount(
            len(header_horizontal_control_ports))

        self.header_vertical_control_ports = ["control_port"]
        self.table_control_ports.setHorizontalHeaderLabels(
            header_horizontal_control_ports
        )
        self.table_control_ports.setVerticalHeaderLabels(
            self.header_vertical_control_ports
        )

        # adjust size columns of horizontal header
        for i, _ in enumerate(header_horizontal_control_ports):
            self.table_control_ports.setColumnWidth(i, 155)

        self.button_add_control_port = QPushButton("Add control_port")
        self.button_add_control_port.clicked.connect(self.new_control_port)

        self.button_delete_control_port = QPushButton(
            "Remove selected control_port")
        self.button_delete_control_port.clicked.connect(
            self.delete_control_port)

        self.button_clear_all = QPushButton("Clear All")
        self.button_clear_all.clicked.connect(self.clear_all)

        # layout_buttons_control_port = QHBoxLayout()

        # layout_buttons_control_port.addWidget(self.button_add_control_port)
        # layout_buttons_control_port.addWidget(self.button_delete_control_port)

        # cell_double = QTableWidget(layout_buttons_control_port)

        self.button_next = QPushButton("Next >")
        self.button_next.clicked.connect(self.next_page)

        self.button_prev = QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        self.layout.addWidget(self.table_control_ports, 1, 0, 1, 3)
        # layout.addWidget(cell_double, 1, 3)
        self.layout.addWidget(self.button_clear_all, 0, 1)
        self.layout.addWidget(self.button_add_control_port, 0, 2, Qt.AlignTop)
        self.layout.addWidget(
            self.button_delete_control_port, 0, 3, Qt.AlignTop)

        self.layout.addWidget(self.button_next, 4, 3)
        self.layout.addWidget(self.button_prev, 4, 2)

        # create navigation list
        self.comboBox = QComboBox()
        self.comboBox.addItems(gui_pages)
        self.comboBox.setCurrentText("add_control_port_page")

        # There is an alternate signal to send the text.
        self.comboBox.currentTextChanged.connect(self.text_changed)
        self.layout.addWidget(self.comboBox, 4, 1)

        self.setLayout(self.layout)

        self.help = Help(self.layout, 3, 3)
        self.table_control_ports.cellClicked.connect(self.update_help)

        # self.new_control_port()

    def update_help(self):
        example = ""
        col = self.table_control_ports.currentColumn()

        if col is not None:
            text = self.table_control_ports.horizontalHeaderItem(col).text()
            print(f"col:{col},text:{text},selection:Control Port")

            self.layout.itemAt(self.layout.count() - 1).widget().show()
            if col == 0:
                description = "Choose a name for the Control Port"

            elif col == 1:
                description = "Choose a letter and its subscript to name the Control."
                example = "U_B"
            elif col == 2:
                description = "Choose a set of words to describe the Control."
                example = "Temperature"

            elif col == 3:
                description = (
                    "Choose a letter and its subscript to name the Observation"
                )
                example = "Y_B"
            elif col == 4:
                description = "Choose a set of words to describe the Observation."
                example = "Normal Heat Flux"

            elif col == 5:
                description = "Choose what is the kind of your Control Port."
                example = """It could be one of the following list:
                \n- scalar-field
                \n- vector-field
                \n- tensor-field"""

            elif col == 6:
                description = "the integer identifying the region in the mesh where the Control Port belong, useful for boundary ports"
                example = "Default is None."

            elif col == 7:
                description = (
                    "Defines wether the Control Port is an 'effort' or a 'flow'."
                )
                example = "Default is effort."

            elif col == 8:
                description = "Choose which is the ID for the mesh that interest your Control Port."
                example = "Default is 0."

            self.help.updateFields(text, description, example)

        else:
            self.help.clear()
            self.layout.itemAt(self.layout.count() - 1).widget().hide()

    def text_changed(self, page):  # s is a str
        self.comboBox.setCurrentText("add_control_port_page")
        self.switch_window.emit(page)
        self.hide()

    def next_page(self):
        """This funciont emit the signal to navigate to the next page."""
        self.switch_window.emit("add_fem_page")
        self.hide()

    def previous_page(self):
        """This funcion emits the signal to navigate to the prvious page."""
        self.switch_window.emit("add_parameter_page")
        self.hide()

    def choice_clicked(self, text):
        def foo():
            print(text)
            description = ""
            example = ""

            if text == "Position":
                description = (
                    "Defines wether the Control Port is an 'effort' or a 'flow'."
                )
                example = "Default is effort."

            elif text == "Kind":
                description = "Choose what is the kind of your Control Port."
                example = """It could be one of the following list:
                \n- scalar-field
                \n- vector-field
                \n- tensor-field"""

            self.help.updateFields(text, description, example)

        return foo

    def new_control_port(self):
        """This function adds 1 row in the table for control_port"""
        count = self.table_control_ports.rowCount()
        self.table_control_ports.insertRow(count)
        self.header_vertical_control_ports += ["control_port"]
        self.table_control_ports.setVerticalHeaderLabels(
            self.header_vertical_control_ports
        )

        controlport_choice_kind = QComboBox()
        if "domain" in self.session.keys() and self.session["domain"] == "Segment":
            controlport_choice_kind.addItems(
                ["scalar-field"]
            )
        else:
            controlport_choice_kind.addItems(
                ["scalar-field", "vector-field", "tensor-field"]
            )
            
        controlport_choice_kind.textHighlighted.connect(
            self.choice_clicked("Kind"))
        self.table_control_ports.setCellWidget(
            count, 5, controlport_choice_kind)

        controlport_choice_poistion = QComboBox()
        controlport_choice_poistion.addItems(["effort", "flow"])
        controlport_choice_poistion.textHighlighted.connect(
            self.choice_clicked("Position"))
        self.table_control_ports.setCellWidget(
            count, 7, controlport_choice_poistion)

        # set defaults

        for i in range(9):
            if i == 6:
                # Region
                new_value = QTableWidgetItem("None")
                self.table_control_ports.setItem(count, 6, new_value)
            elif i == 8:
                # mesh_id
                new_value = QTableWidgetItem("0")
                self.table_control_ports.setItem(count, 8, new_value)
            else:
                # others
                new_value = QTableWidgetItem("")
                self.table_control_ports.setItem(count, i, new_value)

    def delete_control_port(self):
        """This function removes 2 rows in the table (1 for control_port, 1 for co-control_port)"""
        if len(self.header_vertical_control_ports) > 1:
            self.header_vertical_control_ports.pop()
            self.table_control_ports.setVerticalHeaderLabels(
                self.header_vertical_control_ports
            )

            self.table_control_ports.removeRow(
                self.table_control_ports.currentRow())

        else:
            print("not enough element to delete!")

    def clear_all(self):
        self.table_control_ports.setRowCount(0)
        self.new_control_port()
