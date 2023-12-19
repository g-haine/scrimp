from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import (
    QPushButton,
    QGridLayout,
    QTableWidget,
    QTableWidgetItem,
    QComboBox,
    QLabel,
)
from PyQt5.QtCore import Qt
from utils.GUI import gui_pages, gui_width, gui_height, Help, check_black_listed_words


class Window(QtWidgets.QWidget):
    """This class defines the add brick and cobrick page of the GUI and asks to insert
    the bricks and co-brick realted to the new distribuited brick-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self, session):
        QtWidgets.QWidget.__init__(self)
        self.session = session

        self.setWindowTitle("Definition of Brick/s")
        self.setFixedWidth(gui_width)
        self.setFixedHeight(gui_height)

        self.layout = QGridLayout()

        # create a QTableWidget bricks
        self.table_bricks = QTableWidget()

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

        self.button_delete_brick = QPushButton("Remove selected brick")
        self.button_delete_brick.clicked.connect(self.delete_brick)

        self.button_clear_all = QPushButton("Clear All")
        self.button_clear_all.clicked.connect(self.clear_all)

        self.button_next = QPushButton("Next >")
        self.button_next.clicked.connect(self.next_page)

        self.button_prev = QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        self.layout.addWidget(self.table_bricks, 1, 0, 1, 3)
        self.layout.addWidget(self.button_clear_all, 0, 1)
        self.layout.addWidget(self.button_add_brick, 0, 2, Qt.AlignTop)
        self.layout.addWidget(self.button_delete_brick, 0, 3, Qt.AlignTop)

        self.layout.addWidget(self.button_next, 4, 3)
        self.layout.addWidget(self.button_prev, 4, 2)

        # resume session
        self.label_session = QLabel()
        self.layout.addWidget(self.label_session, 3, 1)

        # create navigation list
        self.comboBox = QComboBox()
        self.comboBox.addItems(gui_pages)
        self.comboBox.setCurrentText("add_brick_page")

        # There is an alternate signal to send the text.
        self.comboBox.currentTextChanged.connect(self.text_changed)
        self.layout.addWidget(self.comboBox, 4, 1)

        self.setLayout(self.layout)

        self.help = Help(self.layout, 3, 3)
        self.table_bricks.cellClicked.connect(self.update_help)

    def update_help(self):
        example = ""
        col = self.table_bricks.currentColumn()

        if col is not None:
            text = self.table_bricks.horizontalHeaderItem(col).text()
            print(f"col:{col},text:{text},selection:Brick")

            self.layout.itemAt(self.layout.count() - 1).widget().show()
            if col == 0:
                description = "Choose a name for the Brick, it will be used mainly for plotting purpose."

            elif col == 1:
                description = "Insert the form in GWFL getfem language."

            elif col == 2:
                description = "Choose the the regions of mesh where the form applies. If more than one use a coma ',' to separate them with no spaces in between."
                example = "If multiple: 0,1,2"

            elif col == 3:
                description = (
                    "This is a parameter to help easy identification of linear bricks."
                )
                example = "Defaults is True."

            elif col == 4:
                description = "This is a parameter to help easy identification of matrices applied to time derivative of a variable."
                example = "The mass matrices.\n Defaults to False."

            elif col == 5:
                description = """This is a parameter to help easy identification of 'where' is the form. This serves for both the time-resolution and plotting purposes."""
                example = """ in the Dirac structure ('flow' side or 'effort' side), or in the 'constitutive' relations.
                Defaults to 'constitutive'"""

            elif col == 6:
                description = "Choose the ID of the mesh where the form applies."
                example = "Default is 0."

            self.help.updateFields(text, description, example)

        else:
            self.help.clear()
            self.layout.itemAt(self.layout.count() - 1).widget().hide()

    def text_changed(self, page):  # s is a str
        self.comboBox.setCurrentText("add_brick_page")
        if not check_black_listed_words(self, self.table_bricks, "Bricks"):
            self.switch_window.emit(page)
            self.hide()

    def update_page(self, gui):
        table_states = gui.add_state_costate_page.table_states
        table_costates = gui.add_state_costate_page.table_costates
        table_ports = gui.add_port_page.table_ports
        table_control_ports = gui.add_control_port_page.table_control_ports
        table_parameters = gui.add_parameter_page.table_parameters
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

        s = itemToString(s, table_states, "States")
        s = itemToString(s, table_costates, "Costates")
        s = itemToString(s, table_ports, "Ports")
        s = itemToString(s, table_control_ports, "Control Ports")
        s = itemToString(s, table_parameters, "Parameters")

        self.label_session.setText(s)

    def next_page(self):
        """This funciont emit the signal to navigate to the next page."""
        if not check_black_listed_words(self, self.table_bricks, "Bricks"):
            self.switch_window.emit("add_expression_page")
            self.hide()

    def previous_page(self):
        """This funcion emits the signal to navigate to the prvious page."""
        if not check_black_listed_words(self, self.table_bricks, "Bricks"):
            self.switch_window.emit("add_term_page")
            self.hide()

    def choice_clicked(self, text):
        def foo():
            print(text)
            description = ""
            example = ""

            if text == "Linear":
                description = (
                    "This is a parameter to help easy identification of linear bricks."
                )
                example = "Defaults is True."

            elif text == "dt":
                description = "This is a parameter to help easy identification of matrices applied to time derivative of a variable."
                example = "The mass matrices.\n Defaults to False."

            elif text == "Position":
                description = "A parameter to help easy identification of 'where' is the form: in the Dirac structure('flow' side or 'effort' side), or in the 'constitutive' relations. This serves for both the time-resolution and plotting purposes."
                example = "Defaults to 'constitutive"

            self.help.updateFields(text, description, example)

        return foo

    def new_brick(self):
        """This function adds 1 row in the table for brick"""
        count = self.table_bricks.rowCount()
        self.table_bricks.insertRow(count)
        self.header_vertical_bricks += ["brick"]
        self.table_bricks.setVerticalHeaderLabels(self.header_vertical_bricks)

        brick_choice_linear = QComboBox()
        brick_choice_linear.addItems(["True", "False"])
        brick_choice_linear.textHighlighted.connect(self.choice_clicked("Linear"))
        self.table_bricks.setCellWidget(count, 3, brick_choice_linear)

        brick_choice_dt = QComboBox()
        brick_choice_dt.addItems(["False", "True"])
        brick_choice_dt.textHighlighted.connect(self.choice_clicked("dt"))
        self.table_bricks.setCellWidget(count, 4, brick_choice_dt)

        brick_choice_position = QComboBox()
        brick_choice_position.addItems(["constitutive", "flow", "effort"])
        brick_choice_position.textHighlighted.connect(self.choice_clicked("Position"))
        self.table_bricks.setCellWidget(count, 5, brick_choice_position)

        # set defaults
        # mesh_id
        new_value = QTableWidgetItem("0")
        self.table_bricks.setItem(count, 6, new_value)

        for i in range(self.table_bricks.columnCount()):
            if i not in [3, 4, 5, 6]:
                # others
                new_value = QTableWidgetItem("")
                self.table_bricks.setItem(count, i, new_value)

    def delete_brick(self):
        """This function removes 2 rows in the table (1 for brick, 1 for co-brick)"""
        if len(self.header_vertical_bricks) > 1:
            self.header_vertical_bricks.pop()
            self.table_bricks.setVerticalHeaderLabels(self.header_vertical_bricks)

            self.table_bricks.removeRow(self.table_bricks.currentRow())

        else:
            print("not enough element to delete!")

    def clear_all(self):
        self.table_bricks.setRowCount(0)
        self.new_brick()
