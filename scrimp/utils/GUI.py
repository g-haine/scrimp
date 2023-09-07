from PyQt5.QtWidgets import (
    QHBoxLayout,
    QPushButton,
    QLineEdit,
    QGridLayout,
    QTableWidget,
    QComboBox,
    QLabel,
    QTextEdit,
)
from PyQt5.QtGui import QTextCharFormat, QFont, QTextCursor
from PyQt5 import QtCore, QtGui, QtWidgets

gui_pages = [
    "create_dphs_page",
    "set_domain_page",
    "add_state_costate_page",
    "add_port_page",
    "add_parameter_page",
    "add_control_port_page",
    "add_fem_page",
    "set_hamiltonian_page",
    "add_term_page",
    "add_brick_page",
    "add_expression_page",
    "add_initial_value_page",
    "set_time_scheme_page",
    "generate_code_page",
]

gui_width = 1700
gui_height = 600
button_size = (200, 150, 100, 40)
main_size_font = "+3"
secondary_size_font = "+2"


heading = """# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2022 Ghislain Haine
#
# See the LICENSE file in the root directory for license information.
#
# github: https://github.com/g-haine/scrimp

from scrimp import *
from scrimp.utils.mesh import set_verbose_gf
from itertools import zip_longest


"""


def text_add_loop(self, file, filename):
    loop = f"""    for state, costate, param, fem, port, control_port in zip_longest(
            states, costates, params, FEMs, ports, control_ports
        ):
            # Add a state
            if state is not None:
                {filename}.add_state(state)

        # Add its co-state        	
            if costate is not None:
                {filename}.add_costate(costate)

        # Add a Finite Element Method to the `port`
            if fem is not None:
                {filename}.add_FEM(fem)

        # Add a (possibly space-varying) parameter to the `port`
            if param is not None:
                {filename}.add_parameter(param)
                
        # Add a resistive `port`        	
            if port is not None:
                {filename}.add_port(port)

        # Add a control `port` on the boundary (Neumann, thus position='effort' - default)        	
            if control_port is not None:
                {filename}.add_control_port(control_port)"""

    file.write(loop)


def text_set_hamiltonian(self, file, filename):
    file.write(
        f"""\n\n    ## Set Hamiltonian
    {filename}.hamiltonian.set_name("{self.set_hamiltonian_page.line_edit_hamiltonian_name.text()}")"""
    )


def text_create_class(self, file, filename):
    index = self.create_dphs.comboBox_dphs_type.currentIndex()
    type_dphs = self.create_dphs.comboBox_dphs_type.itemText(index)
    file.write(
        f'    # Init the distributed port-Hamiltonian system\n    {filename} = DPHS("{type_dphs}")\n\n'
    )


def text_set_domain(self, file, filename):
    type_domain = self.set_domain_page.list_widget.currentItem().text()
    file.write(
        f"""    # Set the domain (using the built-in geometry `{type_domain}`)
    {filename}.set_domain(Domain("{type_domain}",{{"""
    )

    table = self.set_domain_page.table
    rows = table.rowCount()
    for row in range(rows):
        item = table.item(row, 0)
        if item is not None:
            text = item.text()
            if row + 1 < rows:
                file.write(f'"{text}":{table.item(row,1).text()},')
            else:
                file.write(f'"{text}":{table.item(row,1).text()}')

    file.write("}))")


def text_add_states(self, file):
    table_states = self.add_state_costate_page.table_states
    rows = table_states.rowCount()
    cols = table_states.columnCount()
    file.write(
        f"""\n\n    # Define State/s`)
    states = [\n"""
    )

    for row in range(rows):
        file.write("    State(")
        for col in range(cols):
            item = table_states.item(row, col)
            if col == 2:
                text = table_states.cellWidget(row, col).currentText()
            elif item is not None:
                text = item.text()

            file.write(f'"{text}"')

            if col + 1 < cols:
                file.write(f",")

        if row + 1 < rows:
            file.write("),\n")
        else:
            file.write(")")

    file.write("\n    ]")


def text_add_costates(self, file):
    table_costates = self.add_state_costate_page.table_costates
    rows = table_costates.rowCount()
    cols = table_costates.columnCount()
    file.write(
        f"""\n\n    # Define Co-State/s`)
    costates = [\n"""
    )

    for row in range(rows):
        file.write("    CoState(")
        for col in range(cols):
            item = table_costates.item(row, col)
            if col == 2:
                text = f"states[{row}]"
            elif col == 3:
                text = table_costates.cellWidget(row, col).currentText()

            elif item is not None:
                text = item.text()

            # wether to add "" or not
            if col not in [2, 3]:
                file.write(f'"{text}"')
            else:
                file.write(f"{text}")

            # wether to add , or not
            if col + 1 < cols:
                file.write(f",")

        if row + 1 < rows:
            file.write("),\n")
        else:
            file.write(")")

    file.write("\n    ]")


def text_add_ports(self, file):
    table_ports = self.add_port_page.table_ports
    rows = table_ports.rowCount()
    cols = table_ports.columnCount()
    file.write(
        f"""\n\n    # Define Ports/s`)
    ports = [\n"""
    )

    for row in range(rows):
        file.write("    Port(")
        for col in range(cols):
            item = table_ports.item(row, col)
            if col in [3, 5, 6]:
                text = table_ports.cellWidget(row, col).currentText()
            elif item is not None:
                text = item.text()

            if col not in [5, 6]:
                file.write(f'"{text}"')
            else:
                file.write(f"{text}")

            if col + 1 < cols:
                file.write(f",")

        if row + 1 < rows:
            file.write("),\n")
        else:
            file.write(")")

    file.write("\n    ]")


def text_add_parameters(self, file):
    table_parameters = self.add_parameter_page.table_parameters
    rows = table_parameters.rowCount()
    cols = table_parameters.columnCount()
    file.write(
        f"""\n\n    # Define parameters/s`)
    parameters = [\n"""
    )

    for row in range(rows):
        file.write("    Parameter(")
        for col in range(cols):
            item = table_parameters.item(row, col)
            if col in [2]:
                text = table_parameters.cellWidget(row, col).currentText()
            elif item is not None:
                text = item.text()

            file.write(f'"{text}"')

            if col + 1 < cols:
                file.write(f",")

        if row + 1 < rows:
            file.write("),\n")
        else:
            file.write(")")

    file.write("\n    ]")


def update_parameters_page(self):
    """This function updates the add parameter page accounting for the existing states and ports already declared."""
    table_states = self.add_state_costate_page.table_states
    rows_states = table_states.rowCount()

    table_ports = self.add_port_page.table_ports
    rows_ports = table_ports.rowCount()

    table_parameters = self.add_parameter_page.table_parameters

    for row in range(rows_states):
        item = table_states.item(row, 0)
        if item is not None:
            table_parameters.setItem(row, 4, item.clone())
            if table_parameters.rowCount() < rows_states + rows_ports:
                self.add_parameter_page.new_parameter()

    for row in range(rows_ports):
        item = table_ports.item(row, 0)
        if item is not None:
            table_parameters.setItem(row + rows_states, 4, item.clone())
            if table_parameters.rowCount() < rows_states + rows_ports:
                self.add_parameter_page.new_parameter()


def text_add_control_ports(self, file):
    table_control_ports = self.add_control_port_page.table_control_ports
    rows = table_control_ports.rowCount()
    cols = table_control_ports.columnCount()
    file.write(
        f"""\n\n    # Define control_ports/s)
    control_ports = [\n"""
    )

    for row in range(rows):
        file.write("    Control_Port(")
        for col in range(cols):
            item = table_control_ports.item(row, col)
            if col in [5, 7]:
                text = table_control_ports.cellWidget(row, col).currentText()
            elif item is not None:
                text = item.text()

            file.write(f'"{text}"')

            if col + 1 < cols:
                file.write(f",")

        if row + 1 < rows:
            file.write("),\n")
        else:
            file.write(")")

    file.write("\n    ]")


def update_control_ports_page(self):
    """This function updates the add control ports page accounting for the specific domain"""

    item = self.set_domain_page.list_widget.currentItem()
    if item is not None:
        type_domain = item.text()
    else:
        return

    table_control_ports = self.add_control_port_page.table_control_ports
    table_control_ports.setRowCount(0)
    self.add_control_port_page.new_control_port()

    if type_domain == "Rectangle":
        where = ["bottom", "right", "top", "left"]
        for row in range(4):
            table_control_ports.setItem(
                row,
                0,
                QtWidgets.QTableWidgetItem(f"Boundary control ({where[row]})"),
            )
            table_control_ports.setItem(
                row, 6, QtWidgets.QTableWidgetItem(f"{10+row}"))
            if table_control_ports.rowCount() < 4:
                self.add_control_port_page.new_control_port()

    elif type_domain == "Disck" or type_domain == "Concentric":
        table_control_ports.setItem(
            row,
            0,
            QtWidgets.QTableWidgetItem(f"Boundary control"),
        )
        table_control_ports.setItem(
            row, 6, QtWidgets.QTableWidgetItem(f"{10+row}"))

    elif type_domain == "Segment":
        where = ["left", "right"]
        for row in range(2):
            table_control_ports.setItem(
                row,
                0,
                QtWidgets.QTableWidgetItem(f"Boundary control ({where[row]})"),
            )
            table_control_ports.setItem(
                row, 6, QtWidgets.QTableWidgetItem(f"{10+row}"))
            if table_control_ports.rowCount() < 2:
                self.add_control_port_page.new_control_port()


def update_FEMs_page(self):
    """This function updates the add FEM page accounting for the existing states and ports already declared."""
    table_states = self.add_state_costate_page.table_states
    rows_states = table_states.rowCount()

    table_ports = self.add_port_page.table_ports
    rows_ports = table_ports.rowCount()

    table_control_ports = self.add_control_port_page.table_control_ports
    rows_control_ports = table_control_ports.rowCount()

    table_FEMs = self.add_fem_page.table_FEMs

    for row in range(rows_states):
        item = table_states.item(row, 0)
        if item is not None:
            table_FEMs.setItem(row, 0, item.clone())
            if table_FEMs.rowCount() < rows_states + rows_ports + rows_control_ports:
                self.add_fem_page.new_FEM()

    for row in range(rows_ports):
        item = table_ports.item(row, 0)
        if item is not None:
            table_FEMs.setItem(row + rows_states, 0, item.clone())
            if table_FEMs.rowCount() < rows_states + rows_ports + rows_control_ports:
                self.add_fem_page.new_FEM()

    for row in range(rows_control_ports):
        item = table_control_ports.item(row, 0)
        if item is not None:
            table_FEMs.setItem(row + rows_states + rows_ports, 0, item.clone())
            if table_FEMs.rowCount() < rows_states + rows_ports + rows_control_ports:
                self.add_fem_page.new_FEM()


def text_add_FEM(self, file):
    table_FEMs = self.add_fem_page.table_FEMs
    rows = table_FEMs.rowCount()
    cols = table_FEMs.columnCount()
    file.write(
        f"""\n\n    # Define FEM/s`)
    FEMs = [\n"""
    )

    for row in range(rows):
        file.write("    FEM(")
        for col in range(cols):
            item = table_FEMs.item(row, col)
            if col == 2:
                text = table_FEMs.cellWidget(row, col).currentText()
            elif item is not None:
                text = item.text()

            file.write(f'"{text}"')

            if col + 1 < cols:
                file.write(f",")

        if row + 1 < rows:
            file.write("),\n")
        else:
            file.write(")")

    file.write("\n    ]")


def text_add_terms(self, file, filename):
    table_terms = self.add_term_page.table_terms
    rows = table_terms.rowCount()
    cols = table_terms.columnCount()
    file.write(
        f"""\n\n    # Define Term/s`)
    terms = [\n"""
    )

    for row in range(rows):
        file.write("    Term(")

        for col in range(cols):
            item = table_terms.item(row, col)
            if item is not None:
                text = item.text()
            if col != 3:
                file.write(f'"{text}"')
            else:
                file.write(f'[')
                for i, region in enumerate(text.split(",")):
                    if i + 1 < len(text.split(",")):
                        file.write(f'{region},')
                    else:
                        file.write(f'{region}]')

            if col + 1 < cols:
                file.write(f",")

        if row + 1 < rows:
            file.write("),\n")
        else:
            file.write(")")

    file.write("\n    ]")

    file.write(
        f"""\n\n    for term in terms:
        {filename}.hamiltonian.add_term(term)\n"""
    )


class Help:
    """This class define an help section in the window of the GUI.
    Each time a field is selected it will show the description of the field and an example
    of how to use it.
    """

    def __init__(self, layout, row: int, col: int):
        self.textEdit_help = QTextEdit()
        self.textEdit_help.setReadOnly(True)
        self.textEdit_help.setFixedHeight(350)
        self.textEdit_help.setFixedWidth(500)
        self.label_help = QLabel("Help:")

        self.cursor = self.textEdit_help.textCursor()
        self.layout = layout
        self.layout.addWidget(self.label_help, row - 1, col)
        self.layout.addWidget(self.textEdit_help, row, col)

        self.layout.itemAt(self.layout.count() - 2).widget().hide()
        self.layout.itemAt(self.layout.count() - 1).widget().hide()

    def clear(self):
        self.textEdit_help.clear()
        self.layout.itemAt(self.layout.count() - 2).widget().hide()
        self.layout.itemAt(self.layout.count() - 1).widget().hide()

    def updateFields(self, name: str, description: str, example: str):
        """This function update the labels of the Help"

        Args:
            name (str): name of the selected field in the GUI
            description (str): description of the field.
            example (str): an example of how to use the field
        """
        self.textEdit_help.clear()
        self.layout.itemAt(self.layout.count() - 2).widget().show()
        self.layout.itemAt(self.layout.count() - 1).widget().show()

        html = None

        if example == "":
            html = f"""
                <font size = {main_size_font}> <b>{name}</b> </font>
                <br>
                <font size = {main_size_font}><b>{"_"*20}</b></font>
                <br><br>
                <font size = {main_size_font}><a>{description}</a></font>
                <br>
                <br>
                """
        else:
            html = f"""
                <font size = {main_size_font}><b>{name}</b></font>
                <br>
                <font size = {main_size_font}><b>{"_" *20}</b></font>
                <br><br>
                <font size = {main_size_font}><a>{description}</a></font>
                <br>
                <br>
                <i><font size = {secondary_size_font}><b>i.e.</b></font></i>
                <br>
                <i><font size = {secondary_size_font}><a>{example}</a></font></i>
                """

        self.textEdit_help.insertHtml(html)
