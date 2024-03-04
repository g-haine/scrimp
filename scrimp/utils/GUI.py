from PyQt5.QtWidgets import (
    QLineEdit,
    QTableWidget,
    QLabel,
    QTextEdit,
    QMessageBox,
)

from PyQt5.QtGui import QColor
from PyQt5 import QtWidgets
import keyword

# list of GUI's pages
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

# setting the size of the GUI windows
gui_width = 1700
gui_height = 600
button_size = (200, 150, 100, 40)
main_size_font = "+3"
secondary_size_font = "+2"

# list of the words that can't be used
black_listed_words = keyword.kwlist + [""]


# heading of each generated script
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


def update_list_variables(list_variables, table, name_table):
    """This function updates the list of the variables that can be plotted either withParaview or Matplotlib.

    Args:
        list_variables (list): list of the variable
        table (obj): the table object
        name_table (str): the name of the table
    """

    rows = table.rowCount()

    items = []
    for row in range(rows):
        if name_table == "state":
            items.append(table.item(row, 0))

        if name_table == "port":
            items.append(table.item(row, 1))
            items.append(table.item(row, 2))

        if name_table == "control_port":
            items.append(table.item(row, 1))
            items.append(table.item(row, 3))

    for item in items:
        if item is not None:
            variable_name = item.text()
            if variable_name not in list_variables and variable_name != "":
                list_variables.append(variable_name)


def check_black_listed_words(self, widget, widget_name):
    """This function checks if the field in the passed table contains a black listed word. It returns True if the there is a violation.

    Args:
        widget (obj): the PyQT5 widget on which the check has to be done
        widget_name (str): the name of the PyQT5 widget


    Returns:
        bool: True if a black listed word was found, else True.
    """

    violantion = False
    list_violations = []
    self.msg = QMessageBox()

    if isinstance(widget, QLineEdit):
        name = widget.text()
        if (
            name is not None
            and name != ""
            and name in self.session["black_listed_words"]
        ):
            print("name not valid")
            violantion = True
            widget.setStyleSheet("background-color:rgb(255,204,203)")
            list_violations.append(f"black listed word: {name}")
        else:
            widget.setStyleSheet("background-color:rgb(255,255,255)")
        self.msg.setWindowTitle(f"Warning in QLineEdit: {widget_name}")

    if isinstance(widget, QTableWidget):
        table = widget
        table_name = widget_name
        cols = table.columnCount()
        rows = table.rowCount()

        for row in range(rows):
            for col in range(cols):
                # if col is not None and row is not None:
                name = ""
                item = table.item(row, col)
                if item is not None:
                    name = item.text()

                    if (
                        name is not None
                        and name != ""
                        and name in self.session["black_listed_words"]
                    ):
                        print("name not valid")
                        violantion = True
                        table.item(row, col).setBackground(QColor(255, 204, 203))
                        list_violations.append(
                            f"Cell {row}:{col} contains black listed word: {name}"
                        )
                    else:
                        item.setBackground(QColor(255, 255, 255))

            self.msg.setWindowTitle(f"Warning in Table {table_name}")

    if len(list_violations) > 0:
        self.msg.setText("\n".join(list_violations))
        self.msg.exec_()
        print(list_violations)

    return violantion


# all the text_  functions are responsible to the generation of the desired script
def text_main(self, file, dphs):
    """This function writes the text for the 'main' part.

    Args:
        file (obj): the file indicating the generated script
        dphs (obj): the discrete port hamiltonian object of the script"""

    file.write(f"\n    return {dphs}\n\n")
    file.write(
        f"""if __name__ == "__main__":
    {dphs} = {dphs}_eq()"""
    )


def text_export_variables(self, file, dphs, export_variable):
    """This function adds the text for the exportation of the variable in Paraview.

    Args:
        file (obj): the file indicating the generated script
        dphs (obj): the discrete port hamiltonian object of the script
        export_variable (str): the name of the variable to export
    """

    if export_variable:
        file.write(f"\n    # export variable for Paraview")
        for variable in self.session["selected_variables"]:
            file.write(f"\n    {dphs}.export_to_pv('{variable}')")


def text_add_loop(self, file, dphs):
    """This function adds the text related to the main loop responsible for the genration of all the principal components of the DPHS

    Args:
        file (obj): the file indicating the generated script
        dphs (obj): the discrete port hamiltonian object of the script
    """

    loop = f"""\n\n    for state, costate, param, fem, port, control_port in zip_longest(
            states, costates, parameters, FEMs, ports, control_ports
        ):
            # Add a state
            if state is not None:
                {dphs}.add_state(state)

        # Add its co-state
            if costate is not None:
                {dphs}.add_costate(costate)

        # Add a Finite Element Method to the `port`
            if fem is not None:
                {dphs}.add_FEM(fem)

        # Add a (possibly space-varying) parameter to the `port`
            if param is not None:
                {dphs}.add_parameter(param)

        # Add a resistive `port`
            if port is not None:
                {dphs}.add_port(port)

        # Add a control `port` on the boundary (Neumann, thus position='effort' - default)
            if control_port is not None:
                {dphs}.add_control_port(control_port)"""

    file.write(loop)


def text_set_hamiltonian(self, file, dphs):
    """This function adds the text related to the definition of the Hamiltonian

    Args:
        file (obj): the file indicating the generated script
        dphs (obj): the discrete port hamiltonian object of the script
    """

    file.write(
        f"""\n\n    ## Set Hamiltonian
    {dphs}.hamiltonian.set_name("{self.set_hamiltonian_page.line_edit_hamiltonian_name.text()}")"""
    )


def text_create_class(self, file, dphs):
    """This function adds the text for the defintion of the instance of the DPHS class for the script.

    Args:
        file (obj): the file indicating the generated script
        dphs (obj): the discrete port hamiltonian object of the script
    """

    index = self.create_dphs.comboBox_dphs_type.currentIndex()
    type_dphs = self.create_dphs.comboBox_dphs_type.itemText(index)
    file.write(
        f'    # Init the distributed port-Hamiltonian system\n    {dphs} = DPHS("{type_dphs}")\n\n'
    )


def text_set_domain(self, file, dphs):
    """This function adds the text related to the definition of the domain for the dphs

    Args:
        file (obj): the file indicating the generated script
        dphs (obj): the discrete port hamiltonian object of the script
    """

    type_domain = self.set_domain_page.list_widget.currentItem().text()

    if type_domain == "Segment":
        type_domain = "Interval"

    file.write(
        f"""    # Set the domain (using the built-in geometry `{type_domain}`)
    {dphs}.set_domain(Domain("{type_domain}",{{"""
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
    """This function adds the text related to the definition of the states.

    Args:
        file (obj): the file indicating the generated script
    """

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
            if col in [3, 4]:
                if text == "":
                    text = None
                file.write(f"{text}")
            else:
                file.write(f'"{text}"')

            if col + 1 < cols:
                file.write(f",")

        if row + 1 < rows:
            file.write("),\n")
        else:
            file.write(")")

    file.write("\n    ]")


def text_add_costates(self, file):
    """This function adds the text related to the definition of the co-states.

    Args:
        file (obj): the file indicating the generated script
    """

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
    """This function adds the text related to the definition of the ports.

    Args:
        file (obj): the file indicating the generated script
    """

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

            if col < 4:
                file.write(f'"{text}"')
            else:
                if text == "":
                    text = None
                file.write(f"{text}")

            if col + 1 < cols:
                file.write(f",")

        if row + 1 < rows:
            file.write("),\n")
        else:
            file.write(")")

    file.write("\n    ]")


def text_add_parameters(self, file):
    """This function adds the text related to the definition of the parameters.

    Args:
        file (obj): the file indicating the generated script
    """

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


def update_intial_values_page(self):
    """This function updates the add parameter page accounting for the existing states and ports already declared."""

    table_initial_values = self.add_initial_value_page.table_initial_values

    table_states = self.add_state_costate_page.table_states
    rows_states = table_states.rowCount()

    # it updates only if it didn't do it yet
    if table_initial_values.rowCount() == 0:
        self.old_n_states = -1

    # it updates only if there is a change in the number of states
    if self.old_n_states != rows_states:
        self.old_n_states = rows_states

        for row in range(rows_states):
            item = table_states.item(row, 0)
            if item is not None:
                self.add_initial_value_page.new_initial_value()
                table_initial_values.setItem(row, 0, item.clone())


def update_parameters_page(self):
    """This function updates the add parameter page accounting for the existing states and ports already declared."""

    table_parameters = self.add_parameter_page.table_parameters

    table_states = self.add_state_costate_page.table_states
    rows_states = table_states.rowCount()

    # it updates only if it didn't do it yet
    if table_parameters.rowCount() == 0:
        self.old_n_states = -1

    # it updates only if there is a change in the number of states
    if self.old_n_states != rows_states:
        self.old_n_states = rows_states

        table_ports = self.add_port_page.table_ports
        rows_ports = table_ports.rowCount()

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
    """This function adds the text related to the definition of the control ports.

    Args:
        file (obj): the file indicating the generated script
    """

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

            if col not in [6, 8]:
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


def update_control_ports_page(self):
    """This function updates the add control ports page accounting for the specific domain"""

    table_control_ports = self.add_control_port_page.table_control_ports

    item = self.set_domain_page.list_widget.currentItem()
    if item is not None:
        type_domain = item.text()
    else:
        return

    # it updates only if it didn't do it yet
    if table_control_ports.rowCount() == 0:
        self.old_domain = -1

    # it updates only if there is a change in the domain
    if self.old_domain != type_domain:
        self.old_domain = type_domain

        # table_control_ports.setRowCount(0)
        self.add_control_port_page.new_control_port()

        if type_domain == "Rectangle":
            where = ["bottom", "right", "top", "left"]

            control = ["U_B", "U_R", "U_T", "U_L"]
            observation = ["Y_B", "Y_R", "Y_T", "Y_L"]
            for row in range(4):
                # set name
                table_control_ports.setItem(
                    row,
                    0,
                    QtWidgets.QTableWidgetItem(f"Boundary control ({where[row]})"),
                )
                # set name control
                table_control_ports.setItem(
                    row,
                    1,
                    QtWidgets.QTableWidgetItem(control[row]),
                )
                # set nameobservation
                table_control_ports.setItem(
                    row,
                    3,
                    QtWidgets.QTableWidgetItem(observation[row]),
                )

                table_control_ports.setItem(
                    row, 6, QtWidgets.QTableWidgetItem(f"{10+row}")
                )
                if table_control_ports.rowCount() < 4:
                    self.add_control_port_page.new_control_port()

        elif type_domain == "Disck" or type_domain == "Concentric":
            where = ["0", "1"]
            control = ["U_0", "U_1"]
            observation = ["Y_0", "Y_1"]
            # set name
            table_control_ports.setItem(
                row,
                0,
                QtWidgets.QTableWidgetItem(f"Boundary control {where[row]}"),
            )
            # set name control
            table_control_ports.setItem(
                row,
                1,
                QtWidgets.QTableWidgetItem(control[row]),
            )
            # set name observation
            table_control_ports.setItem(
                row,
                3,
                QtWidgets.QTableWidgetItem(observation[row]),
            )

            table_control_ports.setItem(row, 6, QtWidgets.QTableWidgetItem(f"{10+row}"))

        elif type_domain == "Segment":
            where = ["left", "right"]

            control = ["U_L", "U_R"]
            observation = ["Y_L", "Y_R"]
            for row in range(2):
                # set name

                table_control_ports.setItem(
                    row,
                    0,
                    QtWidgets.QTableWidgetItem(f"Boundary control ({where[row]})"),
                )
                # set name control
                table_control_ports.setItem(
                    row,
                    1,
                    QtWidgets.QTableWidgetItem(control[row]),
                )
                # set name observation
                table_control_ports.setItem(
                    row,
                    3,
                    QtWidgets.QTableWidgetItem(observation[row]),
                )

                table_control_ports.setItem(
                    row, 6, QtWidgets.QTableWidgetItem(f"{10+row}")
                )
                if table_control_ports.rowCount() < 2:
                    self.add_control_port_page.new_control_port()


def update_FEMs_page(self):
    """This function updates the add FEM page accounting for the existing states and ports already declared."""

    table_FEMs = self.add_fem_page.table_FEMs
    table_states = self.add_state_costate_page.table_states
    rows_states = table_states.rowCount()

    # it updates only if it didn't do it yet
    if table_FEMs.rowCount() == 0:
        self.old_n_states = -1

    # it updates only if there is a change in the number of states
    if self.old_n_states != rows_states:
        self.old_n_states = rows_states

        table_ports = self.add_port_page.table_ports
        rows_ports = table_ports.rowCount()

        table_control_ports = self.add_control_port_page.table_control_ports
        rows_control_ports = table_control_ports.rowCount()

        for row in range(rows_states):
            item = table_states.item(row, 0)
            if item is not None:
                self.add_fem_page.new_FEM()
                table_FEMs.setItem(row, 0, item.clone())
                # if table_FEMs.rowCount() < rows_states + rows_ports + rows_control_ports:

        for row in range(rows_ports):
            item = table_ports.item(row, 0)
            if item is not None:
                self.add_fem_page.new_FEM()
                table_FEMs.setItem(row + rows_states, 0, item.clone())
                # if table_FEMs.rowCount() < rows_states + rows_ports + rows_control_ports:

        for row in range(rows_control_ports):
            item = table_control_ports.item(row, 0)
            if item is not None:
                self.add_fem_page.new_FEM()
                table_FEMs.setItem(row + rows_states + rows_ports, 0, item.clone())
                # if table_FEMs.rowCount() < rows_states + rows_ports + rows_control_ports:


def update_expressions_page(self):
    """This function updates the add expression page accounting for the existing control ports already declared."""

    table_expressions = self.add_expression_page.table_expressions

    table_control_ports = self.add_control_port_page.table_control_ports
    rows_control_ports = table_control_ports.rowCount()

    # it updates only if it didn't do it yet
    if table_expressions.rowCount() == 0:
        self.old_n_control_ports = -1

    # it updates only if there is a change in the number of control ports
    if self.old_n_control_ports != rows_control_ports:
        self.old_n_control_ports = rows_control_ports

        for row in range(rows_control_ports):
            item = table_control_ports.item(row, 0)
            if item is not None:
                self.add_expression_page.new_expression()
                table_expressions.setItem(row, 0, item.clone())
                # if table_expressions.rowCount() < rows_control_ports:


def text_add_FEM(self, file):
    """This function adds the text related to the definition of the FEMs.

    Args:
        file (obj): the file indicating the generated script
    """

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

            if col == 1:
                file.write(f"{text}")
            else:
                file.write(f'"{text}"')

            if col + 1 < cols:
                file.write(f",")

        if row + 1 < rows:
            file.write("),\n")
        else:
            file.write(")")

    file.write("\n    ]")


def text_add_terms(self, file, dphs):
    """This function adds the text related to the definition of the terms.

    Args:
        file (obj): the file indicating the generated script
        dphs (obj): the discrete port hamiltonian object of the script
    """

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
            if col < 2:
                file.write(f'"{text}"')
            elif col == 2:
                file.write(f"[")
                for i, region in enumerate(text.split(",")):
                    if i + 1 < len(text.split(",")):
                        file.write(f"{region},")
                    else:
                        file.write(f"{region}]")
            else:
                file.write(f"{text}")

            if col + 1 < cols:
                file.write(f",")

        if row + 1 < rows:
            file.write("),\n")
        else:
            file.write(")")

    file.write("\n    ]")

    file.write(
        f"""\n\n    for term in terms:
        {dphs}.hamiltonian.add_term(term)\n"""
    )


def text_add_bricks(self, file, dphs):
    """This function adds the text related to the definition of the briks.

    Args:
        file (obj): the file indicating the generated script
        dphs (obj): the discrete port hamiltonian object of the script
    """

    table_bricks = self.add_brick_page.table_bricks
    rows = table_bricks.rowCount()
    cols = table_bricks.columnCount()
    file.write(
        f"""\n\n    # Define Bricks/s`)
    bricks = [\n"""
    )

    for row in range(rows):
        file.write("    Brick(")
        for col in range(cols):
            item = table_bricks.item(row, col)
            if col in [3, 4, 5]:
                text = table_bricks.cellWidget(row, col).currentText()
            elif item is not None:
                text = item.text()

            if col == 2:
                file.write(f"[")
                for i, region in enumerate(text.split(",")):
                    if i + 1 < len(text.split(",")):
                        file.write(f"{region},")
                    else:
                        file.write(f"{region}]")

            elif col not in [3, 4, 6]:
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

    file.write(
        f"""\n\n    for brick in bricks:
        {dphs}.add_brick(brick)\n"""
    )


def text_add_expressions(self, file, dphs):
    """This function adds the text related to the definition of the expressions.

    Args:
        file (obj): the file indicating the generated script
        dphs (obj): the discrete port hamiltonian object of the script
    """

    table_expressions = self.add_expression_page.table_expressions
    rows = table_expressions.rowCount()
    cols = table_expressions.columnCount()
    file.write(
        f"""\n\n    # Define expression/s`)
    expressions = ["""
    )

    for row in range(rows):
        for col in range(cols):
            if col != 0:
                item = table_expressions.item(row, col)
                if item is not None:
                    text = item.text()

                file.write(f'"{text}"')

        if row + 1 < rows:
            file.write(f",")

    file.write("]")

    file.write(
        f"""\n\n    for control_port, expression in zip(control_ports, expressions):
        {dphs}.set_control(control_port.get_name(), expression)\n"""
    )


def text_add_initial_values(self, file, dphs):
    """This function adds the text related to the definition of the initial values.

    Args:
        file (obj): the file indicating the generated script
        dphs (obj): the discrete port hamiltonian object of the script
    """

    table_initial_values = self.add_initial_value_page.table_initial_values
    rows = table_initial_values.rowCount()
    cols = table_initial_values.columnCount()
    file.write(f"""\n\n    # Define initial_value/s\n""")

    for row in range(rows):
        file.write(f"""    {dphs}.set_initial_value(""")
        for col in range(cols):
            item = table_initial_values.item(row, col)
            if item is not None:
                text = item.text()

            file.write(f'"{text}"')

            if col + 1 < cols:
                file.write(f",")
            else:
                file.write(f")\n")


def text_set_time_scheme(self, file, dphs):
    """This function adds the text related to the definition of the set time scheme.

    Args:
        file (obj): the file indicating the generated script
        dphs (obj): the discrete port hamiltonian object of the script
    """

    checkBox_answer = self.set_time_scheme_page.checkBox_answer

    if checkBox_answer.isChecked():
        table = self.set_time_scheme_page.table
        file.write(
            f"""\n    # Solve in time
    {dphs}.set_time_scheme("""
        )

        rows = table.rowCount()
        for row in range(rows):
            item = table.item(row, 0)
            if item is not None:
                text = item.text()
                if text in ["pc_type", "ts_type", "ksp_type"]:
                    file.write(f'{text}="{table.item(row,1).text()}",')
                elif row + 1 < rows:
                    file.write(f"{text}={table.item(row,1).text()},")
                else:
                    file.write(f"{text}={table.item(row,1).text()}")

        file.write(")\n")


def text_solve(self, file, dphs):
    """This function adds the text related to the definition of the solver.

    Args:
        file (obj): the file indicating the generated script
        dphs (obj): the discrete port hamiltonian object of the script
    """

    file.write(f"""\n    # Solve\n    {dphs}.solve()""")


def text_plot(self, file, dphs):
    """This function adds the text related to the definition of plotting part.

    Args:
        file (obj): the file indicating the generated script
        dphs (obj): the discrete port hamiltonian object of the script
    """

    file.write(
        f"""\n\n    # Plot the Hamiltonian with the power supplied at the boundary
    {dphs}.plot_Hamiltonian(save_figure=True)\n"""
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
