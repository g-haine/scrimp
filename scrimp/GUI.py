import re
import sys
from PyQt5 import QtWidgets
import json
from GUI import (
    welcome_page,
    load_page,
    create_dphs_page,
    set_domain_page,
    add_state_costate_page,
    add_port_page,
    add_parameter_page,
    add_control_port_page,
    add_fem_page,
    set_hamiltonian_page,
    add_term_page,
    add_brick_page,
    add_expression_page,
    add_initial_value_page,
    set_time_scheme_page,
    generate_code_page,
    export_variable_page,
)
from utils.GUI import (
    heading,
    black_listed_words,
    text_create_class,
    text_set_domain,
    text_add_states,
    text_add_costates,
    text_add_ports,
    text_add_parameters,
    update_parameters_page,
    text_add_control_ports,
    update_control_ports_page,
    text_add_FEM,
    update_FEMs_page,
    text_add_loop,
    text_set_hamiltonian,
    text_add_terms,
    text_add_bricks,
    update_expressions_page,
    text_add_expressions,
    text_add_initial_values,
    text_set_time_scheme,
    text_main,
    text_solve,
    text_plot,
    update_intial_values_page,
    text_export_variables,
)
import os


class Controller:
    """This class magage the navigation of the GUI."""

    def __init__(self):
        self.session = {}
        self.session["variables"] = []
        self.session["black_listed_words"] = black_listed_words
        self.session["filename"] = ""
        self.session["filepath"] = ""
        self.session["create_dphs_page"] = {"loaded_from_file": False}
        self.session["set_domain_page"] = {"loaded_from_file": False}
        self.session["add_state_costate_page"] = {"loaded_from_file": False}
        self.session["add_port_page"] = {"loaded_from_file": False}
        self.session["add_parameter_page"] = {"loaded_from_file": False}
        self.session["add_control_port_page"] = {"loaded_from_file": False}
        self.session["add_fem_page"] = {"loaded_from_file": False}
        self.session["set_hamiltonian_page"] = {"loaded_from_file": False}
        self.session["add_fem_page"] = {"loaded_from_file": False}
        self.session["add_term_page"] = {"loaded_from_file": False}
        self.session["add_brick_page"] = {"loaded_from_file": False}
        self.session["add_expression_page"] = {"loaded_from_file": False}
        self.session["add_initial_value_page"] = {"loaded_from_file": False}
        self.session["set_time_scheme_page"] = {"loaded_from_file": False}
        self.session["generate_code_page"] = {"loaded_from_file": False}
        self.session["auto_save"] = False
        # the pages that are included in the GUI
        self.welcome_page = welcome_page.Window(self.session)
        self.load_page = load_page.Window(self.session)
        self.create_dphs = create_dphs_page.Window(self.session)
        self.set_domain_page = set_domain_page.Window(self.session)
        self.add_port_page = add_port_page.Window(self.session)
        self.add_state_costate_page = add_state_costate_page.Window(self.session)
        self.add_parameter_page = add_parameter_page.Window(self.session)
        self.add_control_port_page = add_control_port_page.Window(self.session)
        self.add_fem_page = add_fem_page.Window(self.session)
        self.set_hamiltonian_page = set_hamiltonian_page.Window(self.session)
        self.add_term_page = add_term_page.Window(self.session)
        self.add_brick_page = add_brick_page.Window(self.session)
        self.add_expression_page = add_expression_page.Window(self.session)
        self.add_initial_value_page = add_initial_value_page.Window(self.session)
        self.set_time_scheme_page = set_time_scheme_page.Window(self.session)
        self.generate_code_page = generate_code_page.Window(self.session)
        self.export_variable_page = export_variable_page.Window(self.session)

        # the connections
        self.welcome_page.switch_window.connect(self.show_window)
        self.load_page.switch_window.connect(self.show_window)
        self.create_dphs.switch_window.connect(self.show_window)
        self.set_domain_page.switch_window.connect(self.show_window)
        self.add_state_costate_page.switch_window.connect(self.show_window)
        self.add_port_page.switch_window.connect(self.show_window)
        self.add_parameter_page.switch_window.connect(self.show_window)
        self.add_control_port_page.switch_window.connect(self.show_window)
        self.add_fem_page.switch_window.connect(self.show_window)
        self.set_hamiltonian_page.switch_window.connect(self.show_window)
        self.add_term_page.switch_window.connect(self.show_window)
        self.add_brick_page.switch_window.connect(self.show_window)
        self.add_expression_page.switch_window.connect(self.show_window)
        self.add_initial_value_page.switch_window.connect(self.show_window)
        self.set_time_scheme_page.switch_window.connect(self.show_window)
        self.generate_code_page.switch_window.connect(self.show_window)
        self.export_variable_page.switch_window.connect(self.show_window)

    def show_window(self, text: str):
        """This function according to the variable text shows the corresponding window of gui.


        Args:
            text (str): name of the page where the user wants to navigate.
        """
        if text == "welcome_page":
            self.welcome_page.update_page()
            self.welcome_page.show()
        elif text == "load_page":
            self.load_page.update_page()
            self.load_page.show()
        elif text == "create_dphs_page":
            self.create_dphs.update_page()
            self.create_dphs.show()
        elif text == "set_domain_page":
            self.set_domain_page.update_page()
            self.set_domain_page.show()
        elif text == "add_state_costate_page":
            self.add_state_costate_page.update_page()
            self.add_state_costate_page.show()
        elif text == "add_port_page":
            self.add_port_page.update_page()
            self.add_port_page.show()
        elif text == "add_parameter_page":
            update_parameters_page(self)
            self.add_parameter_page.update_page(self)
            self.add_parameter_page.show()
        elif text == "add_control_port_page":
            update_control_ports_page(self)
            self.add_control_port_page.update_page()
            self.add_control_port_page.show()
        elif text == "add_fem_page":
            update_FEMs_page(self)
            self.add_fem_page.update_page()
            self.add_fem_page.show()
        elif text == "set_hamiltonian_page":
            self.set_hamiltonian_page.update_page()
            self.set_hamiltonian_page.show()
        elif text == "add_term_page":
            self.add_term_page.update_page(self)
            self.add_term_page.show()
        elif text == "add_brick_page":
            self.add_brick_page.update_page(self)
            self.add_brick_page.show()
        elif text == "add_expression_page":
            update_expressions_page(self)
            self.add_expression_page.update_page()
            self.add_expression_page.show()
        elif text == "add_initial_value_page":
            update_intial_values_page(self)
            self.add_initial_value_page.update_page()
            self.add_initial_value_page.show()

        elif text == "set_time_scheme_page":
            self.set_time_scheme_page.update_page()
            self.set_time_scheme_page.show()

        elif text == "generate_code_page":
            self.generate_code_page.update_page()
            self.generate_code_page.show()
        elif text == "generate_code":
            self.generate_code()
        elif text == "Matplotlib":
            filename, dphs = self.generate_code()
            self.run_code(text, filename, dphs)
        elif text == "Paraview":
            self.export_variable_page.update_page()
            self.export_variable_page.show()
        elif text == "start_Paraview":
            filename, dphs = self.generate_code(True)
            self.run_code(text, None, None)
        else:
            print("the emitted signal:", text)
            pass

        if self.session["auto_save"]:
            self.save_session()

    def save_session(self):

        # create_dphs_page
        self.session["create_dphs_page"][
            "line_edit_dphs_name"
        ] = self.create_dphs.line_edit_dphs_name.text()

        index = self.create_dphs.comboBox_dphs_type.currentIndex()
        type_dphs = self.create_dphs.comboBox_dphs_type.itemText(index)
        self.session["create_dphs_page"]["comboBox_dphs_type"] = type_dphs

        # set domain page
        if self.set_domain_page.list_widget.currentItem() is not None:
            type_domain = self.set_domain_page.list_widget.currentItem().text()
            self.session["set_domain_page"]["list_widget"] = type_domain
        l = []
        table = self.set_domain_page.table
        for row in range(table.rowCount()):
            item = table.item(row, 1)
            if item is not None:
                text = item.text()
                l.append(text)

        if len(l) > 0:
            self.session["set_domain_page"]["parameters"] = l

        # add_state_page
        table_states = self.add_state_costate_page.table_states
        rows = table_states.rowCount()
        cols = table_states.columnCount()
        l = []
        for row in range(rows):
            l_row = []
            text = ""
            for col in range(cols):
                item = table_states.item(row, col)
                if col == 2:
                    cell = table_states.cellWidget(row, col)
                    if cell is not None:
                        text = cell.currentText()
                elif item is not None:
                    text = item.text()
                if col in [3, 4]:
                    if text == "":
                        text = None
                l_row.append(text)
            l.append(l_row)

        if len(l) > 0:
            self.session["add_state_costate_page"]["states"] = l

        # add_costate_page
        table_costates = self.add_state_costate_page.table_costates
        rows = table_costates.rowCount()
        cols = table_costates.columnCount()
        l = []
        for row in range(rows):
            l_row = []
            text = ""
            for col in range(cols):
                item = table_costates.item(row, col)

                if col == 2:
                    text = f"states[{row}]"
                elif col == 3:
                    cell = table_costates.cellWidget(row, col)
                    if cell is not None:
                        text = cell.currentText()
                elif item is not None:
                    text = item.text()

                l_row.append(text)
            l.append(l_row)

        if len(l) > 0:
            self.session["add_state_costate_page"]["costates"] = l

        # add_port_page
        table_ports = self.add_port_page.table_ports
        rows = table_ports.rowCount()
        cols = table_ports.columnCount()
        l = []
        for row in range(rows):
            l_row = []
            text = ""
            for col in range(cols):
                item = table_ports.item(row, col)

                if col in [3, 5, 6]:
                    cell = table_ports.cellWidget(row, col)
                    if cell is not None:
                        text = cell.currentText()
                elif item is not None:
                    text = item.text()

                l_row.append(text)
            l.append(l_row)

        if len(l) > 0:
            self.session["add_port_page"]["ports"] = l

        # add_parameter_page
        table_parameters = self.add_parameter_page.table_parameters
        rows = table_parameters.rowCount()
        cols = table_parameters.columnCount()

        l = []
        for row in range(rows):
            l_row = []
            text = ""
            for col in range(cols):
                item = table_parameters.item(row, col)
                text = ""
                if col in [2]:
                    cell = table_parameters.cellWidget(row, col)
                    if cell is not None:
                        text = cell.currentText()
                elif item is not None:
                    text = item.text()
                l_row.append(text)
            l.append(l_row)

        if len(l) > 0:
            self.session["add_parameter_page"]["parameters"] = l

        # add_control_port_page
        table_control_ports = self.add_control_port_page.table_control_ports
        rows = table_control_ports.rowCount()
        cols = table_control_ports.columnCount()

        l = []
        for row in range(rows):
            l_row = []
            text = ""
            for col in range(cols):
                item = table_control_ports.item(row, col)
                text = ""
                if col in [5, 7]:
                    cell = table_control_ports.cellWidget(row, col)
                    if cell is not None:
                        text = cell.currentText()
                elif item is not None:
                    text = item.text()
                l_row.append(text)
            l.append(l_row)

        if len(l) > 0:
            self.session["add_control_port_page"]["control_ports"] = l

        # add_fem_page
        table_FEMs = self.add_fem_page.table_FEMs
        rows = table_FEMs.rowCount()
        cols = table_FEMs.columnCount()

        l = []
        for row in range(rows):
            l_row = []
            for col in range(cols):
                item = table_FEMs.item(row, col)
                if col == 2:
                    cell = table_FEMs.cellWidget(row, col)
                    if cell is not None:
                        text = cell.currentText()
                elif item is not None:
                    text = item.text()
                l_row.append(text)
            l.append(l_row)

        if len(l) > 0:
            self.session["add_fem_page"]["fems"] = l

        # set_hamiltonian_page
        self.session["set_hamiltonian_page"][
            "name_hamiltonian"
        ] = self.set_hamiltonian_page.line_edit_hamiltonian_name.text()

        # add_term_page
        table_terms = self.add_term_page.table_terms
        rows = table_terms.rowCount()
        cols = table_terms.columnCount()

        l = []
        for row in range(rows):
            l_row = []
            l_regions = []
            for col in range(cols):
                item = table_terms.item(row, col)
                text = ""
                if item is not None:
                    text = item.text()
                elif col == 2:
                    for i, region in enumerate(text.split(",")):
                        l_regions.append(region)
                    l_row.append(l_regions)
                l_row.append(text)
            l.append(l_row)

        if len(l) > 0:
            self.session["add_term_page"]["terms"] = l

        # add_brick_page
        table_bricks = self.add_brick_page.table_bricks
        rows = table_bricks.rowCount()
        cols = table_bricks.columnCount()
        l = []
        for row in range(rows):
            l_row = []
            l_regions = []
            for col in range(cols):
                item = table_bricks.item(row, col)
                text = ""
                if col in [3, 4, 5]:
                    cell = table_bricks.cellWidget(row, col)
                    if cell is not None:
                        text = cell.currentText()
                elif item is not None:
                    text = item.text()

                if col == 2:
                    for i, region in enumerate(text.split(",")):
                        l_regions.append(region)
                    l_row.append(l_regions)
                else:
                    l_row.append(text)
            l.append(l_row)

        if len(l) > 0:
            self.session["add_brick_page"]["bricks"] = l

        # add_expression_page
        table_expressions = self.add_expression_page.table_expressions
        rows = table_expressions.rowCount()
        cols = table_expressions.columnCount()
        l = []
        for row in range(rows):
            l_row = []
            for col in range(cols):
                item = table_expressions.item(row, col)
                if item is not None:
                    text = item.text()
                    l_row.append(text)
            l.append(l_row)

        if len(l) > 0:
            self.session["add_expression_page"]["expressions"] = l

        # add_initial_value_page
        table_initial_values = self.add_initial_value_page.table_initial_values
        rows = table_initial_values.rowCount()
        cols = table_initial_values.columnCount()
        l = []
        for row in range(rows):
            l_row = []
            for col in range(cols):
                item = table_initial_values.item(row, col)
                if item is not None:
                    text = item.text()
                l_row.append(text)
            l.append(l_row)

        if len(l) > 0:
            self.session["add_initial_value_page"]["initial_values"] = l

        # set_time_scheme_page
        checkBox_answer = self.set_time_scheme_page.checkBox_answer

        if checkBox_answer.isChecked():

            table = self.set_time_scheme_page.table

            rows = table.rowCount()
            cols = table.columnCount()
            l = []
            for row in range(rows):
                l_row = []
                for col in range(cols):
                    item = table.item(row, col)
                    text = ""
                    if item is not None:
                        text = item.text()
                        l_row.append(text)
                l.append(l_row)

            if len(l) > 0:
                if self.set_time_scheme_page.list_widget.currentItem() is not None:
                    time_scheme = (
                        self.set_time_scheme_page.list_widget.currentItem().text()
                    )
                    self.session["set_time_scheme_page"]["time_scheme"] = time_scheme

            self.session["set_time_scheme_page"]["parameters"] = l

        # generate_code_page

        # store session
        filename = self.create_dphs.line_edit_filname.text()
        file_path = self.create_dphs.file_path

        if filename is None or filename == "":
            filename = "last"

        with open(os.path.join(file_path, filename + ".session.json"), "w") as f:
            json.dump(self.session, f, indent=4)

    def generate_code(self, export_variables=False):
        """This function retrieves the field from the GUI and create a python script.
        It return the filename and the name of Discrete Port Hamiltonia System declared.

        Returns:
            filename (str): name of the python script just generated.
            dphs (str): name of the dphs used in the script.
        """
        # create folder
        folder = "generated_scripts"
        if not os.path.exists(folder):
            os.makedirs(folder)

        # create file
        dphs = self.create_dphs.line_edit_dphs_name.text()
        filename = self.generate_code_page.line_edit_filname.text()
        folder = self.generate_code_page.line_edit_directory.text()
        file = open(os.path.join(folder, filename + ".py"), "w")

        # write heading and imports
        file.write(heading)

        # write func definition
        file.write(f"def {dphs}_eq():\n")

        # # write verbosity
        # file.write("    set_verbose_gf(0)\n\n")

        # create class
        text_create_class(self, file, dphs)

        # set domain and its parameters
        text_set_domain(self, file, dphs)

        # define the cariables and their dicretizations
        file.write("    ## Define the variables and their discretizations")

        # define State/s
        text_add_states(self, file)

        # define Co-State/s
        text_add_costates(self, file)

        # define Port/s
        text_add_ports(self, file)

        # # define Parameter/s
        text_add_parameters(self, file)

        # define Control Ports
        text_add_control_ports(self, file)

        # define FEMs
        text_add_FEM(self, file)

        # write loop for variables insertion (states,costates,...)
        text_add_loop(self, file, dphs)

        # set Hamiltoniam
        text_set_hamiltonian(self, file, dphs)

        # define Term/s
        text_add_terms(self, file, dphs)

        # define Brick/s
        text_add_bricks(self, file, dphs)

        # define Expression/s
        text_add_expressions(self, file, dphs)

        # define Intial Value/s
        text_add_initial_values(self, file, dphs)

        # set time scheme and its parameters
        text_set_time_scheme(self, file, dphs)

        # solve
        text_solve(self, file, dphs)

        # export solutions for ParaView
        text_export_variables(self, file, dphs, export_variables)

        # plot hamiltonian
        text_plot(self, file, dphs)

        # # define main
        text_main(self, file, dphs)

        file.close()

        print(f"created {filename}")

        return filename, dphs

    def run_code(self, text, filename, dphs):
        """This function execute the script just generated.

        Args:
            text (str): the signal emitted from generate_code_page('Matplotlib' or 'GMesh')
            filename (str): the name of the python script just generated.
            dphs (str): the name of the dphs used in the script.
        """
        try:
            if text == "Matplotlib":
                exec(f"import {filename}\n{filename}.{dphs}_eq()")

            elif text == "start_Paraview":
                print(
                    f"Selected Variable to Export: {self.session['selected_variables']}"
                )

                for variable in self.session["selected_variables"]:
                    path = self.session["path_to_export_varaible"]
                    if path is not None:
                        os.system(
                            f"paraview {os.path.join(path,variable,variable)}.pvd"
                        )
                    else:
                        os.system(f"paraview {os.path.join(variable,variable)}.pvd")
        except Exception as e:
            print(e)


def main():
    """The main function that launches the GUI."""
    app = QtWidgets.QApplication(sys.argv)
    # screen = app.primaryScreen()
    # screen_size = screen.availableSize()
    # rect = screen.availableGeometry()
    controller = Controller()
    controller.show_window("welcome_page")
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
