import sys
from PyQt5 import QtWidgets
from GUI import (
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
)
from utils.GUI import (
    heading,
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
    text_add_expressions

)
import os


class Controller:
    """This class magage the navigation of the GUI."""

    def __init__(self):
        # the pages that are included in the GUI
        self.create_dphs = create_dphs_page.Window()
        self.set_domain_page = set_domain_page.Window()
        self.add_port_page = add_port_page.Window()
        self.add_state_costate_page = add_state_costate_page.Window()
        self.add_parameter_page = add_parameter_page.Window()
        self.add_control_port_page = add_control_port_page.Window()
        self.add_fem_page = add_fem_page.Window()
        self.set_hamiltonian_page = set_hamiltonian_page.Window()
        self.add_term_page = add_term_page.Window()
        self.add_brick_page = add_brick_page.Window()
        self.add_expression_page = add_expression_page.Window()
        self.add_initial_value_page = add_initial_value_page.Window()
        self.set_time_scheme_page = set_time_scheme_page.Window()
        self.generate_code_page = generate_code_page.Window()

        # the connections
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

    def show_window(self, text: str):
        """This function according to the variable text shows the corresponding window of gui.


        Args:
            text (str): name of the page where the user wants to navigate.
        """
        if text == "set_domain_page":
            self.set_domain_page.show()
        elif text == "add_state_costate_page":
            self.add_state_costate_page.show()
        elif text == "add_port_page":
            self.add_port_page.show()
        elif text == "create_dphs_page":
            self.create_dphs.show()
        elif text == "add_parameter_page":
            update_parameters_page(self)
            self.add_parameter_page.show()
        elif text == "add_control_port_page":
            update_control_ports_page(self)
            self.add_control_port_page.show()
        elif text == "add_fem_page":
            update_FEMs_page(self)
            self.add_fem_page.show()
        elif text == "set_hamiltonian_page":
            self.set_hamiltonian_page.show()
        elif text == "add_term_page":
            self.add_term_page.show()
        elif text == "add_brick_page":
            self.add_brick_page.show()
        elif text == "add_expression_page":
            update_expressions_page(self)
            self.add_expression_page.show()
        elif text == "add_initial_value_page":
            self.add_initial_value_page.show()
        elif text == "set_time_scheme_page":
            self.set_time_scheme_page.show()
        elif text == "generate_code_page":
            self.generate_code_page.show()
        elif text == "generate_code":
            self.generate_code()
        else:
            print("the emitted signal:", text)
            pass

    def generate_code(self):
        # create folder
        folder = "generated_scripts"
        if not os.path.exists(folder):
            os.makedirs(folder)

        # create file
        filename = self.create_dphs.line_edit_dphs_name.text()
        file = open(os.path.join(folder, filename + ".py"), "w")

        # write heading and imports
        file.write(heading)

        # write func definition
        file.write(f"def {filename}_eq():\n")

        # write verbosity
        file.write("    set_verbose_gf(0)\n\n")

        # create class
        text_create_class(self, file, filename)

        # # set domain and its parameters
        # text_set_domain(self, file, filename)

        # # define the cariables and their dicretizations
        # file.write("    ## Define the variables and their discretizations")

        # # define State/s
        # text_add_states(self, file)

        # # define Co-State/s
        # text_add_costates(self, file)

        # # define Port/s
        # text_add_ports(self, file)

        # # define Parameter/s
        # text_add_parameters(self, file)

        # define Control Ports
        text_add_control_ports(self, file)

        # # define FEMs
        # text_add_FEM(self, file)

        # # write loop for variables insertion (states,costates,...)
        # text_add_loop(self, file, filename)

        # # set Hamiltoniam
        # text_set_hamiltonian(self, file, filename)

        # # define Term/s
        # text_add_terms(self, file, filename)

        # # define Brick/s
        # text_add_bricks(self, file, filename)

        # # define Expression/s
        text_add_expressions(self, file, filename)

        file.close()
        print(f"created {filename}")


def main():
    """The main function that launches the GUI."""
    app = QtWidgets.QApplication(sys.argv)
    # screen = app.primaryScreen()
    # screen_size = screen.availableSize()
    # rect = screen.availableGeometry()
    controller = Controller()
    controller.show_window("create_dphs_page")
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
