import sys
from PyQt5 import QtWidgets
from GUI import (
    home_page,
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


class Controller:
    """This class magage the navigation of the GUI."""

    def __init__(self):
        # the pages that are included in the GUI
        self.create_dphs = home_page.Window()
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
            self.add_parameter_page.show()
        elif text == "add_control_port_page":
            self.add_control_port_page.show()
        elif text == "add_fem_page":
            self.add_fem_page.show()
        elif text == "set_hamiltonian_page":
            self.set_hamiltonian_page.show()
        elif text == "add_term_page":
            self.add_term_page.show()
        elif text == "add_brick_page":
            self.add_brick_page.show()
        elif text == "add_expression_page":
            self.add_expression_page.show()
        elif text == "add_initial_value_page":
            self.add_initial_value_page.show()
        elif text == "set_time_scheme_page":
            self.set_time_scheme_page.show()
        elif text == "generate_code_page":
            self.generate_code_page.show()
        else:
            print("the emitted signal:", text)
            pass


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
