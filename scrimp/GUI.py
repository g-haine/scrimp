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
from utils.GUI import heading
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
        index = self.create_dphs.comboBox_dphs_type.currentIndex()
        type_dphs = self.create_dphs.comboBox_dphs_type.itemText(index)
        file.write(
            f'    # Init the distributed port-Hamiltonian system\n    {filename} = DPHS("{type_dphs}")\n\n'
        )

        # set domain and its parameters
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
