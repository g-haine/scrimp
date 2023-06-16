import sys
from PyQt5 import QtWidgets
from GUI import (
    add_control_port_page,
    add_parameter_page,
    add_port_page,
    add_state_costate_page,
    home_page,
    set_domain_page,
)


class Controller:
    """This class magage the navigation of the GUI."""

    def __init__(self):
        # the pages that are included in the GUI
        self.create_DPHS = home_page.Window()
        self.set_domain_page = set_domain_page.Window()
        self.add_port_page = add_port_page.Window()
        self.add_state_costate_page = add_state_costate_page.Window()
        self.add_parameter_page = add_parameter_page.Window()
        self.add_control_port_page = add_control_port_page.Window()

        # the connections
        self.create_DPHS.switch_window.connect(self.show_window)
        self.set_domain_page.switch_window.connect(self.show_window)
        self.add_state_costate_page.switch_window.connect(self.show_window)
        self.add_port_page.switch_window.connect(self.show_window)
        self.add_parameter_page.switch_window.connect(self.show_window)
        self.add_control_port_page.switch_window.connect(self.show_window)

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
        elif text == "create_DHPS_page":
            self.create_DPHS.show()
        elif text == "add_parameter_page":
            self.add_parameter_page.show()
        elif text == "add_control_port_page":
            self.add_control_port_page.show()
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
    controller.show_window("create_DHPS_page")
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
