import sys
from PyQt5 import QtWidgets
from GUI import add_port_page, add_state_costate_page, home_page, set_domain_page

class Controller:

    def __init__(self):
        pass

    def show_create_DPHS(self):
        self.create_DPHS = home_page.DPHS_window()
        self.create_DPHS.switch_window.connect(self.show_set_domain)
        self.create_DPHS.show()

    def show_set_domain(self):
        previous_page = self.create_DPHS
        self.window = set_domain_page.MainWindow(previous_page)
        self.window.switch_window.connect(self.show_add_state_costate)
        self.create_DPHS.hide()
        self.window.show()

    def show_add_state_costate(self):
        previous_page = self.window
        self.window.hide()
        self.window = add_state_costate_page.MainWindow(previous_page)
        self.window.switch_window.connect(self.show_add_ports)
        self.window.show()

    def show_add_ports(self, text):
        previous_page = self.window
        self.window.hide()
        self.window = add_port_page.WindowTwo(previous_page)
        self.window.show()


def main():
    app = QtWidgets.QApplication(sys.argv)
    screen = app.primaryScreen()
    screen_size = screen.availableSize()
    rect = screen.availableGeometry()
    controller = Controller()
    controller.show_create_DPHS()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()