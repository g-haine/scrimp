import sys
from PyQt5 import QtWidgets
from GUI import add_port_page, add_state_costate_page, home_page, set_domain_page

class Controller:

    def __init__(self):
        self.create_DPHS = home_page.DPHS_window()
        self.set_domain_page = set_domain_page.MainWindow()
        self.add_port_page = add_port_page.WindowTwo()
        #self.add_state_costate_page = add_state_costate_page.MainWindow()
        

    def show_create_DPHS(self):
        #self.create_DPHS.switch_window.connect(self.show_set_domain)
        next_page = self.set_domain_page
        self.create_DPHS.next_page = next_page
        self.create_DPHS.show()

    def show_set_domain(self):
        previous_page = self.create_DPHS
        next_page = self.add_port_page
        self.set_domain_page.previous_page = previous_page
        self.set_domain_page.next_page = next_page

        self.window.switch_window.connect(self.show_add_state_costate)
        self.create_DPHS.hide()
        self.window.show()

    def show_add_state_costate(self):
        previous_page = self.set_domain_page
        self.window.hide()
        self.window.switch_window.connect(self.show_add_ports)
        self.window.show()

    def show_add_ports(self, text):
        previous_page = self.add_port_page
        self.window.hide()
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