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
