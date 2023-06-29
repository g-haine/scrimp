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
        s1 = name + "\n"
        s2 = "_" * 20 + "\n"
        s3 = description + "\n\n"
        s4 = ""
        if example != "":
            s4 = "i.e.\n"
        s5 = example
        self.textEdit_help.setPlainText(s1 + s2 + s3 + s4 + s5)

        # set font 1st string
        format_1 = QTextCharFormat()
        format_1.setFont(QFont("Arial", 16, QFont.Bold))
        self.cursor.setPosition(0)
        self.cursor.movePosition(QTextCursor.EndOfLine, QTextCursor.KeepAnchor)
        self.cursor.mergeCharFormat(format_1)

        # skip one line/string
        self.cursor.movePosition(QTextCursor.NextBlock)

        # set font 3rd string
        format_2 = QTextCharFormat()
        format_2.setFont(QFont("Arial", 16))
        self.cursor.movePosition(QTextCursor.NextBlock)
        self.cursor.movePosition(QTextCursor.EndOfLine, QTextCursor.KeepAnchor)
        self.cursor.mergeCharFormat(format_2)

        # skip one line/string
        self.cursor.movePosition(QTextCursor.NextBlock)

        # set font 4th string
        format_3 = QTextCharFormat()
        format_3.setFont(QFont("Courier", 12, QFont.Bold, italic=True))
        self.cursor.movePosition(QTextCursor.NextBlock)
        self.cursor.movePosition(QTextCursor.EndOfLine, QTextCursor.KeepAnchor)
        self.cursor.mergeCharFormat(format_3)

        # set font 5th string
        format_3 = QTextCharFormat()
        format_3.setFont(QFont("Courier", 12, italic=True))
        self.cursor.movePosition(QTextCursor.NextBlock)
        self.cursor.movePosition(QTextCursor.EndOfLine, QTextCursor.KeepAnchor)
        self.cursor.mergeCharFormat(format_3)
