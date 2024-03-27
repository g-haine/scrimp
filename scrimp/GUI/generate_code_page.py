from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import (
    QMessageBox,
    QFileDialog,
    QComboBox,
    QLabel,
    QLineEdit,
    QGridLayout,
)
from utils.GUI import gui_pages, gui_width, gui_height


class Window(QtWidgets.QWidget):
    """This class defines the main page of the GUI and ask the details for the definition of the new distribuited port-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self, session):
        QtWidgets.QWidget.__init__(self)
        self.session = session
        self.setWindowTitle(
            "Generate the Python script for the defined Discrete Port Hamiltonian System"
        )
        self.setFixedWidth(gui_width)
        self.setFixedHeight(gui_height)

        self.layout = QGridLayout()

        label_filename = QLabel(
            '<font size="4"> Name for your script based on SCRIMP: </font>'
        )
        self.line_edit_filname = QLineEdit()
        self.line_edit_filname.setPlaceholderText(
            "Please insert the name for your script"
        )

        self.button_file_dialog = QtWidgets.QPushButton("Select Folder", self)
        self.button_file_dialog.clicked.connect(self.get_path)

        label_directory = QLabel('<font size="4"> The selected directory is: </font>')
        self.line_edit_directory = QLineEdit()
        self.line_edit_directory.setPlaceholderText(
            "Insert manually your filepath or click the button on the right to select the directory"
        )

        self.file_path = ""

        self.msg = QMessageBox()
        self.msg.setWindowTitle("Select Output:")
        self.msg.setText(
            "Select your prefered output by clicking on the related button:"
        )
        self.msg.setIcon(QMessageBox.Question)
        self.button_paraview = self.msg.addButton("Paraview", QMessageBox.NoRole)
        self.button_matplotlib = self.msg.addButton("Matplotlib", QMessageBox.YesRole)
        self.button_cancel = self.msg.addButton("Cancel", QMessageBox.RejectRole)
        self.button_paraview.clicked.connect(self.popup_button)
        self.button_matplotlib.clicked.connect(self.popup_button)

        self.button_generate = QtWidgets.QPushButton("Generate")
        self.button_generate.clicked.connect(self.generate_code)

        self.button_generate_run = QtWidgets.QPushButton("Generate and Run")
        self.button_generate_run.clicked.connect(self.show_message)

        self.button_prev = QtWidgets.QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        self.layout.addWidget(label_filename, 1, 0)
        self.layout.addWidget(self.line_edit_filname, 1, 1)
        self.layout.addWidget(label_directory, 2, 0)
        self.layout.addWidget(self.line_edit_directory, 2, 1)
        self.layout.addWidget(self.button_file_dialog, 2, 2)
        self.layout.addWidget(self.button_generate, 3, 4)
        self.layout.addWidget(self.button_generate_run, 3, 5)
        self.layout.addWidget(self.button_prev, 3, 3)

        # create navigation list
        self.comboBox = QComboBox()
        self.comboBox.addItems(gui_pages)
        self.comboBox.setCurrentText("generate_code_page")

        # There is an alternate signal to send the text.
        self.comboBox.currentTextChanged.connect(self.text_changed)
        self.layout.addWidget(self.comboBox, 3, 2)

        self.setLayout(self.layout)

    def show_message(self):
        self.msg.exec_()

    def popup_button(self, i):
        button = self.msg.clickedButton()
        text = button.text()
        print(text)
        self.switch_window.emit(text)

    def get_path(self):
        self.file_path = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        print(self.file_path)
        self.line_edit_directory.setText(self.file_path)

    def text_changed(self, page):  # s is a str
        """This function allows the navigation trhough the navigation list.
        After checking the presence of black listed words, the function hides the current page for showing the selected one.

        Args:
            page (str): the name of the page.
        """
        self.update_session()
        self.comboBox.setCurrentText("generate_code_page")
        self.switch_window.emit(page)
        self.hide()

    def generate_code(self):
        print("generate code...")
        self.update_session()
        self.switch_window.emit("generate_code")

    def previous_page(self):
        """This funcion emits the signal to navigate to the prvious page."""
        self.update_session()
        self.switch_window.emit("set_time_scheme_page")
        self.hide()

    def update_page(self):
        """This function manages the update of the current page."""
        self.line_edit_filname.setText(self.session["filename"])
        self.line_edit_directory.setText(self.session["filepath"])

    def update_session(self):
        self.session["filename"] = self.line_edit_filname.text()
        self.session["filepath"] = self.line_edit_directory.text()
