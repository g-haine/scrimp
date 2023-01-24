from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import (QApplication, QWidget, QPushButton, QLabel, QLineEdit, QGridLayout, QMessageBox)

class DPHS_window(QtWidgets.QWidget):
    """This class defines the main page of the GUI and ask the detaitil for the definition of the new distribuited port-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal()
    previous_page = None
    next_page = None

    def __init__(self):
        QtWidgets.QWidget.__init__(self)
        self.setWindowTitle('Definition of new distributed port-Hamiltonian system')
        self.setFixedWidth(600)
        self.setFixedHeight(250)

        layout = QtWidgets.QGridLayout()

        label_dpHs_name = QLabel('<font size="4"> Name for your dpHs</font>')
        self.lineEdit_dpHs_name= QLineEdit()
        self.lineEdit_dpHs_name.setPlaceholderText('Please enter the name')
        layout.addWidget(label_dpHs_name, 1, 0)
        layout.addWidget(self.lineEdit_dpHs_name, 1, 1)


        label_dpHs_type = QLabel('<font size="4"> Type of dpHS </font>')
        self.lineEdit_dpHs_type = QLineEdit()
        self.lineEdit_dpHs_type.setPlaceholderText('Please enter the type of the system')
        layout.addWidget(label_dpHs_type, 2, 0)
        layout.addWidget(self.lineEdit_dpHs_type, 2, 1)

        self.button = QtWidgets.QPushButton('Next >')
        self.button.clicked.connect(self.nextPage)

        layout.addWidget(self.button,3,2)

        self.setLayout(layout)
        
    def nextPage(self):
        self.hide()
        self.next_page.show()
        # self.switch_window.emit()
