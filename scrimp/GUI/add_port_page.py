from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import (QApplication, QWidget, QPushButton, QLabel, QLineEdit, QGridLayout, QMessageBox)


class WindowTwo(QtWidgets.QWidget):

    previous_page = None
    next_page = None

    def __init__(self):
        QtWidgets.QWidget.__init__(self)

        self.setWindowTitle('Window Two')

        layout = QtWidgets.QGridLayout()

        self.label = QtWidgets.QLabel()
        layout.addWidget(self.label)

        self.button = QtWidgets.QPushButton('Close')
        self.button.clicked.connect(self.close)

        layout.addWidget(self.button)

        self.setLayout(layout)



