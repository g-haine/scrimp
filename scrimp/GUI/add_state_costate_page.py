from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import (QApplication, QWidget, QPushButton, QLabel, QLineEdit, QGridLayout, QMessageBox)

class MainWindow(QtWidgets.QWidget):

    switch_window = QtCore.pyqtSignal(str)
    previous_page = None
    next_page = None

    def __init__(self):
        QtWidgets.QWidget.__init__(self)

        self.setWindowTitle('Code generator APP for Scrimp')

        layout = QtWidgets.QGridLayout()

        self.line_edit = QtWidgets.QLineEdit()
        layout.addWidget(self.line_edit)

        self.button = QtWidgets.QPushButton('Next >')
        self.button.clicked.connect(self.switch)
        layout.addWidget(self.button)

        self.setLayout(layout)

    def switch(self):
        self.switch_window.emit(self.line_edit.text())
