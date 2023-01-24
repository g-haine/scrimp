from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import (QListWidget,QListWidgetItem,QAbstractItemView, QApplication, QWidget, QPushButton, QLabel, QLineEdit, QGridLayout, QMessageBox)

class MainWindow(QtWidgets.QWidget):

    switch_window = QtCore.pyqtSignal(str)
    previous_page = None
    next_page = None

    def __init__(self):
        QtWidgets.QWidget.__init__(self)
        
        self.setWindowTitle('Definition of the domain for the dpHs')
        self.setFixedWidth(600)
        self.setFixedHeight(250)
        
        
        layout = QtWidgets.QGridLayout()

        label_type_domain = QLabel('<font size="4">' +
         '<p>Select a domain for your dpHs from the built-in<\p>' +
         '<p> list or define your own one:<\p>' +
         '</font>')
        
        layout.addWidget(label_type_domain, 1, 0)

        # creating a QListWidget
        list_widget = QListWidget(self)
 
        # list widget items
        item1 = QListWidgetItem("Rectangle")
        item2 = QListWidgetItem("Circle")
        item3 = QListWidgetItem("Other")
 
        # adding items to the list widget
        list_widget.addItem(item1)
        list_widget.addItem(item2)
        list_widget.addItem(item3)
 
        # setting selection mode property
        list_widget.setSelectionMode(QAbstractItemView.SingleSelection)

        layout.addWidget(list_widget, 2, 1)


        label_dpHs_type = QLabel('<font size="4"> Type of dpHS </font>')
        self.lineEdit_dpHs_type = QLineEdit()
        self.lineEdit_dpHs_type.setPlaceholderText('Please enter the type of the system')
        layout.addWidget(label_dpHs_type, 3, 0)
        layout.addWidget(self.lineEdit_dpHs_type, 3, 1)

        self.button_next = QtWidgets.QPushButton('Next >')
        self.button_next.clicked.connect(self.nextPage)

        layout.addWidget(self.button_next,3,2)

        self.button_prev = QtWidgets.QPushButton('< Prev')
        self.button_prev.clicked.connect(self.previousPage)

        layout.addWidget(self.button_prev,3,1)

        self.setLayout(layout)

    def nextPage(self):
        #self.switch_window.emit("next")
        self.hide()
        self.next_page.show()

    def previousPage(self):
        self.hide()
        self.previous_page.show()
