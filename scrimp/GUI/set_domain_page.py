from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import (
    QListWidget,
    QListWidgetItem,
    QAbstractItemView,
    QPushButton,
    QLabel,
    QLineEdit,
    QGridLayout,
    QTableWidget,
)


class Window(QtWidgets.QWidget):
    """This class defines the set domain page of the GUI and asks to define domain of the new distribuited port-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self):
        QtWidgets.QWidget.__init__(self)

        self.setWindowTitle("Define the domain for the dpHs")
        self.setFixedWidth(600)
        self.setFixedHeight(600)

        layout = QGridLayout()

        label_type_domain = QLabel(
            '<font size="4">'
            + "<p>Select a domain for your dpHs from the built-in<\p>"
            + "<p> list or define your own one:<\p>"
            + "</font>"
        )

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

        label_parameter = QLabel('<font size="4"> Parameters</font>')
        layout.addWidget(label_parameter, 3, 0)

        self.button_next = QPushButton("Next >")
        self.button_next.clicked.connect(self.next_page)

        layout.addWidget(self.button_next, 5, 2)

        self.button_prev = QPushButton("< Prev")
        self.button_prev.clicked.connect(self.previous_page)

        layout.addWidget(self.button_prev, 5, 1)

        # create a QTableWidget
        self.table = QTableWidget()
        self.table.setRowCount(2)
        self.table.setColumnCount(2)
        self.table.setGeometry(50, 100, 200, 300)

        # adding header to the table
        header_horizontal = ["Name", "Value"]

        self.table.setHorizontalHeaderLabels(header_horizontal)

        # adjust size columns of horizontal header
        for i, _ in enumerate(header_horizontal):
            self.table.setColumnWidth(i, 100)

        layout.addWidget(self.table, 4, 0)

        self.button_new = QPushButton("New")
        self.button_new.clicked.connect(self.new_rows)

        layout.addWidget(self.button_new, 3, 1)

        self.setLayout(layout)

    def new_rows(self):
        """This function adds 2 rows in the table (1 for state, 1 for co-state)"""
        count = self.table.rowCount()
        self.table.insertRow(count)

    def next_page(self):
        """This funciont emit the signal to navigate to the next page."""
        self.switch_window.emit("add_state_costate_page")
        self.hide()

    def previous_page(self):
        """This funciont emit the signal to navigate to the previous page."""
        self.switch_window.emit("create_DHPS_page")
        self.hide()
