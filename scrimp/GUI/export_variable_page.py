from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QLabel,  QGridLayout ,QCheckBox
from utils.GUI import gui_pages, gui_width, gui_height


class Window(QtWidgets.QWidget):
    """This class defines the main page of the GUI and ask the details for the definition of the new distribuited port-Hamiltonian system.

    Args:
        QtWidgets (_type_): _description_
    """

    switch_window = QtCore.pyqtSignal(str)

    def __init__(self,session):
        QtWidgets.QWidget.__init__(self)
        self.session  = session
        self.setWindowTitle("Selection desired variable to export in Paraview")
        self.setFixedWidth(int(gui_width/2))
        self.setFixedHeight(int(gui_height/2))

        self.layout = QGridLayout()

        label = QLabel('<font size="4"> List of available Variables in Getfem: </font>')
        self.label_hidden = QLabel('')       

        self.button_export = QtWidgets.QPushButton("Export")
        self.button_export.clicked.connect(self.export_selected_variables)

        self.layout.addWidget(label, 1, 0)



        self.setLayout(self.layout)


    def update_page(self):
        self.label_hidden.setText("")
        self.listCheckBox = [QCheckBox("all")]

        if len(self.session["variables"])>0 :
            self.layout.addWidget(self.listCheckBox[0],2,0)
        
        for i,v in enumerate(self.session["variables"]):
            self.listCheckBox.append(QCheckBox(v))
            self.layout.addWidget(self.listCheckBox[i+1],i+3,0)
        
        self.layout.addWidget(self.label_hidden, len(self.listCheckBox)+2, 0)
        self.layout.addWidget(self.button_export, len(self.listCheckBox)+3, 0)

    def export_selected_variables(self):
        if len(self.session["variables"])>0:
            self.session["selected_variables"] = []
            if self.listCheckBox[0].checkState(): # if all selected
                self.session["selected_variables"] =  self.session["variables"].copy()
            else:
                for i,checkbox in enumerate(self.listCheckBox[1:]):
                    if checkbox.checkState():
                        self.session["selected_variables"].append(self.session["variables"][i])
            self.switch_window.emit("start_Paraview")
            self.hide()
        else:
            self.label_hidden.setText("No Variables added")
