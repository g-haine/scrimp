# from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QVBoxLayout
# app = QApplication([])
# window = QWidget()
# layout = QVBoxLayout()
# layout.addWidget(QPushButton('Top'))
# layout.addWidget(QPushButton('Bottom'))
# window.setLayout(layout)
# window.show()
# app.exec()

"""T9 application based on PyQt5."""

import os
import os.path as op
import platform
from PyQt5.QtWidgets import (QApplication, QWidget, QPushButton, QLabel, QLineEdit, QGridLayout, QMessageBox)
import signal
import subprocess
import sys


class Code_Generator_gui:

    # Create an instance of QApplication
    app = QApplication(sys.argv)

    # Retireving screen information
    screen = app.primaryScreen()
    screen_size = screen.availableSize()
    rect = screen.availableGeometry()

    # Create an instance of your application's GUI
    window = QWidget()

    # Fixing layout
    layout = QGridLayout()

    # # Add widget
    # button_calibration = QPushButton(window)

    # # Add widget
    # button_testing = QPushButton(window)

    # list of pids
    processes = []

    def __init__(self) -> None:

        self.window.setWindowTitle("Code generator APP for Scrimp")
        self.window.setGeometry(0, 0,self.screen_size.width(), self.screen_size.height())
        #self.window.move(60, 15)
        helloMsg = QLabel("<h1>Welcome to the Code Generator for SCRIMP !</h1>", parent=self.window)
        self.layout.addWidget(helloMsg,0,0)
        # helloMsg.move(60, 15)

        label_dpHs_name = QLabel('<font size="4"> Name for your distributed port-Hamiltonian system (dpHs)</font>')
        self.lineEdit_username = QLineEdit()
        self.lineEdit_username.setPlaceholderText('Please enter the name')
        self.layout.addWidget(label_dpHs_name, 1, 0)
        self.layout.addWidget(self.lineEdit_username, 1, 1)


        label_dpHs_type = QLabel('<font size="4"> Type of dpHS </font>')
        self.lineEdit_username = QLineEdit()
        self.lineEdit_username.setPlaceholderText('Please enter the type of the system')
        self.layout.addWidget(label_dpHs_type, 2, 0)
        self.layout.addWidget(self.lineEdit_username, 2, 1)

        # self.button_calibration.setText("Calibration")
        # self.button_calibration.setGeometry(200, 100, 300, 100)
        # self.button_calibration.clicked.connect(self.button_calibration_clicked)

        # self.button_testing.setText("T9 Online")
        # self.button_testing.setGeometry(200, 300, 300, 100)
        # self.button_testing.clicked.connect(self.button_testing_clicked)

        # Show GUI
        self.window.setLayout(self.layout)
        self.window.show()

        #  Run application's event loop (or main loop)
        sys.exit(self.app.exec_())

    def run_process(self, processes):

        for p in processes:
            if platform.system() == "Linux":
                try:
                    self.processes.append(
                        subprocess.Popen(
                            ["gnome-terminal -- python " + f"{p}"], shell=True
                        )
                    )
                except Exception as e:
                    self.processes.append(
                        subprocess.Popen(["gnome-terminal -- python3 " + p], shell=True)
                    )

            else:  # windows
                cmd = [sys.executable, p]
                print("command sent to subprocess: ", cmd[1:])
                try:
                    self.processes.append(subprocess.Popen(cmd))
                except Exception as e:
                    print("can't open the terminal!")

    def button_calibration_clicked(self):
        if platform.system() == "Linux":
            processes = [
                # op.join("./presentation","lsl_EEG_stream_reader.py"),
                op.join("./presentation", "T9_bci_calibration.py"),
                op.join("./classification", "online_T9_calibration.py"),
                # "./gui_LSL_streams_simulator.py"
            ]
        else:  # windows
            processes = [
                # op.join("./presentation","lsl_EEG_stream_reader.py"),
                op.join("presentation", "T9_bci_calibration.py"),
                op.join("classification", "online_T9_calibration.py"),
                # "gui_LSL_streams_simulator.py"
            ]

        print("Calibration Button clicked")
        self.run_process(processes)

    def button_testing_clicked(self):
        processes = [
            op.join("./classification", "online_T9_testing.py"),
            op.join("./presentation", "T9_bci_testing.py"),
        ]
        print("T9 Online Button clicked")
        self.run_process(processes)

        processes = ["./eegSimulator.py"]
        # run process for marker stream
        # run process
        self.run_process2(
            processes, nameStream="MyEEGStream", typeStream="EEG", n_chans=32, srate=500
        )
        # self.run_process2(processes,nameStream="MyMarkerStream",typeStream="Markers",n_chans=1,srate=0,amp=1,n_class=1,n_train=10,refresh_rate=60)

    def run_process2(
        self,
        processes,
        nameStream,
        typeStream,
        n_chans,
        srate,
        amp=0,
        n_class=0,
        n_train=0,
        refresh_rate=0,
    ):

        for p in processes:

            if typeStream == "EEG":
                cmd = [
                    sys.executable,
                    p,
                    "-n",
                    nameStream,
                    "-t",
                    typeStream,
                    "-nch",
                    str(n_chans),
                    "-sf",
                    str(srate),
                ]
            else:
                cmd = [
                    sys.executable,
                    p,
                    "-n",
                    nameStream,
                    "-t",
                    typeStream,
                    "-nch",
                    str(n_chans),
                    "-sf",
                    str(srate),
                    "-amp",
                    str(amp),
                    "-ncl",
                    str(n_class),
                    "-ntr",
                    str(n_train),
                    "-rf",
                    str(refresh_rate),
                    "-vT",
                    "string",
                ]

            if platform.system() == "Linux":
                try:
                    self.processes.append(
                        subprocess.Popen(
                            [
                                "gnome-terminal -- python "
                                + f"{p} -n {nameStream} -t {typeStream} -nch {n_chans} -sf {srate} "
                            ],
                            shell=True,
                        )
                    )
                except Exception as e:
                    self.processes.append(
                        subprocess.Popen(["gnome-terminal -- python3 " + p], shell=True)
                    )

            else:  # windows
                print("command sent to subprocess: ", cmd[1:])
                try:
                    self.processes.append(subprocess.Popen(cmd))
                except Exception as e:
                    print("can't open the terminal!")

    def __del__(self):
        # kills all the subprocess created
        for p in self.processes:
            os.kill(p.pid, signal.SIGTERM)


if __name__ == "__main__":
    Code_Generator_gui()