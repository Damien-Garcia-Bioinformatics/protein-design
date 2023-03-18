#! /usr/bin/env python3

import sys
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QPushButton, QWidget, QApplication, QGridLayout, QLabel, QLineEdit, QFileDialog

class Window(QWidget):
    def __init__(self):
        super().__init__()
        #setting up the layout of the window
        layout = QGridLayout()
        self.setWindowTitle("4ORSA")
        self.setLayout(layout)

        generate = QLabel("Generate")
        layout.addWidget(generate, 0,0)
        new = QLabel("New")
        layout.addWidget(new,0,1)
        sequences = QLabel("Sequences")
        layout.addWidget(sequences,0,2)
        file_pdb = QLabel("PDB File :")
        layout.addWidget(file_pdb,2,0)
        thr = QLabel("Mutation Mask Setup :")
        layout.addWidget(thr,3,0)
        rdm = QLabel("Should we mutate inside (Y)\nor outside (N) the\nmasked sequence?")
        layout.addWidget(rdm,4,0)

        #Setting up inputs
        #idéalement un navigateur ici mais marche pas
        #self.input1 = QFileDialog()
        self.input1 = QLineEdit("PDB File Name")
        layout.addWidget(self.input1, 2,1)
        self.minmask = QLineEdit("minimum")
        layout.addWidget(self.minmask,3,1)
        self.maxmask = QLineEdit("maximum")
        layout.addWidget(self.maxmask,3,2)
        
        #establishing the action when the form is complete
        btn1 = QPushButton("Launch")
        btn1.clicked.connect(self.initialize)
        layout.addWidget(btn1,4,2)

        #list of inputs for random sequence generation
        self.maskStatus = QLineEdit("Y/N")
        layout.addWidget(self.maskStatus,4,1)
        

    def initialize(self):
        #calls the AlignClass library to construct and show the graph
        #self.input1.getOpenFileName(caption='Select PDB File')[0] si input1 QFileDialogue
        
        with open("param.txt", 'w+') as file :
            file.write(f"{self.input1.text()} {int(self.minmask.text())} {int(self.maxmask.text())} {self.maskStatus.text()}")
