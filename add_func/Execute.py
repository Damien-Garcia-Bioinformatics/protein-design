#! /usr/bin/env python3
import sys
from PyQt6.QtWidgets import QApplication
import Interface
from Interface import Window


#app launch
app = QApplication(sys.argv)
window = Window()
window.show()
sys.exit(app.exec())

