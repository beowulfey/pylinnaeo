import logging
import sys

from PyQt5.QtCore import QFile
from PyQt5.QtWidgets import QApplication

from linnaeo.main import Linnaeo
from classes.views import LinnaeoApp


def main():
    logging.basicConfig(level=logging.DEBUG)  # , format="%(asctime)s:%(levelname)s:%(message)s")
    app = LinnaeoApp(sys.argv)

    file = QFile(":/qss/linnaeo.qss")
    file.open(QFile.ReadOnly)
    style = str(file.readAll())
    app.setStyleSheet(style)
    window = Linnaeo()
    app._window = window
    window.show()


    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
