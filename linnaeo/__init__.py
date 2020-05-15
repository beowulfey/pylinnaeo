import logging
import sys

from PyQt5.QtCore import QFile
from PyQt5.QtWidgets import QApplication

from linnaeo.main import Linnaeo
from linnaeo.classes.views import LinnaeoApp

import time

start_time = time.clock()

def main():
    #start_time = time.clock()
    logging.basicConfig(level=logging.DEBUG)  # , format="%(asctime)s:%(levelname)s:%(message)s")
    appLogger = logging.getLogger("INIT")

    app = LinnaeoApp(sys.argv)
    appLogger.debug("Trying stylesheet")
    #file = QFile(":/qss/linnaeo.qss")
    #file.open(QFile.ReadOnly)
    #style = str(file.readAll())
    #app.setStyleSheet(style)
    appLogger.debug("Creating Window")

    window = Linnaeo()
    app._window = window
    appLogger.debug("Window created; displaying")
    window.show()
    appLogger.debug("Startup took %f seconds" % float(time.clock()-start_time))


    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
