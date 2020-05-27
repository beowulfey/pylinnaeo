import logging
import sys
import time

from linnaeo.main import LinnaeoApp, Linnaeo

start_time = time.perf_counter()


def main():
    logging.basicConfig(level=logging.DEBUG)
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
    window.show()
    appLogger.info("Startup took %f seconds" % float(time.perf_counter()-start_time))


    sys.exit(app.exec_())


if __name__ == '__main__':
    main()

