__version__ = 'v0.2.2'
import logging
import sys
import time

from PyQt5.QtGui import QIcon

from linnaeo.main import LinnaeoApp, Linnaeo

start_time = time.perf_counter()
conout = None


def main():
    logging.basicConfig(level=logging.DEBUG)
    appLogger = logging.getLogger("INIT")
    #handler = logging.StreamHandler(sys.stderr)
    #appLogger.addHandler(handler)

    app = LinnaeoApp(sys.argv)
    app.setApplicationName("Linnaeo")
    app.setApplicationVersion(__version__)
    #console = PythonConsole()
    #console.show()
    #console.eval_in_thread()
    appLogger.debug("Trying stylesheet")
    #file = QFile(":/qss/linnaeo.qss")
    #file.open(QFile.ReadOnly)
    #style = str(file.readAll())
    #app.setStyleSheet(style)
    appLogger.debug("Creating Window")

    window = Linnaeo()
    app._window = window
    window.setWindowTitle("Linnaeo")
    icon = QIcon(":/icons/linnaeo.ico")
    window.setWindowIcon(icon)
    window.show()
    #del icon
    appLogger.info("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    appLogger.info("\t\t          STARTUP COMPLETE!")
    appLogger.info(f"\t\tLoaded linnaeo-%s in {float(time.perf_counter()-start_time):.3f} seconds" % __version__)
    appLogger.info("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()

