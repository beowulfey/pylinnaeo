__version__ = 'v0.4.0'
import logging
import sys
import time

from PyQt5.QtCore import QFile, QByteArray, Qt
from PyQt5.QtGui import QIcon, QPalette, QColor
from PyQt5.QtWidgets import QStyle

from linnaeo.main import LinnaeoApp, Linnaeo

start_time = time.perf_counter()
conout = None


def main():
    logging.basicConfig(level=logging.DEBUG, stream=sys.stderr)
    appLogger = logging.getLogger("INIT")
    #handler = logging.StreamHandler(sys.stderr)
    #appLogger.addHandler(handler)

    app = LinnaeoApp(sys.argv)
    app.setApplicationName("Linnaeo")
    app.setApplicationVersion(__version__)

    #app.setStyle("Fusion")
    #app.setMode("dark")

    if app.mode == "dark":
        dark_palette = QPalette()
        dark_palette.setColor(QPalette.Window, QColor(53, 53, 53))
        dark_palette.setColor(QPalette.WindowText, Qt.white)
        dark_palette.setColor(QPalette.Base, QColor(25, 25, 25))
        dark_palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
        dark_palette.setColor(QPalette.ToolTipBase, Qt.white)
        dark_palette.setColor(QPalette.ToolTipText, Qt.white)
        dark_palette.setColor(QPalette.Text, Qt.white)
        dark_palette.setColor(QPalette.Button, QColor(53, 53, 53))
        dark_palette.setColor(QPalette.ButtonText, Qt.white)
        dark_palette.setColor(QPalette.BrightText, Qt.red)
        dark_palette.setColor(QPalette.Link, QColor(42, 130, 218))
        dark_palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
        dark_palette.setColor(QPalette.HighlightedText, Qt.black)
        app.setPalette(dark_palette)
        app.setStyleSheet("QToolTip { font-family: Default-Noto; color: #000000; background-color: #efefef; border: 1px solid white; }")

    else:
        app.setStyleSheet("QToolTip { font-family: Default-Noto; }")
    appLogger.debug("Creating Window")

    window = Linnaeo()
    app.mainWindow = window
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
