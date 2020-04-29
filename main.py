#!/usr/bin/python3

from PyQt5.QtWidgets import QApplication
import sys
import logging
import sherlock


def main():
    logging.basicConfig(level=logging.DEBUG)#, format="%(asctime)s:%(levelname)s:%(message)s")
    app = QApplication(sys.argv)
    window = sherlock.Sherlock()
    window.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
