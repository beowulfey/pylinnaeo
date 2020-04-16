#!/usr/bin/python3

from PyQt5.QtWidgets import QApplication
import sys

import bioglot

# Main function for running ProAlign/Bioglot/Whatever this is called.


def main():
    app = QApplication(sys.argv)
    window = bioglot.BioGlot()
    window.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
