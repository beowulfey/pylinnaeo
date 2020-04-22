#!/usr/bin/python3

from PyQt5.QtWidgets import QApplication
import sys

import sherlock

# Main function for running ProAlign/Bioglot/Whatever this is called.


def main():
    app = QApplication(sys.argv)
    window = sherlock.Sherlock()
    window.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
