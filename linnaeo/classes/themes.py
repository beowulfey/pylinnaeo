from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor


class AbstractTheme:
    pos = QColor(Qt.white)
    neg = QColor(Qt.white)
    pol = QColor(Qt.white)
    aro = QColor(Qt.white)
    phb = QColor(Qt.white)
    cys = QColor(Qt.white)
    gly = QColor(Qt.white)
    pro = QColor(Qt.white)

    def __init__(self):
        self.theme = {}
        self.initTheme()

    def initTheme(self):
        self.theme = {
            # Charged; positive
            "R": self.pos, "K": self.pos, "H": self.aro,
            # Charged, negative
            "D": self.neg, "E": self.neg,
            # Misc
            "C": self.cys, "G": self.gly, 'P': self.pro,
            # Aromatic
            "W": self.aro, "F": self.aro, "Y": self.aro,
            # Hydrophobic
            "A": self.pro, 'I': self.phb, 'L': self.phb,
            'M': self.phb, 'V': self.phb,
            # Polar
            'S': self.pol, 'T': self.pol, 'N': self.pol,
            'Q': self.pol
        }


class FirstTheme(AbstractTheme):
    """ My first attempt. Not recommended. Kept for historical purposes """
    pos = QColor(100, 140, 255)
    neg = QColor(255, 70, 90)
    cys = QColor(255, 255, 85)
    aro = QColor(145, 255, 168)
    gly = QColor(255, 255, 0)


class Theme2(AbstractTheme):
    """ Based on that lyft tool, loosely clustalX ish """
    pos = QColor(196, 88, 90)  # red
    neg = QColor(209, 97, 210)  # magenta
    pol = QColor(165, 242, 139)  # green
    aro = QColor(121, 234, 202)  # cyan
    phb = QColor(108, 126, 223)  # blue
    cys = QColor(247, 245, 161)  # yellow
    gly = QColor(254, 239, 238)  # light tan
    pro = QColor(251, 215, 188)  # light orange


class PaleByType(AbstractTheme):
    """ Paler version of Theme2. My favorite -- default colors."""
    pos = QColor(219, 138, 139) # red
    neg = QColor(225, 144, 226)  # magenta
    pol = QColor(190, 241, 172) # green
    aro = QColor(160, 237, 216)  # cyan
    phb = QColor(151, 164, 232) # blue
    cys = QColor(244, 242, 186) # yellow
    gly = QColor(247, 237, 236) # light tan
    pro = QColor(246, 222, 204) # light orange


class Mono(AbstractTheme):
    """ Crappy blue to red mono theme. Total junk"""
    phb = QColor(174, 98, 204)
    gly = QColor(157, 106, 216)
    pro = QColor(137, 116, 227)
    aro = QColor(127, 133, 236)
    pos = QColor(142, 167, 243)
    cys = QColor(161, 195, 248)
    neg = QColor(210, 234, 254)
    pol = QColor(183, 217, 252)