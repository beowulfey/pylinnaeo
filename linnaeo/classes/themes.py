from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor


class AbstractTheme:
    """
    Here are the approximate residue percentages based on category.
    Use this to derive the order from lightest to darkest.
    phb 22.9   ala 8.3
    pol 20.5   gly 7.9
    pos 12.7   pro 4.7
    neg 12     cys 1.5
    aro 9.2
    from Baud and Kaulin, 1999 (10.1073/pnas.96.22.12494)
    """

    phb = QColor(Qt.white)
    pol = QColor(Qt.white)
    pos = QColor(Qt.white)
    neg = QColor(Qt.white)
    aro = QColor(Qt.white)
    ala = QColor(Qt.white)
    gly = QColor(Qt.white)
    pro = QColor(Qt.white)
    cys = QColor(Qt.white)

    def __init__(self):
        self.theme = {}
        self.initTheme()

    def initTheme(self):
        """ My standard categorizing. Can be changed per theme. """
        self.theme = {
            # Hydrophobic
            'I': self.phb, 'L': self.phb, 'M': self.phb,
            'V': self.phb,
            # Polar
            'S': self.pol, 'T': self.pol, 'N': self.pol,
            'Q': self.pol,
            # Charged; positive
            "R": self.pos, "K": self.pos, "H": self.pos,
            # Charged, negative
            "D": self.neg, "E": self.neg,
            # Aromatic
            "W": self.aro, "F": self.aro, "Y": self.aro,
            # Misc
            "A": self.ala, "G": self.gly, 'P': self.pro,
            "C": self.cys,
        }

    """def getDescr(self):
        descrdict = {
            'W':'Aromatic', 'D':'Negative', 'R':'Positive', 'S':'Polar','I':'Hydrophobic'
        }
        tcolor = '#FFFFFF' if color.getHsl()[2] / 255 * 100 <= 50 else '#000000'
        char = '<span style=\"background-color: %s; color: %s\">%s</span>' % (color.name(), tcolor, char)
        test1 = ['I','S','R',]
        for 
        string = 'Color categories:\n"""


class Comments(AbstractTheme):
    pass

class ColorSafe(AbstractTheme):
    """ Colorblind-friendly. See it here https://bit.ly/2U7cSr4 """
    phb = QColor('#A1C8EC')  # blue 1
    pol = QColor('#C85200')  # orange 1
    pos = QColor('#FF800F')  # orange 2
    neg = QColor('#FBBD77')  # orange 3
    aro = QColor('#5F9ED1')  # blue 2
    ala = QColor('#0167AB')  # blue 3
    gly = QColor('#EAE9E9')  # gray 1
    pro = QColor('#CFCFCF')  # gray 2
    cys = QColor('#898989')  # gray 3




class Rainbow(AbstractTheme):
    """ All same value, so I don't like this one much. """
    # ff5122 (orange) #ffbf22 (light-orange) #d0ff22 (yellow)
    # 62ff22 (li-grn) #22ff51 (green)        #22ffc0 (li-blu)
    ##22cfff (medblu)# 2261ff (blue)  #c022ff (purp)
    orange = QColor('#ffbf22')
    red = QColor('#ff5122')
    mblu = QColor('#22cfff')
    blue = QColor('#2261ff')
    green = QColor('#d0ff22')
    cys = QColor('#c022ff')
    white = QColor('#FFFFFF')
    ala = QColor('#EAE9E9')

    def __init__(self):
        super().__init__()
        self.theme = {}
        self.initTheme()

    def initTheme(self):
        """ My standard categorizing. Can be changed per theme. """
        self.theme = {
            # Hydrophobic
            'I': self.green, 'L': self.green, 'M': self.green,
            'V': self.green,
            # Polar
            'S': self.orange, 'T': self.orange, 'N': self.mblu,
            'Q': self.mblu,
            # Charged; positive
            "R": self.red, "K": self.red, "H": self.red,
            # Charged, negative
            "D": self.mblu, "E": self.mblu,
            # Aromatic
            "W": self.blue, "F": self.blue, "Y": self.blue,
            # Misc
            "A": self.ala, "G": QColor('#EAE9E9'), 'P': self.orange,
            "C": self.cys,
        }

    """
    phb = QColor('#ff5122')
    pol = QColor('#62ff22')
    pos = QColor('#22ffc0')
    neg = QColor('#ff5122')
    aro = QColor('#2261ff')
    ala = QColor('#22ff51')
    gly = QColor('#22cfff')
    pro = QColor('#d0ff22')
    cys = QColor('#c022ff')
    
    phb = QColor('#ff5122')
    pol = QColor('#ffbf22')
    pos = QColor('#d0ff22')
    neg = QColor('#62ff22')
    aro = QColor('#22ff51')
    ala = QColor('#22ffc0')
    gly = QColor('#22cfff')
    pro = QColor('#2261ff')
    cys = QColor('#c022ff')
    """


class PaleByType(AbstractTheme):
    """ Paler version of Theme2. My favorite -- default colors."""
    phb = QColor('#97A4E8')
    pol = QColor('#BEF1AC')
    pos = QColor('#DB8A8B')
    neg = QColor('#E190E2')
    aro = QColor('#A0EDD8')
    ala = QColor('#F7EDEC')
    gly = QColor('#EAE9E9')
    pro = QColor('#F6DECC')
    cys = QColor('#F4F2BA')


class Bold(AbstractTheme):
    """ Based on that colorbox.io, loosely clustalX ish """
    phb = QColor('#6C7EDF')
    pol = QColor('#A5F28B')
    pos = QColor('#C4585A')
    neg = QColor('#D161D2')
    aro = QColor('#79EACA')
    ala = QColor('#FEEFEE')
    gly = QColor('#EAE9E9')
    pro = QColor('#FBD7BC')
    cys = QColor('#F7F5A1')


    """
    pos = QColor(196, 88, 90)  # red #C4585A
    neg = QColor(209, 97, 210)  # magenta #D161D2
    pol = QColor(165, 242, 139)  # green #A5F28B
    aro = QColor(121, 234, 202)  # cyan #79EACA
    phb = QColor(108, 126, 223)  # blue #6C7EDF
    cys = QColor(247, 245, 161)  # yellow #F7F5A1
    gly = QColor(254, 239, 238)  # light tan #FEEFEE
    pro = QColor(251, 215, 188)  # light orange #FBD7BC
    """

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

class Grayscale(AbstractTheme):
    phb = QColor('#D6D6D6')  # blue 1
    pol = QColor('#B3B3B3')  # orange 1
    pos = QColor('#8F8F8F')  # orange 2
    neg = QColor('#6B6B6B')  # orange 3
    aro = QColor('#000000')  # blue 2
    ala = QColor('#E8E8E8')  # blue 3
    gly = QColor('#FAFAFA')  # gray 1
    pro = QColor('#474747')  # gray 2
    cys = QColor('#242424')  # gray 3


    """
    phb = QColor('#FAFAFA')  # blue 1    D6D6D6
    pol = QColor('#E8E8E8')  # orange 1  B3B3B3
    pos = QColor('#D6D6D6')  # orange 2  8F8F8F
    neg = QColor('#B3B3B3')  # orange 3  6B6B6B
    aro = QColor('#8F8F8F')  # blue 2    000000
    ala = QColor('#6B6B6B')  # blue 3    E8E8E8
    gly = QColor('#474747')  # gray 1    FAFAFA
    pro = QColor('#242424')  # gray 2    474747
    cys = QColor('#000000')  # gray 3    242424"""

    """
    gly = QColor(255, 255, 255)
    phb = QColor(223, 223, 223)
    pro = QColor(191, 191, 191)
    aro = QColor(159, 159, 159)
    pol = QColor(128, 128, 128)
    cys = QColor(96, 96, 96)
    pos = QColor(64, 64, 64)
    neg = QColor(32, 32, 32)"""


class FirstTheme(AbstractTheme):
    """ My first attempt. Not recommended. Kept for historical purposes """
    pos = QColor(100, 140, 255)
    neg = QColor(255, 70, 90)
    cys = QColor(255, 255, 85)
    aro = QColor(145, 255, 168)
    gly = QColor(255, 255, 0)
