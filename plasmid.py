#! /usr/bin/python3

import math

from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *

window = 0                                             # glut window number
width, height = 500, 400                               # window size
aspect = width/height


def draw_map(x, y, radius):
    ntris = 1000
    twicepi = math.pi*2.0
    #glEnable(GL_LINE_SMOOTH)
    glLineWidth(5.0)
    glBegin(GL_LINE_LOOP)
    for i in range(ntris):
        #glLineWidth(2.0)
        glVertex2f(x + radius * math.cos(i * twicepi / ntris), \
                   y + radius * math.sin(i * twicepi / ntris))
    glEnd()

def refresh2d(width, height):
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    adjwidth = height * aspect
    left = (width - adjwidth) / 2
#    print(left)
    glViewport(0, 0, width, height)
    glOrtho(0.0, width, 0.0, height, 0.0, 1.0)
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()

def draw():                                            # ondraw is called all the time
    glClearColor(1.0, 1.0, 1.0, 1.0)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) # clear the screen
    glLoadIdentity()                                   # reset position
    refresh2d(width, height)

    # draw Rectangle
    glColor3f(0.0, 0.0, 0.0)                           # blue?
    draw_map(250, 200, 150)


    glutSwapBuffers()                                  # important for double buffering
# initialization
glutInit()                                             # initialize glut
glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH)
glutInitWindowSize(width, height)                      # set window size
glutInitWindowPosition(0, 0)                           # set window position
window = glutCreateWindow("OMG PLASMID")              # create window with title
glutDisplayFunc(draw)                                  # set draw function callback
glutIdleFunc(draw)                                     # draw all the time
glutMainLoop()
