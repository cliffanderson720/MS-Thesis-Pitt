# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 11:52:10 2018

@author: cdhig
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import numpy as np

height = 4*2.54 #inches
width  = 10*2.54 # inches
x = np.linspace(0,width)
def line(axobject,x,m,thickness,alph):
    y = m*x/width*height
    axobject.plot(x,y,'k',linewidth=thickness,alpha=alph)
    

fig,ax = plt.subplots(1,figsize=(width/2.54,height/2.54))
for i in range(101):
    line(ax,x,i/100,0.1,0.5)
    if i % 5 == 0:
        line(ax,x,i/100,0.4,0.75)
    if i % 10 == 0:
        line(ax,x,i/100,0.6,1)


ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlim(xmin=7.5,xmax=width)
ax.set_ylim(ymin=0,ymax=height)
ax.grid(which='both',alpha=0.5,linewidth=0.2)
plt.subplots_adjust(left=0.03, right=0.99, top=0.99, bottom=0.05)
plt.savefig('HCT.png',dpi=600)
plt.show()

#%% 
'''
for original, see URL:
    https://stackoverflow.com/questions/29400116/using-matplotlib-how-can-i-print-something-actual-size
'''
left_margin = 1.   # cm
right_margin = 1.  # cm
figure_width = 25.94 # cm: 10"*2.54 = 27.94 minus left and right margins
figure_height = 10. # cm
top_margin = 1.    # cm
bottom_margin = 1. # cm

box_width = left_margin + figure_width + right_margin   # cm
box_height = top_margin + figure_height + bottom_margin # cm

cm2inch = 1/2.54 # inch per cm

# specifying the width and the height of the box in inches
fig,ax = plt.subplots(1,figsize=(box_width*cm2inch,box_height*cm2inch))
xmin = 8 # cm
x = np.linspace(xmin,figure_width)
def line(axobject,x,m,thickness,alph,xmin=xmin):
    y = m*x/figure_width*figure_height
    axobject.plot(x-xmin,y,'k',linewidth=thickness,alpha=alph)
    
for i in range(101):
    line(ax,x,i/100,0.1,0.5)
    if i % 5 == 0:
        line(ax,x,i/100,0.4,0.75)
    if i % 10 == 0:
        line(ax,x,i/100,0.6,1)

fig.subplots_adjust(left   = left_margin / box_width,
                    bottom = bottom_margin / box_height,
                    right  = 1. - right_margin / box_width,
                    top    = 1. - top_margin   / box_height,
                    )
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlim(xmin=0.0,xmax=figure_width-xmin)
ax.set_ylim(ymin=0.0,ymax=figure_height)
ax.grid(which='both',alpha=0.5,linewidth=0.2)
fig.savefig('HCT.png', dpi=600)