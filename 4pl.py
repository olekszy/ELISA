import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from tkinter import Tk
from tkinter.filedialog import askopenfilename
import pandas as pd
import glob,os
from tkinter import *
from tkinter import ttk
workspace = "/home/fagi/Pulpit/ELISA"
dir = os.chdir(workspace)
Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file
df = pd.read_csv(filename, sep = ';', index_col=0)
df = df.apply(lambda x: x.str.replace(',','.'))
plate = pd.read_csv("Plateorder.csv", sep = ";", index_col=0)
linear_df=pd.DataFrame({'well':plate.values.ravel(), 'OD':df.values.ravel()})

def logistic4(x, A, B, C, D):
    """4PL lgoistic equation."""
    return ((A-D)/(1.0+((x/C)**B))) + D

def residuals(p, y, x):
    """Deviations of data from fitted 4PL curve"""
    A,B,C,D = p
    err = y-logistic4(x, A, B, C, D)
    return err

def peval(x, p):
    """Evaluated value at x with current parameters."""
    A,B,C,D = p
    return logistic4(x, A, B, C, D)

top = Tk()
standards = Listbox(top)
standards.insert(END, linear_df[linear_df.columns[0]])
standards.pack()
top.mainloop()