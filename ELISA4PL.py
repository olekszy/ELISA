import functionselisa as fel
import sys
import pandas as pd
import glob,os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from scipy.interpolate import interp1d
from pick import pick
from tabulate import tabulate
from os import listdir
from os.path import isfile, join

workdir = "workfolder"
sheets = "samplesheets"
samples = [f for f in listdir(workdir) if isfile(join(workdir, f))] #find folders for samples

templates = [f for f in listdir(sheets) if isfile(join(sheets, f))] #find folders for templates
title = 'Please choose files to analyse'

chosensamples = pick(samples, title, multiselect=True, min_selection_count=0)

print (chosensamples)

samplelist = []
for samples in chosensamples:
    title = "Attach template to " + str(samples[0])
    sheet, index = pick(templates, title)
    print(sheet +" attached to " + str(samples[0]))
    samplelist.append(samples[0]+":"+sheet)

finalresults = pd.DataFrame()
standards = [] # Create standards
#Standards input
while True:
    try:
        countstandards = int(input("How many standards did you use? "))
    except ValueError:
        print("Incorrect Input") # try again
        continue
    else:
        #Numeber was successfully typed
        #we're ready to exit the loop.
        break

##### Define Concentrations
while True:
    try:
        value = input("Starting Concentration ")
    except ValueError:
        print("Incorrect Input") # try again
        continue
    else:
        #Numeber was successfully typed
        #we're ready to exit the loop.
        break #Input number of standards
#Input Dilution factor
while True:
    try:
        dil = float(input("Dilution factor "))
    except ValueError:
        print("Incorrect Input") # try again
        continue
    else:
        #Numeber was successfully typed
        #we're ready to exit the loop.
        break #Input number of standards
     #Change to float

for sample in samplelist:
    csv,platetemplate = sample.split(":")
    csv = workdir+"/"+csv
    platetemplate = sheets+"/"+platetemplate
    fn = str(csv) 
    print("Start analysis " + csv + " with template " + platetemplate)
    final, r2 = fel.analysis(csv,platetemplate,countstandards,value,dil,fn)
    final["R^2"] = round(r2,4)
    print(tabulate(final, headers='keys', tablefmt='psql'))
    finalresults = finalresults.append(final)

name = input("How to save your analysis?")
FirstStandards = finalresults[finalresults['Sample'].str.match("1Standard")]
print("Printing 1 Standards")
StandardsMean = FirstStandards["Concentration"].mean()

FinalConcentrations = finalresults['Concentration'].tolist() #change column to list
listConc = [] #create list for append
for i in FinalConcentrations:    
    x = (1000*i)/StandardsMean #formula over all rows
    listConc.append(x)
finalresults["ProperELISA units"] = listConc #create column from list
#print(finalresults)
finalresults = finalresults[~finalresults['Sample'].isin(["B", "no antigen"])]
finalresults.to_csv("results/"+name+".csv", sep = "\t")