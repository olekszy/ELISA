import sys
import pandas as pd
import glob,os
import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from scipy.interpolate import interp1d
from pick import pick
from tabulate import tabulate

platetemplate = sys.argv[1] #Read template of plate
csv = sys.argv[2] #Read csv file from reader

#dir = os.chdir(workspace)

#wells = pd.read_csv("../Plateorder.csv", sep = ";", index_col=0) #read wells 

template = pd.read_csv(platetemplate, sep = ";", index_col=0) # Read sample sheet

def import1d(x):
    df = pd.read_csv(x, sep = ';', index_col=0)
    out = np.empty(df.shape[0], dtype=object)
    out[:] = df.values.tolist()
    return out
def importdata(x):
    df = pd.read_csv(x, sep = ';', index_col=0)
    df = df.apply(lambda x: x.str.replace(',','.'))
    return df.astype("float")
def calcc(y):
    a = []
    const = (max(y)-min(y))/2
    for i in y :
        z = np.absolute(i-const)
        a.append(z)
    return min(a)
#Develop 4pl curve

def calcC(y):
    return np.median(y)

def logistic4(x, A, B, C, D):
    """4PL lgoistic equation."""
    return ((A-D)/(1.0+((C/x)**B)))+D # Logistic equation

def residuals(p, y, x): # Deviations of Data
    """Deviations of data from fitted 4PL curve"""
    A,B,C,D = p
    err = y-logistic4(x, A, B, C, D) 
    return err

def peval(x, p):
    """Evaluated value at x with current parameters."""
    A,B,C,D = p
    return logistic4(x, A, B, C, D)

def concentration(y,D,B,C,A): #exchanged D and A paramters in scipy
    x = C*(((A-D)/(y-D))-1)**(1/B)
    return x  

def cuttable(x, wave):
    z = pd.read_csv(x, sep=";",header = None).fillna("")
    row = z[z[1].str.contains(wave)].index.values.astype(int)[0]
    name = z.iloc[row,1]
    table = z.iloc[row+2:row+10]
    table = pd.DataFrame(table).set_index([0])
    table = table.drop(columns=table.columns[(table == '').any()])
    table = table.apply(lambda x: x.str.replace(',','.'))
    return table.astype("float")

def r2value(y, y_predicted):
    ss_res = np.sum((y - y_predicted) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - (ss_res / ss_tot)
    return r2
def rsqtable(x, y):
    pd.DataFrame(columns = ["x", "y", "y - yavg", "(y-yavg)^2", "ypred", "ypred - yavg", "(ypred-yavg)^2"])
    rsq["x"] = x
    rsq["y"] = y
    yavg = y.mean()
    rsq["y - yavg"] = rsq["y"] - yavg
    rsq["(y-yavg)^2"] = rsq["y - yavg"]**2
    rsq["ypred"] = y_true
    rsq["ypred - yavg"] = rsq["ypred"] - yavg
    rsq["(ypred-yavg)^2"] = rsq["ypred - yavg"]**2
def readsample(temp,df,S):
    z = []
    #Read Standards 
    for column in temp.columns: # search  in columns
        for row in temp.index: #search in rows
            y = temp.at[row,column] #find specific place
            if S in y:
                z.append(row+":"+column) #fill list with locations of samples
    value = [] #create list for extracted values
    for i in z: 
        x = i[0] # take rows
        y = int(i[2:]) # take column
        value.append(df.at[x,y]) #fill list with extracted values of samples
        avgvalue = np.mean(value) #count averages of repetitions
    return value, avgvalue.item()
def xsamples(temp,df,popt):
    results = pd.DataFrame({"Sample":[], "Concentration":[]}) #create new dataframe
    samplelist = np.unique(temp.iloc[:].values) #extract sample names
    for i in samplelist: #loop over sample names extracted before
        a,b = readsample(temp, df, i) #Read samples from template
        z = concentration(b, *popt.tolist()) #Calculate concentration from 4PL Curve
        temporary = pd.DataFrame({"Sample":[i], "Concentration":[z.real]}) #create temporary table for better appending
        #print(temporary) #Check
        results = results.append(temporary, ignore_index=True) #append results to table
        #print(z) #Check
    return results
def ODbyELISA(results):
    z = results[results['Sample'].str.contains("1Standard")].values # extract 1Standard to 1000
    const = z.flat[1] #extract value
    conc = results['Concentration'].tolist() #change column to list
    lista = [] #create list for append
    #print(const)
    for i in conc:    
        x = (1000*i)/const #formula over all rows
        lista.append(x)
    results["ELISA units"] = lista #create column from list
    return results
def omitstandards(standards,measured,diluted):
    title = 'Please choose Standards to omit'
    options = standards
    option, index = pick(options, title)
    selected = pick(options, title, multiselect=True, min_selection_count=0)
    f = standards.index(selected[0][0])
    standards.remove(selected[0][0])
    y = np.delete(measured,f)
    x = np.delete(diluted,f)
    print ("Omiting " + str(selected[0][0]))
    return f, y, x
##################################PROGRAM STARTS HERE#################################################
#csv = "IgA_29-05-2020__prot_S1_spec_plytka6_B.csv" #CSV from Reader
def analysis(csv,platetemplate):   
    print(csv)
    df450 = cuttable(csv, "(450)") #Cut 450 nm table
    df570 = cuttable(csv, "(570)") #Cut 570 nm table
    print("Properly imported plates")
    #Background substraction 
    dfavg  = df450.sub(df570) #Cut background
    print("Plates substracted")
    #Read Standards 
    standards = [] # Create standards
    countstandards = int(input("How many standards did you use")) #input number of standards

    for i in range(1,countstandards+1):
        x = str(i)+"Standard"
        standards.append(x)# append list of S1..S2.. etc
    np.unique(standards.append("B"))

    ##### Define Concentrations
    value = input("Starting Concentration") #Input number of standards
    #Input Dilution factor
    dil = float(input("Dilution factor")) #Change to float

    list = [] # create list
    list.append(value) #Add starting value

    for i in range(countstandards-1):
        value = float(value)/dil #Create series for standarization
        list.append(value) # append series values
    list.append(0)
    a = np.asarray(list) # save as array
    a = np.loadtxt(a, dtype='float') #delete dtype at the end
    x = a
    print("Standards created")
    col = []

    for i in standards:   
        z,avg = readsample(template, dfavg, i)
        col.append(avg)
    col = np.asarray(col)
    #Attach Standards OD
    y = col

    x=a #crucial!!!!!!!
    print("Standards attached")
    print(x)

    null, y, x = omitstandards(standards,col,a)

    from scipy.optimize import curve_fit
    popt, pcov = curve_fit(logistic4, x, y)
    y_true = logistic4(x, *popt)

    f = interp1d(x, y_true)
    f2 = interp1d(x, y_true, kind='cubic')

    xnew = np.linspace(min(x), max(x))

    print("ABCD paramteres completed")

    plt.plot(x, y, 'o', xnew, f(xnew), '-', xnew, f2(xnew), '--')
    plt.title(csv+" "+"4PL curve", fontdict=None, loc='center', pad=None)
    plt.legend(['data', 'linear', 'cubic'], loc='best')
    plt.xlabel("Concentration")
    plt.ylabel("OD")
    print("R2 value" +" "+ str(r2value(y, y_true)))

    print("Saving results to files")
    filename = csv.split("/")[-1]
    print(filename[0:-4])

    plt.savefig(filename[0:-4])

    results = xsamples(template,dfavg,popt)
    print(tabulate(results, headers='keys', tablefmt='psql'))
    final = ODbyELISA(results)

    final.to_csv(csv[0:-4]+"_final.csv", sep = ";", header = True)

analysis(csv,platetemplate)