#This program uses a lot of packages.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy.polynomial.polynomial as polynomial
from copy import deepcopy
from scipy.optimize import curve_fit
from scipy.stats import sem

#Function Definitions

#This function is needed to convert time data from Excel files into decimals. 
def convert2seconds(inputstr, minute=False):
    """Converts time in format HH:MM:SS to decimal point seconds. Choice of minutes or seconds (default).
    str --> float(min) OR int(s)"""
    (h, m, s) = inputstr.split(':')
    if minute==False:
        result = int(h) * 3600 + int(m) * 60 + int(s)
    elif minute==True:
        result = float(h) * 60 + float(m) + (float(s)/60)
    return result

#This function fits a linear curve to set of data points. It's needed to find the initial rate of enzymatic activity.
def linear_fit(time_axis, data_axis, drange = 'None', yrange=None):
    """This function fits a linear curve to set of data points. Time_axis are the x-axis time data; data_axis are
    the y-axis absorbance data. The optional parameter drange can be used to specify the x-axis index range over 
    which to perform the fitting (Note:works by index, not by time values.)
    (list, list, tuple) --> np.array object"""
    if drange == 'None':
        linfit = np.polyfit(time_axis,data_axis , 1)
    elif drange != 'None':
        linfit = np.polyfit(time_axis[drange[0]:drange[1]],data_axis[drange[0]:drange[1]] , 1)
    
    func = np.poly1d(linfit)  
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    plt.plot(time_axis, data_axis, 'bo', label="Data")
    plt.plot(time_axis,func(time_axis), 'b-',label="Polyfit")
    if yrange!=None:
        ax.set_ylim(bottom=yrange[0], top=yrange[1])
    
    plt.show()
    return linfit

#A simple function to get concentration given A, epsilon, and pathlength.
def BeerLambertTransform(absorbance, epsilon, pathlength):
    """Given an absorbance, an epsilon, and a pathlength, returns the concentration
    float --> float"""
    concentration = absorbance/(epsilon*pathlength)
    return concentration



def triplicateEditor(triplicatesList, omitRepl):
    """This function takes in a list triplicates (list of lists), and a list of 2-tuples, where the first element
    of the tuple is the index (important, INDEX, not the concetration itself) of the triplicate at concentration x
    that is to be edited, and where the second element of the tuple is a list containing the indexes of the 
    replicates in that index(concetration x) that are to be omitted. Returns a list with the same format as the input
    list triplicatesList.
    (list, list) --> list
    (list of lists, lists of tuples ([(concindex, list(replicate indexes to be removed))])) --> list of lists"""
    
    edit_concs = [i[0] for i in omitRepl]
    edit_repls = [i[1] for i in omitRepl]
    omission_triplicate_slopes = []

    for triplicate in triplicatesList:
        if triplicatesList.index(triplicate) in edit_concs:
            edited_triplicate = []
            for j in edit_concs:
                if j == triplicatesList.index(triplicate):
                    index_of_editconc_hit = edit_concs.index(j)

            replicates_to_omit = edit_repls[index_of_editconc_hit]

            for replicate in triplicate:
                if triplicate.index(replicate) not in replicates_to_omit:
                    print (replicate)
                    edited_triplicate.append(replicate)
            omission_triplicate_slopes.append(edited_triplicate)
        else:
            omission_triplicate_slopes.append(triplicate)

    print ('This is triplicate slopes:')
    print (triplicatesList)
    print ('This is omission slopes:')
    print (omission_triplicate_slopes)
    return omission_triplicate_slopes

#This is the big function, a lot is happening here. 
def get_MM_data_points(excel_file, dranges=None, omitRepl=False, minute=False, yrange=None):
    """This function reads an Excel file, orders the data according to
    well, invokes both convert2seconds() and linear_fit(), performs initial statistics on slopes by replicate, and 
    returns a tuple of lists with the data needed.
    (inputfile(.xlsx), list) --> tuple of lists (list, list) """
    #Note: omitRepl should be a list of tuples, where tuple[0] defines the concentration index (in range(12))and tuple[1]
    # defines the replicate index (0,1, or 2)
    
    #read Excel data
    df=pd.read_excel(excel_file)         
   

    #here we extract the data in each well (absorbance as a function of time). We add the data in each well to 
    #data_tuples as a tuple, where tuple[0] gives the string (well name), and tuple[1] is a list with the
    #(absorbance as a function of time).    
    data_tuples = []
    columns = list(df)
    well_names = columns[1:]
    
    #range(37) is used because I am assuming 12 concetration values in triplicate. Change if different format
    #is used. 
    well_data = [df.iloc[:,i].tolist() for i in range(37) if type(df.iloc[:,i].tolist()[0]) == float] 
    for i in range(len(well_names)):
        data_tuples.append((well_names[i], well_data[i]))

    #Because the Excel data is not sorted according to the rows in a 96-well plate, we now perform this sorting.
    #In this way, data from wells will appear sequentially by 96-well plate row and column, e.g. D1,D2,...,D12,E1,...,E12 etc. 
    #The data will be a list of tuples, where the first index of the tuple (tuple[0]) is the well name (type==str), e.g. D1, and the 
    #second index of the tuple (tuple[1]) is a list with the absorbance data. 
    #NOTE: len(data_tuples) should = len(range_list) = number of wells in experiment. 

    sorted_data_tuples = []
    alphabet = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    for row in alphabet:
        for well in data_tuples:
            if well[0][0]==row:
                sorted_data_tuples.append(well)

    #Here we extract time data. This will give us the x-axis data for the MM plot. 
    time_data = df.iloc[:,0].tolist()
    timeDataSeconds = [convert2seconds(t.strftime('%H:%M:%S'), minute) for t in time_data]
    
     
    #Here, we iterate over the sorted well data. If no range for linear_fit is specified, the all the data is used;
    #else, the only the range-specified absorbance data is used. The output of the linear_fit is appended to the list fit_data. 
    fit_data = []
    
    for i in range(len(sorted_data_tuples)):
        if dranges == None:
            print (sorted_data_tuples[i][0])
            print (i)
            well_fit = linear_fit(timeDataSeconds, sorted_data_tuples[i][1], yrange)
            print (well_fit[0])
        else:
            if len(dranges)==len(sorted_data_tuples):
                
                print (data_tuples[i][0])
                print (i)
                well_fit = linear_fit(timeDataSeconds, sorted_data_tuples[i][1], (dranges[i][0], dranges[i][1]), yrange)
                print (well_fit[0])
            else:
                print ('ERROR: LENGTH OF SUPPLIED RANGES != LENGTH OF DATA ')

        fit_data.append(well_fit)

    #we now extract only the slopes from the fit_data elements.This makes things a little easier.  
    slopes_only = []
    for i in fit_data:
        slopes_only.append(i[0])

    #here we take every 12th element in slopes_only and combine them together. This gives a list of lists wherein
    #each element is contains data from triplicates at the same concentration. Note that this assumes D1, E1 and F1
    #are replicates of each other, just as D5, E6 and F6. Ensure that if you want to use this program, you set up
    #your 96 well plates in this format. 
    
    
    #need to make this robust to changes in order of columns in excel file. Doing this will reduce the chance of error in the 
    #MM data analysis, and will also bypass human error (hopefully). 
    triplicate_slopes = []
    for i in range(12):
        repl1 = slopes_only[i]
        repl2 = slopes_only[i+12]
        repl3 = slopes_only[i+12+12]
        outlist = [repl1, repl2, repl3]
        triplicate_slopes.append(outlist)

    
    if omitRepl != False:
        edited_triplicate_slopes = triplicateEditor(triplicate_slopes, omitRepl)
    else:
        edited_triplicate_slopes = triplicate_slopes
    

            
                        
    #Apply Beer Lambert Law (also reflection about x-axis). Here we convert absorbance data into concentration data. 
    BeerLambertTransformed = []
    
    for i in edited_triplicate_slopes:
        transform = []
        for j in i:
            transform.append(-j*60000)
        BeerLambertTransformed.append(transform)
    
                
                
    #Here we do some statistics, obtaining the mean and standard error of the mean for the triplicates. 
    meanOfTriplicates = [np.mean(i) for i in BeerLambertTransformed]
    semOfTriplicates = [sem(i) for i in BeerLambertTransformed]    

    
    #Normalize to 0
    #we don't touch the SEM because normalisation doesn't affect this. 
    normalized_meansOfTriplicates = [(i-meanOfTriplicates[-1]) for i in meanOfTriplicates]
   

    return (normalized_meansOfTriplicates, semOfTriplicates)


#We are getting to the business end of things, time to define the MM equation:
def MMEquation(SubstrateConc, VMax, K_m): 
    """The MM equation."""
    velocity = (VMax*SubstrateConc)/(K_m + SubstrateConc)
    return velocity

#We need to fit the data we obtained in get_MM_data_points to the MM equation. The function below does this. 
def fitMichaelisMenten(subsConcentrations, MMDataPoints): 
    """Fits data to the Michaelis Menten equation. Returns the parameters Km and Vmax. 
    (substrate_Concentrations, MMDataPoints(a tuple of format (points, stdevs))) --> (Vmax, Km)
    (list, 2-tuple) --> 2-tuple"""
    params, covar = curve_fit(MMEquation, subsConcentrations,MMDataPoints[0], p0=[2, 2]) 
    K_m = params[1]
    V_max = params[0]
    print ('The Km is: ' + str(K_m))
    print ('The Vmax is: ' + str(V_max))
    return params

#Unfortunately, we can't just plot an equation in Python...We need to generate data points etc. This function does
#all of this:
def plotMichaelisMenten(params, subsConcentrations, MMDataPoints, excel_file): 
    """Given MM parameters (tuple(Vmax, Km))(generated by fitMichaelisMenten(), substrate concetrations (i.e. x-axis),
    and MM data points, plots a nice-looking MM plot."""
    substrate_fit = np.linspace(0.15,0,1000)
    velocities = [MMEquation(i, params[0], params[1]) for i in substrate_fit]
    uM_concs = [i*1000 for i in subsConcentrations]
    plt.plot(substrate_fit, velocities, color='0.4')
    plt.errorbar(subsConcentrations, MMDataPoints[0] ,xerr=None,yerr=MMDataPoints[1],fmt="none", capsize=5, ecolor='0.4', )
    plt.scatter(subsConcentrations, MMDataPoints[0])
    plt.ylabel('V (uM/min)')
    plt.xlabel('[S] (mM)')
    plt.title(excel_file[:-5])
    plt.savefig(str(excel_file)+'.png', dpi=100)



#Alright! Let's put all of the above functions together into 1 mega-function that you can call with a simple
#1-liner. This is probably quite amenable to developing into a command-line script, and I will probably pursue
#this direction sometime in the future. If this approach is taken, a control file (.ctl) will probably be required
#to specify all of the parameters.

def do_Michaelis_Menten(excel_file, drange, subsConcentrations,  omitConc, omitRepl=False, minute=False, yrange=None):
    """Combines all of the above functions to generate a user-friendly 1-liner function. Also has the option
    omitConc which allows the user to omit problematic data points in the MM (this is concentration values, 
    not indexes). OmitRepl allows the user to omit specific replicates at a particular concentration (concentration
    index, not values here). There is redundancy between these 2 functionalities."""
    funcSubConcs = deepcopy(subsConcentrations)
    MM_data_points = get_MM_data_points(excel_file, drange, omitRepl, minute, yrange)
    print ("And behold, the Michaelis-Menten Curve:")
    if len(omitConc)==0:
        MMparams = fitMichaelisMenten(subsConcentrations, MM_data_points)
    else:
        for i in omitConc:
            ind = funcSubConcs.index(i)
            del funcSubConcs[ind]
            del MM_data_points[0][ind]
            del MM_data_points[1][ind]
  
        MMparams = fitMichaelisMenten(funcSubConcs, MM_data_points)
    plotMichaelisMenten(MMparams, funcSubConcs, MM_data_points, excel_file)
    return (MM_data_points, MMparams)

#I think this is working. 