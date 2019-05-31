#Mahakaran Sandhu, 2019
#This program uses a lot of packages.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy.polynomial.polynomial as polynomial
from copy import deepcopy
from scipy.optimize import curve_fit
from scipy.stats import sem
from scipy.stats.distributions import t
import fpdf
import os

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
def linear_fit(time_axis, abs_data, well_names,savedir, drange = 'None', yrange=None):
    """This function fits a linear curve to set of data points. Time_axis are the x-axis time data; abs_data are
    the y-axis absorbance data. The optional parameter drange can be used to specify the x-axis index range over 
    which to perform the fitting (Note:works by index, not by time values.)
    (list, list, tuple) --> np.array object"""
    if drange == 'None':
        linfit = np.polyfit(time_axis,abs_data , 1)
    elif drange != 'None':
        linfit = np.polyfit(time_axis[drange[0]:drange[1]],abs_data[drange[0]:drange[1]] , 1)
    
    func = np.poly1d(linfit)  
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    plt.plot(time_axis, abs_data, 'bo', label=well_names, markersize=2)
    plt.plot(time_axis,func(time_axis), 'b-',label="Linear fit",color='k' )
    plt.xlabel('time(s)')
    plt.ylabel('Abs(340nm)')
    plt.legend(loc='best')
    plt.text(250,0.250, 'slope='+str(linfit[0]))

    if yrange!=None:
        ax.set_ylim(bottom=yrange[0], top=yrange[1])
    
    plt.savefig(str(well_names)+'.png', dpi=100)

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
def get_MM_data_points(excel_file, savedir, dranges=None, omitRepl=False, minute=False, yrange=None):
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
            well_fit = linear_fit(timeDataSeconds, sorted_data_tuples[i][1], sorted_data_tuples[i][0], yrange)
            print (well_fit[0])
        else:
            if len(dranges)==len(sorted_data_tuples):
                
                print (data_tuples[i][0])
                print (i)
                well_fit = linear_fit(timeDataSeconds, sorted_data_tuples[i][1],sorted_data_tuples[i][0], savedir, (dranges[i][0], dranges[i][1]), yrange)
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

    
    


    #Apply Beer Lambert Law (also reflection about x-axis). Here we convert absorbance data into concentration data. 
    BeerLambertTransformed = []
    
    for i in triplicate_slopes:
        transform = []
        for j in i:
            transform.append(-j*60000)
        BeerLambertTransformed.append(transform)
    
    print('THESE ARE BEER-LAMBERT TRANSFORMED SLOPES:')
    print(BeerLambertTransformed)


    outfile = open(excel_file[:-5]+'_SLOPES.csv', 'w+')

    for i in BeerLambertTransformed:
        for j in i:
            if i.index(j) != 2:
                outfile.write(str(j))
                outfile.write(',')
            else:
                outfile.write(str(j))
        outfile.write('\n')
    outfile.close()
    




    if omitRepl != False:
        edited_triplicate_slopes = triplicateEditor(BeerLambertTransformed, omitRepl)
    else:
        edited_triplicate_slopes = BeerLambertTransformed
                
    #Here we do some statistics, obtaining the mean and standard error of the mean for the triplicates. 
    meanOfTriplicates = [np.mean(i) for i in edited_triplicate_slopes]
    semOfTriplicates = [sem(i) for i in edited_triplicate_slopes]    

    
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

    #Now compute confidence intervals on kinetic parameters
    alpha = 0.05
    num_data_pts = len(MMDataPoints[0])
    n_params = 2

    dof = max(0,num_data_pts-n_params)

    tval = t.ppf(1.0-alpha/2, dof)






   

    alpha = 0.05 # 95% confidence interval

    n = len(MMDataPoints[0])    # number of data points
    p = len(params) # number of parameters

    dof = max(0, n-p) # number of degrees of freedom

    tval = t.ppf(1.0 - alpha / 2.0, dof) # student-t value for the dof and confidence level

    for i, p,var in zip(range(n), params, np.diag(covar)):
        sigma = var**0.5
        print ('c{0}: {1} [{2}  {3}]'.format(i, p,
                                      p - sigma*tval,
                                      p + sigma*tval))






    K_m = params[1]
    V_max = params[0]


    print ('The Km is: ' + str(K_m))
    print ('The Vmax is: ' + str(V_max))


    print ('THIS IS PARAMS:\n')
    print (str(params))



    return (params)

#Unfortunately, we can't just plot an equation in Python...We need to generate data points etc. This function does
#all of this:
def plotMichaelisMenten(params, savedir, subsConcentrations, MMDataPoints, excel_file, MM_yrange): 
    """Given MM parameters (tuple(Vmax, Km))(generated by fitMichaelisMenten(), substrate concetrations (i.e. x-axis),
    and MM data points, plots a nice-looking MM plot."""
    substrate_fit = np.linspace(0.15,0,1000)
    velocities = [MMEquation(i, params[0], params[1]) for i in substrate_fit]
    uM_concs = [i*1000 for i in subsConcentrations]
    plt.plot(substrate_fit, velocities, color='k')
    plt.errorbar(subsConcentrations, MMDataPoints[0] ,xerr=None,yerr=MMDataPoints[1],fmt="none", capsize=5, ecolor='k', )
    plt.scatter(subsConcentrations, MMDataPoints[0], color='k')
    plt.ylabel('V (uM/min)')
    plt.xlabel('[S] (mM)')
    if MM_yrange != None:
        plt.ylim(MM_yrange)
    plt.title(excel_file[:-5])
    plt.savefig('MM_plot.png', dpi=300)
    plt.savefig('MM_plot.svg')



#Alright! Let's put all of the above functions together into 1 mega-function that you can call with a simple
#1-liner. This is probably quite amenable to developing into a command-line script, and I will probably pursue
#this direction sometime in the future. If this approach is taken, a control file (.ctl) will probably be required
#to specify all of the parameters.

def do_Michaelis_Menten(excel_file, savedir, drange, subsConcentrations,  omitConc, omitRepl=False, minute=False, yrange=None, MM_yrange=None):
    """Combines all of the above functions to generate a user-friendly 1-liner function. Also has the option
    omitConc which allows the user to omit problematic data points in the MM (this is concentration values, 
    not indexes). OmitRepl allows the user to omit specific replicates at a particular concentration (concentration
    index, not values here). There is redundancy between these 2 functionalities."""
    funcSubConcs = deepcopy(subsConcentrations)
    MM_data_points = get_MM_data_points(excel_file, savedir, drange, omitRepl, minute, yrange)
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
    plotMichaelisMenten(MMparams, savedir, funcSubConcs, MM_data_points, excel_file, MM_yrange)
    return (MM_data_points, MMparams)

#I think this is working. 

#This last function is added functionality for writing the data analysed above into a PDF. It's pretty expensive to run so be mindful. Pretty slow because of all the PNG files. Reducing DPI of PNG files likely to increase speed. 

def MM_makePDF(MMgraphImName,KM, Vmax):


    cwd = os.getcwd()
    current_experiment = cwd.split('/')[-1]
    current_experiment
    kineticReads = [i for i in os.listdir('./') if '_' not in i]
    kineticReads.sort()
    finalMMplot  = [i for i in os.listdir('./') if ('_' in i) and ('png' in i)]

    splitKRead = [kineticReads[x:x+12] for x in range(0, len(kineticReads), 12)]
    splitKRead


    betterSorted = []
    for i in splitKRead:
        out = []
        to_append= []
        for j in i:
            if ('10' in j) or ('11' in j) or ('12' in j):
                to_append.append(j)
            else:
                out.append(j)
        output = out+to_append
        betterSorted.append(output)
        
        
    pdf = fpdf.FPDF()


    pdf.add_page()
    pdf.set_font('Times','', 14)

    pdf.cell(100,100, 'MICHAELIS-MENTEN DATA ANALYSIS FOR EXPERIMENT '+current_experiment, border=0, align='L')
    #pdf.cell(100,10,betterSorted[0][1],border=0, align='L')
    for i in range(len(betterSorted)):
        pdf.add_page()
        pdf.image(betterSorted[i][0], 20,0,75)
        pdf.image(betterSorted[i][1], 110,0,75)
        pdf.image(betterSorted[i][2], 20,49,75)
        pdf.image(betterSorted[i][3], 110,49,75)
        pdf.image(betterSorted[i][4], 20,98,75)
        pdf.image(betterSorted[i][5], 110,98,75)
        pdf.image(betterSorted[i][6], 20,147,75)
        pdf.image(betterSorted[i][7], 110,147,75)
        pdf.image(betterSorted[i][8], 20,196,75)
        pdf.image(betterSorted[i][9], 110,196,75)
        pdf.image(betterSorted[i][10], 20,245,75)
        pdf.image(betterSorted[i][11], 110,245,75)

    pdf.add_page()
    pdf.image(MMgraphImName,30,100,140)
    pdf.cell(100,10, 'MICHAELIS-MENTEN PLOT FOR EXPERIMENT ' + current_experiment, border=0, align='L')
    pdf.ln(h=5)
    pdf.set_font('Times','', 10)

    pdf.cell(100,10, 'The KM is '+ KM, border=0, align='L')
    pdf.ln(h=5)
    pdf.cell(100,10, 'The Vmax is '+Vmax, border=0, align='L')
    
    pdf.output(current_experiment+'.pdf', 'F')
    
