#Written by Mahakaran Sandhu, 2018-2019. 
#Honours Year at the RSC, ANU. 
import matplotlib.pyplot
def plot_xvg(xvg_file, x_name, y_name):
    """Reads a GROMACS xvg file containing data from gmx energy. Plots the resulting data. You will need to provide labels for
    axes. (infile, str, str) --> matplotlib object."""
    data = open(xvg_file)
    data = data.readlines()
    data = [i.strip() for i in data]
    data = [i for i in data if ('@' not in i) and ('#' not in i)]
    x_axis=[]
    y_axis = []
    for i in data: 
        r = i.split()
        x_axis.append(float(r[0]))
        y_axis.append(float(r[1]))
    matplotlib.pyplot.plot(x_axis, y_axis)
    matplotlib.pyplot.xlabel(x_name)
    matplotlib.pyplot.ylabel(y_name)
