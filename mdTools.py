import matplotlib.pyplot
def plot_xvg(xvg_file, x_name, y_name):
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
