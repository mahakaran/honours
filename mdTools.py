def plot_xvg(xvg_file):
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
    plt.plot(x_axis, y_axis)
