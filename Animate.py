import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import Animation, FuncAnimation

fig, ax = plt.subplots()
xs = []
ys = []
times = []
xdata = []
ydata = []
ln, = ax.plot([], [], 'ro')
timetext = ax.text(0.05, 0.9,"", transform=ax.transAxes)
trace = []
tracelength = 5000
scale = 2*10**11 # m

def getData(filename):
    times = [] # 1D
    kinetic = [] # 1D
    gravitational = [] # 1D
    total = [] # 1D

    objs = [] # 3D
    with open(filename, "r") as f:
        for line in f:
            time, energies, positions = line.split(" # ")
            times.append(float(time))

            energies = list(map(float, energies.split(" , ")))
            kinetic.append(energies[0])
            gravitational.append(energies[1])
            total.append(energies[2])#

            pos = positions.split(" * ")
            pos = [list(map(float, q[1:-2].split(", "))) for q in pos]
            objs.append(pos)
    return times, kinetic, gravitational, total, objs

def init():
    global trace
    ax.set_xlim(-scale, scale)
    ax.set_xlabel("x, m")
    ax.set_ylim(-scale, scale)
    ax.set_ylabel("y, m")
    trace = [ax.plot([], [],'-', lw=1)[0] for _ in range(tracenum)]
    return ln, *trace, timetext,

def animate(i):  
    # x, y values to be plotted 
    i = int(i)
    xdata = xs[i]
    ydata = ys[i]
    ln.set_data(xdata, ydata)
    if i > tracelength:
        for j, t in enumerate(trace):
            t.set_data(xs[i-tracelength:i, j], ys[i-tracelength:i, j])
    else:
        for j, t in enumerate(trace):
            t.set_data(xs[:i, j], ys[:i, j])
    #timetext.set_text(str(times[i])+"s")
    timetext.set_text(str(times[i] / (1000 * 365.25 * 24 * 60 * 60))+"k.yr")
    100 * 365.25 * 24 * 60 * 60
    return ln, *trace, timetext

def animateFile(filename, framesskip=1, repeat=True, scale_=2*10**11):
    global xs, ys, times, tracenum, scale
    scale = scale_
    times, kin, grav, total, objs = getData(filename)
    tracenum = len(objs[0])
    pos = np.array(objs)
    
    xs = pos[:, :, 0]
    ys = pos[:, :, 1]    

    ani = FuncAnimation(fig, animate, frames=range(0, len(pos), framesskip),
                        init_func=init, blit=True, repeat=repeat)
    
    plt.show()

def graphEnergies(filename, change=False):
    fig, ax = plt.subplots()
    times, kin, grav, total, _ = getData(filename)
    if not change:
        ax.plot(times, kin, label="KE")
        ax.plot(times, grav, label="GPE")
        ax.plot(times, total, label="Total Energy")
    else:
        te = np.array(total)
        initenergy = te[0]
        change = te - initenergy
        percchange = change / initenergy

        ax.plot(times, percchange, label="PercentChange in Total")
    plt.legend()
    plt.show()


if __name__=="__main__":
    f = "PlanetsPositionsEnergy57-18-04-10-22.txt"
    animateFile(f, framesskip=50, repeat=False, scale_=500 * 10**16 * 1.3)
    graphEnergies(f)
    graphEnergies(f, True)