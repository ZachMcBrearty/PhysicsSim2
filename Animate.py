import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import Animation, FuncAnimation

#from Constants import tracelength, scale

fig, ax = plt.subplots()
fig.set_figheight(5)
fig.set_figwidth(5)
fig.set_tight_layout(True)
xs = []
ys = []
times = []
masses = []
xdata = []
ydata = []
ln, = ax.plot([], [], 'ro', ms=5)
timetext = ax.text(0.02, 0.95,"", transform=ax.transAxes)
trace = []

widthheight = 1.2
scale = (1, "m")
timescale = (1, "s")
tracelength = -1

def setAnimate(widthheight_=None, scale_=None, 
               timescale_=None, tracelength_=None):
    global widthheight, scale, timescale, tracelength
    if widthheight_ is not None:
        widthheight = widthheight_
    if scale_ is not None:
        scale = scale_
    if timescale_ is not None:
        timescale = timescale_
    if tracelength_ is not None:
        tracelength = tracelength_

def getData(filename):
    times = [] # 1D
    kinetic = [] # 1D
    gravitational = [] # 1D
    total = [] # 1D

    objs = [] # 3D
    with open(filename, "r") as f:
        for line in f:
            time, energies, positions = line.strip().split(" # ")
            times.append(float(time))

            energies = list(map(float, energies.split(" , ")))
            kinetic.append(energies[0])
            gravitational.append(energies[1])
            total.append(energies[2])#

            pos = positions.split(" * ")
            pos = [list(map(float, q[1:-1].split(", "))) for q in pos]
            objs.append(pos)
    return times, kinetic, gravitational, total, objs

def getDataBinary(filename, p):
    '''p: expected number of particles'''
    n = 3 + 7 * p
    objs = []
    with open(filename, "rb") as f:
        data = np.fromfile(f)
    data = np.reshape(data, (-1, n))
    times = data[:, 0]
    kinetic = data[:, 1]
    gravi = data[:, 2]
    total = kinetic + gravi
    masses = data[:, 3:3+p]
    objs = data[:, 3+p:]
    objs = np.reshape(objs, (-1, p, 6))
    return times, kinetic, gravi, total, masses, objs

def init():
    global trace
    wh = widthheight/scale[0]
    ax.set_xlim(-wh, wh)
    ax.set_xlabel("x, " + scale[1])
    ax.set_ylim(-wh, wh)
    ax.set_ylabel("y, " + scale[1])
    trace = [ax.plot([], [],'-', lw=1)[0] for _ in range(tracenum)]
    return ln, *trace, timetext,

def animate(i):  
    # x, y values to be plotted 
    i = int(i)
    xdata = xs[i]
    ydata = ys[i]
    ln.set_data(xdata-xdata[0]-xdata[0], ydata-ydata[0]-ydata[0])
    if tracelength == -1 or i <= tracelength:
        for j, t in enumerate(trace):
            t.set_data(xs[:i, j]-xdata[0]-xdata[0], ys[:i, j]-ydata[0]-ydata[0])
    else:
        for j, t in enumerate(trace):
            t.set_data(xs[i-tracelength:i+1, j]-xdata[0]-xdata[0], ys[i-tracelength:i+1, j]-ydata[0]-ydata[0])
        
    #timetext.set_text(str(times[i])+"s")
    timetext.set_text(f"{round(times[i], 2)} / {round(times[-1], 2)}" + " " +timescale[1])
    return ln, *trace, timetext

def convMasstoColour(masses):
    return masses

def animateFile(filename, p=2, frameskip=1, ax=(0,1), save=False, **kws):
    global xs, ys, times, tracenum
    if filename.endswith(".txt"):
        times, kin, grav, total, objs = getData(filename)
    elif filename.endswith(".bin"):
        times, kin, grav, total, masses, objs = getDataBinary(filename, p)
    tracenum = len(objs[0])
    pos = np.array(objs)
    
    times = times / timescale[0]
    masses = convMasstoColour(masses)
    xs = pos[:, :, ax[0]] / scale[0]
    ys = pos[:, :, ax[1]] / scale[0]

    ani = FuncAnimation(fig, animate, frames=frameMaker(0, len(pos), frameskip),
                        init_func=init, blit=True, **kws)
    if save:
        ani.save(filename.split(".")[0] + ".gif")
    plt.show()

def graphEnergies(filename, change=False, p=2):
    fig, ax = plt.subplots()
    if filename.endswith(".txt"):
        times, kin, grav, total, objs = getData(filename)
    elif filename.endswith(".bin"):
        times, kin, grav, total, masses, objs = getDataBinary(filename, p)
    t = times / timescale[0]
    if not change:
        ax.plot(t, kin, label="KE")
        ax.plot(t, grav, label="GPE")
        ax.plot(t, total, label="Total Energy")
        ax.set_ylabel("Energy, J")
        ax.set_xlabel("Time, " + timescale[1])
    else:
        te = np.array(total)
        initenergy = te[0]
        change = (te - initenergy)
        percchange = abs(100 * change / initenergy)
        ax.plot(t, percchange, label="Percent change in Energy")
        ax.set_ylabel("Percent change in Energy, %")
        ax.set_xlabel("Time, "+timescale[1])
    plt.legend()
    #plt.get_current_fig_manager().window.raise_()
    fig.canvas.manager.window.raise_()
    plt.show()

def graphR(filename):
    fig, ax = plt.subplots()
    times, *_, objs = getDataBinary(filename, 3)
    t = times / timescale[0]
    objs = objs / scale[0]
    sun = objs[:, 0]
    earth = objs[:, 1] - sun
    Rearth = np.sqrt(np.average(earth**2, axis=1))
    NormREarth = Rearth / Rearth[0]
    print(f"Max: {max(NormREarth)}, Min: {min(NormREarth)}, diff: {max(NormREarth)-min(NormREarth)}")
    ax.plot(t, NormREarth, label="Earth Orbital Distance")
    ax.set_ylabel("Distance, " + scale[1])
    ax.set_xlabel("Time, " + timescale[1])

    plt.legend()
    #plt.get_current_fig_manager().window.raise_()
    fig.canvas.manager.window.raise_()
    plt.show()
def graphRMSr(filename, p=100, percent=False):
    fig, ax = plt.subplots()
    if filename.endswith(".txt"):
        times, *_, objs = getData(filename)
    elif filename.endswith(".bin"):
        times, *_, objs = getDataBinary(filename, p)
    t = times / timescale[0]
    objs = objs / scale[0]
    RMS = np.sqrt(np.average(objs**2, axis=(1, 2)))
    if percent:
        NormRMS = abs(1 - RMS / RMS[0])
        ax.plot(t, NormRMS, label="% change in RMS")
        ax.set_ylabel("Distance, % of start")
    else:
        ax.plot(t, RMS, label="RMS Distance")
        ax.set_ylabel("Distance, " + scale[1])
    ax.set_xlabel("Time, " + timescale[1])

    plt.legend()
    #plt.get_current_fig_manager().window.raise_()
    fig.canvas.manager.window.raise_()
    plt.show()

def frameMaker(start, stop, step):
    for x in range(start, stop, step):
        yield x
    yield stop-1

if __name__=="__main__":
    # f = "PlanetsSim.txt"
    # fbin = "GravitySim.bin"
    # animateFile(fbin, p=2, frameskip=1, repeat=False)
    #graphEnergies(f, True)
    
    sem_ear = 1.4960 * 10**11 # m
    scale = (sem_ear, "A.U.")
    timescale = (365.25*24*60*60, "years")

    
    fs=365
    tracelength = 2

    frameskip = 1
    p=102
    DEFAULTFILE = "SolColJup20.bin"
    setAnimate(widthheight_=scale[0]*1.2, scale_=scale, 
                timescale_=timescale, tracelength_=tracelength)
    animateFile(DEFAULTFILE, p, save=False)
    graphEnergies(DEFAULTFILE, False, p)
    graphEnergies(DEFAULTFILE, True, p)