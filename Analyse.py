import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

mscale = m_earth = 5.9723 * 10**24 / 100

def getData(filename, p):
    '''p: expected number of particles'''
    n = 3 + 7 * p
    objs = []
    with open(filename, "rb") as f:
        data = np.fromfile(f)
    data = np.reshape(data, (-1, n))
    times = data[:, 0] / (365.25 * 24 * 60 * 60)
    kinetic = data[:, 1]
    gravi = data[:, 2]
    total = kinetic + gravi
    masses = data[:, 3:3+p] / mscale
    objs = data[:, 3+p:]
    objs = np.reshape(objs, (-1, p, 6))
    return [times, kinetic, gravi, total, masses, objs]

# withJupiter = [getData(f"SolColJup2{n}.bin", p=102) for n in range(5)]
# # withJupiter = [getData(f"SolColJup{n}.bin", p=102) for n in range(5)]
# withoutJupiter =[getData(f"SolColNoJup2{n}.bin", p=101) for n in range(5)]
# # withoutJupiter = [getData(f"SolColNoJup{n}.bin", p=101) for n in range(5)]
withJupiter = [getData("SolColJupL10000.bin", p=102)]
withoutJupiter =[getData(f"SolColNoJupL10000.bin", p=101)]

sem_ear = 1.4960 * 10**11 # m
scale = (sem_ear, "A.U.")
timescale = (365.25*24*60*60, "years")

def movingAverage(x):
    ...

def averageRagainstT():
    J = withJupiter[0]
    JT = J[0] / 1000
    JM = J[4][:, 2:]
    JPos = J[5][:, 2:, 0:3] / (1.5*10**11)
    JXYZ = np.sqrt(np.sum(JPos**2, axis=2))
    JR = np.average(JXYZ, axis=1, weights=JM)
    
    NJ = withoutJupiter[0]
    NJT = NJ[0] / 1000
    NJM = NJ[4][:, 1:]
    NJPos = NJ[5][:, 1:, 0:3] / (1.5*10**11)
    NJXYZ = np.sqrt(np.sum(NJPos**2, axis=2))
    NJR = np.sqrt(np.average(NJXYZ**2, axis=1, weights=NJM))

    plt.plot(JT[::1000], JR[::1000], label="Jupiter")
    plt.plot(NJT[::1000], NJR[::1000], label="No Jupiter")
    plt.xlabel("Time, kyr")#
    plt.ylabel("Average Orbital Distance, AU")
    plt.legend()
    plt.show()
averageRagainstT()


def plotframe(i=0):
    fig, ax = plt.subplots()
    fig.set_figheight(5)
    fig.set_figwidth(5)
    fig.set_tight_layout(True)

    wh = 1.2
    ax.set_xlim(-wh, wh)
    ax.set_xlabel("x, " + scale[1])
    ax.set_ylim(-wh, wh)
    ax.set_ylabel("y, " + scale[1])

    time = withJupiter[0][0]
    objs = withJupiter[0][5]

    pos = np.array(objs)

    xs = pos[:, :, 1] / scale[0]
    ys = pos[:, :, 2] / scale[0]

    xdata = xs[i]
    ydata = ys[i]
    ax.plot(xdata-xdata[0], ydata-ydata[0], ".")
    ax.text(0.02, 0.95, f"{time[i] / timescale[0]} {timescale[1]}", transform=ax.transAxes)
    plt.show()

def plotAngMoOverTime():
    time = withoutJupiter[0][0]
    masses = withoutJupiter[0][4]
    objs = withoutJupiter[0][5]

    pos = np.array(objs)

    r = pos[:, :, 0:3]
    v = pos[:, :, 3:6]
    p = masses[:,:,None] * v
    L = np.cross(r, p)

    Lx = np.sum(L[:,:, 0], axis=1)
    Ly = np.sum(L[:,:, 1], axis=1)
    Lz = np.sum(L[:,:, 2], axis=1)

    plt.plot(time, Lx, label="x")
    plt.plot(time, Ly, label="y")
    plt.plot(time, Lz, label="z")

    print(f"x: {np.min(Lx)} {np.max(Lx)} {np.max(Lx)-np.min(Lx)} {100*(np.max(Lx)-np.min(Lx)) / np.average(Lx)}")
    print(f"y: {np.min(Ly)} {np.max(Ly)} {np.max(Ly)-np.min(Ly)} {100*(np.max(Ly)-np.min(Ly)) / np.average(Ly)}")
    print(f"z: {np.min(Lz)} {np.max(Lz)} {np.max(Lz)-np.min(Lz)} {100*(np.max(Lz)-np.min(Lz)) / np.average(Lz)}")

    plt.legend()
    plt.show()

# plotAngMoOverTime()
def plotMomentumOverTime():
    time = withoutJupiter[0][0]
    masses = withoutJupiter[0][4]
    objs = withoutJupiter[0][5]

    pos = np.array(objs)

    xs = pos[:, :, 0] / scale[0]
    ys = pos[:, :, 1] / scale[0]
    zs = pos[:, :, 2] / scale[0]
    vxs = pos[:, :, 3] / scale[0]
    pxs = masses * vxs
    tpxs = np.sum(pxs, axis=1)
    vys = pos[:, :, 4] / scale[0]
    pys = masses * vys
    tpys = np.sum(pys, axis=1)
    vzs = pos[:, :, 5] / scale[0]
    pzs = masses * vzs
    tpzs = np.sum(pzs, axis=1)

    plt.plot(time, tpxs, label="x")
    plt.plot(time, tpys, label="y")
    plt.plot(time, tpzs, label="z")

    print(f"x: {np.min(tpxs)} {np.max(tpxs)} {np.max(tpxs)-np.min(tpxs)}")
    print(f"y: {np.min(tpys)} {np.max(tpys)} {np.max(tpys)-np.min(tpys)}")
    print(f"z: {np.min(tpzs)} {np.max(tpzs)} {np.max(tpzs)-np.min(tpzs)}")

    plt.legend()
    plt.show()

# plotMomentumOverTime()

## rate of merging
# number of particles against time
def numParticles(lst: tuple, i=0) -> list:
    '''Returns (times, numparticles) from the getData tuple'''
    #lst: times kin grav total masses objs
    # mass is set to 0 if coupled / no longer a particle
    # -> count non zero masses
    return np.count_nonzero(lst, axis=1) - i

# numPartiJup = np.array([numParticles(J[4], 2) for J in withJupiter])
# avgNumJup = np.average(numPartiJup, axis=0)

# numPartiNoJup = np.array([numParticles(NJ[4], 1) for NJ in withoutJupiter])
# avgNumNoJup = np.average(numPartiNoJup, axis=0)

def RateOfMerging():
    fig, ax = plt.subplots()
    fig.set_figheight(5)
    fig.set_figwidth(5)
    fig.set_tight_layout(True)
    
    TJ = withJupiter[0][0] / 1000
    TW = withoutJupiter[0][0] / 1000
    ax.plot(TJ[::200], avgNumJup[::200], 
             "+r", label="With Jupiter", markersize=1)
    ax.plot(TW[::200], avgNumNoJup[::200], 
             "xb", label="Without Jupiter", markersize=1)

    #   -> curvefit
    def expCurveJ(x, b, c):
        return 100 * np.exp(-b * x - c*x**0.5) # 100 * np.exp(-b * x)
    def expCurveNJ(x, b, c):
        return 100 * np.exp(-b * x - c*x**0.5) # 100 * np.exp(-b* x)
    popt_Jup, pcov = curve_fit(f=expCurveJ, xdata=TJ,
                               ydata=avgNumJup, p0=(1, 1), absolute_sigma=True)
    perr_Jup = np.sqrt(np.diag(pcov))
    popt_NoJup, pcov = curve_fit(f=expCurveNJ, xdata=TW, 
                                 ydata=avgNumNoJup, p0=(1, 1))
    perr_NoJup = np.sqrt(np.diag(pcov))
    ts = np.linspace(0, withJupiter[0][0][-1] / 1000, 1000)
    plt.plot(ts, expCurveJ(ts, *popt_Jup), "-g", linewidth=2,
             label=f"w/ Jup: \nb={popt_Jup[0]:.3f}  c={popt_Jup[1]:.3f}") #
    print(f"Jup: b = {popt_Jup[0]} +/- {perr_Jup[0]}, c = {popt_Jup[1]} +/- {perr_Jup[1]}")
    plt.plot(ts, expCurveNJ(ts, *popt_NoJup), "-m", linewidth=2,
             label=f"w/o Jup: \nb={popt_NoJup[0]:.3f} c={popt_NoJup[1]:.3f}")
    print(f"NoJup: b = {popt_NoJup[0]} +/- {perr_NoJup[0]}, c = {popt_NoJup[1]} +/- {perr_NoJup[1]} ") # 
    plt.xlabel("Time, 1000 years")
    plt.ylabel("Number of particles")
    # plt.yscale("log")
    # plt.xscale("log")
    plt.legend()
    plt.autoscale() 
    plt.show()
    b = popt_Jup[0]
    c = popt_Jup[1]
    tJup1 = (-c + np.sqrt(c**2 + 8 * b * np.log(10)) / (2 * b)) ** 2
    tJup2 = (-c - np.sqrt(c**2 + 8 * b * np.log(10)) / (2 * b)) ** 2
    b = popt_NoJup[0]
    c = popt_NoJup[1] 
    tNoJup1 = (-c + np.sqrt(c**2 + 8 * b * np.log(10)) / (2 * b)) ** 2
    tNoJup2 = (-c - np.sqrt(c**2 + 8 * b * np.log(10)) / (2 * b)) ** 2
    print(f"tJup = {tJup1}, {tJup2};  tNoJup = {tNoJup1} , {tNoJup2}")

# RateOfMerging()
# plotframe(0)