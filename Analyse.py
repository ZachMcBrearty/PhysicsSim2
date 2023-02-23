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
withJupiter = [getData("SolColJupL.bin", p=102)]
withoutJupiter =[getData(f"SolColNoJupL.bin", p=101)]

# mass distribution 
# multiple lines, one for each mass? (discrete in M_earth/100)
# massesJup = np.array([J[4] for J in withJupiter])
# avgMassesJup = np.average(massesJup, axis=0)

# massesNoJup = np.array([NJ[4] for NJ in withoutJupiter])
# avgMassesNoJup = np.average(massesNoJup, axis=0)

# plt.plot(withJupiter[0][0][::100], avgMassesJup[::100], 
#          "+r", label="With Jupiter", )
# plt.plot(withoutJupiter[0][0][::100], avgMassesNoJup[::100], 
#          "xb", label="Without Jupiter")
# plt.ylim(0, 20)
# plt.show()
# # # spread over time
#     # calc orbital dist for each particle against time
#     #   -> statistics for spread against time

    



## rate of merging
# number of particles against time
def numParticles(lst: tuple, i=0) -> list:
    '''Returns (times, numparticles) from the getData tuple'''
    #lst: times kin grav total masses objs
    # mass is set to 0 if coupled / no longer a particle
    # -> count non zero masses
    return np.count_nonzero(lst, axis=1) - i

numPartiJup = np.array([numParticles(J[4], 2) for J in withJupiter])
avgNumJup = np.average(numPartiJup, axis=0)

numPartiNoJup = np.array([numParticles(NJ[4], 1) for NJ in withoutJupiter])
avgNumNoJup = np.average(numPartiNoJup, axis=0)

def RateOfMerging():
    fig, ax = plt.subplots()
    fig.set_figheight(5)
    fig.set_figwidth(5)
    fig.set_tight_layout(True)
    # for e in range(5):
    #     ax.plot(withJupiter[0][0][::100], numPartiJup[e, ::100], 
    #             "+r", label="With Jupiter", )
    #     ax.plot(withoutJupiter[0][0][::100], numPartiNoJup[e, ::100], 
    #             "xb", label="Without Jupiter")
    ax.plot(withJupiter[0][0][::100], avgNumJup[::100], 
             "+r", label="With Jupiter", )
    ax.plot(withoutJupiter[0][0][::100], avgNumNoJup[::100], 
             "xb", label="Without Jupiter")

    #   -> curvefit
    # def expCurveJ(x, a, b, c):
    #     return a * np.exp(-b * x) + c
    # def expCurveNJ(x, a, b, c):
    #     return a * np.exp(-b * x) + c
    # popt_Jup, pcov = curve_fit(f=expCurveJ, xdata=withJupiter[0][0],
    #                            ydata=avgNumJup, p0=(100, 1, 3), absolute_sigma=True)
    # popt_NoJup, pcov = curve_fit(f=expCurveNJ, xdata=withoutJupiter[0][0], 
    #                              ydata=avgNumNoJup, p0=(100, 1, 2))
    # ts = np.linspace(0, withJupiter[0][0][-1], 1000)
    # plt.plot(ts, expCurveJ(ts, *popt_Jup),
    #          label=f"w/ Jup: a={popt_Jup[0]:.3f} \nb={popt_Jup[1]:.3f} c={popt_Jup[2]:.3f}")
    # plt.plot(ts, expCurveNJ(ts, *popt_NoJup),
    #          label=f"w/o Jup: a={popt_NoJup[0]:.3f} \nb={popt_NoJup[1]:.3f} c={popt_NoJup[2]:.3f}")
    plt.xlabel("Time, years")
    plt.ylabel("Number of particles")
    plt.legend()
    plt.autoscale()
    plt.show()

RateOfMerging()