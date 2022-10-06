import matplotlib.pyplot as plt
from numpy import array
fig, ax = plt.subplots()

#i = []
ke = []
gpe = []
te = []
step = []
with open("Energy.txt") as f:
    for it, line in enumerate(f):
        #i.append(0.1*it)
        d = line.split(",")
        ke.append(float(d[0]))
        gpe.append(float(d[1]))
        te.append(float(d[2]))
        if len(step) > 0:
            step.append(step[-1] + float(d[3]))
        else:
            step.append(float(d[3]))
te = array(te)
# ax.plot(step, ke, label="KE")
# ax.plot(step, gpe, label="GPE")
# ax.plot(step, te, label="Total Energy")
initenergy = te[0]
change = te - initenergy
percchange = change / initenergy

ax.plot(step, percchange, label="PercentChange in Total")
#ax.legend()
plt.legend()
plt.show()