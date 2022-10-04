import matplotlib.pyplot as plt

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


ax.plot(step, ke, label="KE")
ax.plot(step, gpe, label="GPE")
ax.plot(step, te, label="Total Energy")
#ax.legend()
plt.legend()
plt.show()
