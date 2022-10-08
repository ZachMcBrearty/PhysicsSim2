from math import log10
systemtype = ""

if systemtype == "Solar":
    timestep = 1 * 24 * 60 * 60 # seconds, 1 day
    maxtime = 1 * 365.25 * 24 * 60 * 60 # seconds, 5 year#
    soften = 10**6 # used to ensure close interacting bodies do not accelerate infinitely
    frameskip = 1
    G = 6.67 * 10**-11 # N m^2 kg^-2
    scale = 2 * 10**11 # m
    tracelength = 100
elif systemtype == "Globular":
    timestep = 365.25 * 24 * 60 * 60 # seconds, 1 year
    maxtime = 1_000_000_000_000 * 365.25 * 24 * 60 * 60 # seconds, 100,000 years
    soften = 10**8 # used to ensure close interacting bodies do not accelerate infinitely
    frameskip = 25
    G = 6.67 * 10**-11 # N m^2 kg^-2
    scale = 5 * 10 ** 18 # m
    tracelength = 10
else:
    timestep = 0.1
    maxtime = 100
    soften = 0
    frameskip = 10
    G = 1 # Arb Units
    scale = 100
    tracelength = 100

stepcount = maxtime / timestep + 1

Bconst = 1
Aconst = timestep / log10(scale/Bconst + 1)
Cconst = -Aconst*log10(Bconst)
print(stepcount, "steps")