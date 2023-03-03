# extra
from GravitySimulation import *

def solarSystemTests():
    m_other   = 0.33011 * 10**24 # mercury
    sem_other = 57.91   * 10**9
    # m_other   = 4.8675 * 10**24 # venus
    # sem_other = 108.21 * 10**9
    # m_other   = 0.64171 * 10**24 # mars
    # sem_other = 227.92  * 10**9
    # m_other   = 1898.19 * 10**24 # Jup
    # sem_other = 778.57  * 10**9
    # m_other   = 568.34  * 10**24 # Saturn
    # sem_other = 1433.53 * 10**9
    # m_other   = 86.813  * 10**24 # Uranus
    # sem_other = 2872.46 * 10**9
    # m_other   = 102.413 * 10**24 # Neptune
    # sem_other = 4495.06 * 10**9
    a = System(2 * 24 * 60 * 60, file="SOLARSYSTEST.bin") # 2 days
    v_other = vel(sem_other, sem_other)
    a.AddParticle(M, 0.0, 0.0, 0.0, 0.0, 0, 0.0) # Sun
    a.AddParticle(m_other, sem_other, 0.0, 0.0, 0.0, +v_other, 0.0) # Other
    v_earth = vel(sem_ear, sem_ear)
    a.AddParticle(m_earth, -sem_ear, 0, 0, 0, -v_earth, 0) # Earth
    a.Leapfrog()
    a.Record()
    a.Update()
    dur = 100 * 365.25 / 2
    for n in range(int(dur)):
        a.Record()
        a.doTimestep()
        if n % 1000 == 0:
            print(n, end=" ", flush=True)
    a.Record()
    a.File.close()
    setAnimate(widthheight_=scale[0]*1.2, scale_=scale, 
                timescale_=timescale, tracelength_=tracelength)
    animateFile("SOLARSYSTEST.bin", p=3, frameskip=frameskip, repeat=False, ax=(0,1))
    graphR("SOLARSYSTEST.bin")

def testEnv():
    te = System(1*24*60*60, file="TEST.bin")

    v = vel(aph_ear, aph_ear)
    v2 = vel(aph_ear, 0.75*aph_ear)

    te.AddParticle(M, 0,0,0, 0,-m_earth/2 * v / M + m_earth/2 * v2 / M,0)
    
    te.AddParticle(m_earth/2, aph_ear,0,0, 0,v,0)
    te.AddParticle(m_earth/2, -aph_ear,0,0, 0,-v2,0)
    # te.AddParticle(M, -20*eps,0,0, 0.002*eps,0,0)
    # te.AddParticle(M, -5*eps,0,0, 0,0,0)
    for x in range(365*15):
        te.Record() 
        te.doTimestep()
        #te.doTimestep(0.01*60*60)
    te.Record()
    te.File.close()

    p=3

    setAnimate(widthheight_=scale[0]*1.2, scale_=scale, 
                timescale_=timescale, tracelength_=5)
    animateFile("TEST.bin", p=p, frameskip=10)
    graphEnergies("TEST.bin", change=False, p=p)
    graphEnergies("TEST.bin", change=True, p=p)