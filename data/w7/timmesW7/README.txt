 http://cococubed.asu.edu/research_pages/sn1a_w7.shtml
 
 code snippet to read:
 
 masses = np.genfromtxt('massShells.dat')
total = []
xmass = []
for m in masses:
    total.append(m)
    xmass.append(sum(total)/(solMass*1e3))

timesteps = 290
cells = 175
cellsperchunk = 5

def tstepDens(step):
    densf = np.genfromtxt('timmesW7/densEvo.dat', comments=' NSTG', skip_footer=1, skip_header=1)
    denscells = []
    for i in range(0, cells/cellsperchunk):
        chunk = densf[np.ix_(np.arange(timesteps*i + 1, timesteps*(i+1)), np.arange(2,7))]
        denscells = np.hstack((denscells,chunk[step]))
    return [np.power(10, x) for x in denscells]

def tstepTemp(step):
    tempf = np.genfromtxt('timmesW7/tempEvo.dat', comments=' NSTG',  skip_footer=1, skip_header=1)
    tempcells = []
    for i in range(0, cells/cellsperchunk):
        chunk = tempf[np.ix_(np.arange(timesteps*i + 1, timesteps*(i+1)), np.arange(2,7))]
        tempcells = np.hstack((tempcells,chunk[step]))
    return [np.power(10, x) for x in tempcells]

def stepTstamp(step):
    data = np.genfromtxt('timmesW7/tempEvo.dat', comments=' NSTG',  skip_footer=1, skip_header=1)
    chunk = data[np.ix_(np.arange(0, timesteps), [1])]
    return chunk[step][0]

for step in range(0, timesteps, 30):
    plt.semilogy(xmass[3:], tstepDens(step), label='{:2.3f}s'.format(stepTstamp(step)))
plt.legend()

for step in range(0, timesteps, 30):
    plt.semilogy(xmass[3:], tstepTemp(step), label='{:2.3f}s'.format(stepTstamp(step)))
plt.legend()