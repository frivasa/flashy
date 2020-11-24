D. Jeffery @ UNLV https://www.nhn.ou.edu/~jeffery/astro/supernovae/spectra/model/w7/w7.html

snippet to read:

import flashy.nuclear as cgs

filename = 'fullW7.dat'
# read only composition

unks = np.genfromtxt('jefferyW7/unks.dat', dtype=None)

labels = [unks[i][0] for i in range(len(unks))]
values = [list(unks[i])[1:] for i in range(len(unks))]

xvars = dict(zip(labels, values))
xvars.keys()

xmass = []
total = []
# for m in np.flipud(xvars['xmass_cell']):
for m in xvars['xmass_cell']:
    total.append(np.power(10, m))
    xmass.append(sum(total))

r = []
# for m in np.flipud(xvars['xmass_cell']):
for m in xvars['radius']:
    r.append(np.power(10, m))

As, Sp = zip(*[cgs.elemSplit(n) for n in elems])
Zs = np.array([massdict[sp]['z'] for sp in Sp])
Ns = As - Zs
wgts = [massdict[sp]['n'][n] for (n, sp) in zip(Ns, Sp)]
atot = sum(wgts)

chars=14
elems = {}
with open('compcells.dat', 'r') as f:
    for line in f:
        if line.startswith('O='):
            continue
        else:
            raw = line.rstrip()
            snips = [raw[i:i+chars] for i in range(0, len(raw), chars)]
            comps = [(s[:4].lower().strip().capitalize(), float(s[5:])) for s in snips]
            for i, (n, v) in enumerate(comps):
                if n in elems:
                    elems[n].append(v)
                else:
                    elems[n.strip()] = []
                    elems[n].append(v)

print 'Species: {}, zones: {}'.format(len(elems), len(elems['Mg26']))
print 'Adding Fluff...'
fluffcomp = {'C12':0.4875, 'O16':0.4875, 'Ne22': 0.025}
flufflen = len(xmass)-len(elems['Mg26'])
for i in range(flufflen):
    for k in elems.keys():
        if k in fluffcomp:
            elems[k].append(fluffcomp[k])
        else:
            elems[k].append(0.0)
print 'Species: {}, zones: {}'.format(len(elems), len(elems['Mg26']))

%matplotlib notebook
plt.semilogy(total, elems['Ni56'])
plt.semilogy(total, elems['C12'])
plt.semilogy(total, elems['O16'])
plt.axhline(1.0, ls='--', alpha=0.4)
plt.ylim([1e-5,2])