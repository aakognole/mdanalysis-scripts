import MDAnalysis as mda
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

runs = <runs>
cycles = <cycles>

#syspdb = ('./prep_system/system.pdb')
#t = mda.Universe(syspdb)
#nowater = t.select_atoms("all and not resname SOL")
#nowater.write('./prep_system/system.nowater.pdb')

for i in range(1,runs+1):
    print "for run =",i
    trj = ('./%d/solution.%d.merged.dcd' % (i,i))
    pdb = ('./prep_system/system.nowater.pdb')
    u = mda.Universe(pdb, trj)
    #ref = mda.Universe(pdb)
    #align.alignto(u, ref, select="nucleic and (name P or name O3' or name C3' or name C4' or name C5' or name O5')", weights="mass")
    omgs = u.select_atoms("name O1P or name O2P or name O3' or name O5' or name O4' or name O6 or name O4")
    weight1 = []
    weight2 = []
    for ts in u.trajectory:
        cnt1 = 0
        cnt2 = 0
        time = (ts.frame)/100.0
        for j in omgs.atoms.ids:
            direct = u.select_atoms("(sphlayer 0 2.5 (bynum %d)) and name MG" % (j))
            indirect = u.select_atoms("(sphlayer 2.5 5.5 (bynum %d)) and name MG" % (j))
            if direct.atoms.n_atoms >= 1:
               cnt1 = cnt1 + direct.atoms.n_atoms
            if indirect.atoms.n_atoms >= 1:
               cnt2 = cnt2 + indirect.atoms.n_atoms
        weight1.append((time, cnt1))
        weight2.append((time, cnt2))
    weight1 = np.array(weight1)
    weight2 = np.array(weight2)
    if i == 1:
        WT1 = np.array(weight1)
        WT2 = np.array(weight2)
    else:
        second = np.array([weight1[:,1]])
        WT1 = np.hstack((WT1, second.T))
        second = np.array([weight2[:,1]])
        WT2 = np.hstack((WT2, second.T))
    print WT1.shape, WT2.shape

np.savetxt('O-MG-direct.txt', WT1, fmt='%s')
np.savetxt('O-MG-indirect.txt', WT2, fmt='%s')
WT1 = np.loadtxt('O-MG-direct.txt')
WT2 = np.loadtxt('O-MG-indirect.txt')
print WT1.shape, WT2.shape

fig, axes = plt.subplots(runs*2, 1, sharex=True, figsize=(6,runs*3.5))
fig.subplots_adjust(hspace=1)
plt.suptitle('RNA-MG-contact at reaction coordinate XXX', fontsize=12)
for i in range(1,runs+1):
    jj = (i-1)*2
    axes[jj].plot(WT1[:,0], WT1[:,i])
    axes[jj+1].plot(WT2[:,0], WT2[:,i])
    ax1 = axes[jj]
    ax1.set_title("O-Mg direct for Run: %d" % (i), fontsize=10)
    ax2 = axes[jj+1]
    ax2.set_title("O-Mg indirect for Run: %d" % (i), fontsize=10)
    ax2.set_xlabel('Time (ns)', fontsize=10)
    fig.autofmt_xdate()
plt.savefig('O-MG-contacts_at_umb_XXX.png')
print "RNA-O-MG-contacts plots done!"
