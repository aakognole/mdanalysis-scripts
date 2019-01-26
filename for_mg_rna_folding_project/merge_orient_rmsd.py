import MDAnalysis as mda
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

runs = <runs>
cycles = <cycles>

syspdb = ('./prep_system/system.pdb')
t = mda.Universe(syspdb)
nowater = t.select_atoms("all and not resname SOL")
nowater.write('./prep_system/system.nowater.pdb')

for i in range(1,runs+1):
    print "for run =",i
    trjlist = []
    jj = 0
    for ii in range(2,cycles+1,2):
        trjlist.append(('./%d/solution.%d.prod.%d.nowater.dcd' % (i,i,ii)))
    pdb = ('./prep_system/system.nowater.pdb')
    u = mda.Universe(pdb, trjlist)
    ref = mda.Universe(pdb)
    alignment = align.AlignTraj(u, ref, select="nucleic and (name P or name O3' or name C3' or name C4' or name C5' or name O5')", filename="./%d/solution.%d.merged.dcd" % (i,i))
    alignment.run()
    align.alignto(u, ref, select="nucleic and (name P or name O3' or name C3' or name C4' or name C5' or name O5')", weights="mass")
    rna = u.select_atoms("nucleic and (name P or name O3' or name C3' or name C4' or name C5' or name O5')")
    rnaref = ref.select_atoms("nucleic and (name P or name O3' or name C3' or name C4' or name C5' or name O5')")
    R = []
    rgyr = []
    for ts in u.trajectory:
        r = rmsd(rna.positions, rnaref.positions, superposition=True)
        rog = rna.radius_of_gyration()
        time = (ts.frame)/100.0
        R.append((time, r))
        rgyr.append((time, rog))
    R = np.array(R)
    rgyr = np.array(rgyr)
    if i == 1:
        RMSD = np.array(R)
        RGYR = np.array(rgyr)
    else:
        second = np.array([R[:,1]])
        RMSD = np.hstack((RMSD, second.T))
        second = np.array([rgyr[:,1]])
        RGYR = np.hstack((RGYR, second.T))
    print RMSD.shape, RGYR.shape

np.savetxt('RNA_rmsd.txt', RMSD, fmt='%s')
np.savetxt('RNA_rgyr.txt', RGYR, fmt='%s')
RMSD = np.loadtxt('RNA_rmsd.txt')
print RMSD.shape
fig, axes = plt.subplots(runs, 1, sharex=True, figsize=(6,runs*2.5))
fig.subplots_adjust(hspace=0.5)
plt.suptitle('RNA RMSD at reaction coordinate XXX', fontsize=12)
for i in range(1,runs+1):
    jj = i-1
    axes[jj].plot(RMSD[:,0], RMSD[:,i])
    ax0 = axes[jj]
    ax0.set_title("RMSD for Run: %d" % (i), fontsize=10)
    if i == runs:
        ax0.set_xlabel('Time (ns)')
    fig.autofmt_xdate()
plt.savefig('RNA_RMSD_at_umb_XXX.png')

RGYR = np.loadtxt('RNA_rgyr.txt')
print RGYR.shape
fig, axes = plt.subplots(runs, 1, sharex=True, figsize=(6,runs*2.5))
fig.subplots_adjust(hspace=0.5)
plt.suptitle('RNA RGYR at reaction coordinate XXX', fontsize=12)
for i in range(1,runs+1):
    jj = i-1
    axes[jj].plot(RGYR[:,0], RGYR[:,i])
    ax0 = axes[jj]
    ax0.set_title("RGYR for Run: %d" % (i), fontsize=10)
    if i == runs:
        ax0.set_xlabel('Time (ns)')
    fig.autofmt_xdate()
plt.savefig('RNA_RGYR_at_umb_XXX.png')

print "RGYR and RGYR plots done!"
