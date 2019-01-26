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

fig, axes = plt.subplots(cycles/2, runs, sharex=True, figsize=(8*runs,cycles/2*2.0))
fig.subplots_adjust(hspace=0.5)
plt.suptitle('Nucleotide-MG-contact at_umb_XXX', fontsize=12)
ii = 0
for i in range(1,runs+1):
    print "run =", i
    jj = 0
    for j in range(2,cycles+1,2):
        print "... cycle =", j/2
        trj = ('./%d/solution.%d.prod.%d.nowater.dcd' % (i,i,j))
        pdb = ('./prep_system/system.nowater.pdb')
        u = mda.Universe(pdb, pdb, trj)
        ref = mda.Universe(pdb)
        align.alignto(u, ref, select="nucleic and (name P or name O3' or name C3' or name C4' or name C5' or name O5')", weights="mass")
        rna = u.select_atoms("nucleic and (name P or name O3' or name C3' or name C4' or name C5' or name O5')")
        rnaref = ref.select_atoms("nucleic and (name P or name O3' or name C3' or name C4' or name C5' or name O5')")
        nucleotides = u.select_atoms("nucleic")
        normalizer = u.trajectory.n_frames
        weight = np.zeros((2,len(nucleotides.residues)))
        for ts in u.trajectory:
            ires = 0
            time = (ts.frame)/100.0
            for res in nucleotides.residues:
                total = u.select_atoms("(sphlayer 0 5.5 (nucleic and resid %d)) and name MG" % (res.resid))
                if total.atoms.n_atoms >= 1:
                    weight[1,ires] += 1
                    #weight[1,ires] = weight[1,ires] + total.atoms.n_atoms
                if ts.frame == 0:
                    weight[0,ires] = res.resid
                ires = ires + 1
#            if ts.frame > 100:
#                break
        print weight.shape
        weight[1,:] = weight[1,:]/normalizer
        axes[jj,ii].bar(weight[0,:], weight[1,:])
        ax0 = axes[jj,ii]
        ax0.set_title("total for Run: %d & Cycle: %d" % (i,j/2), fontsize=10)
        ax0.set_ylim(0,1.05)
        if j == cycles:
            ax0.set_xlabel('Residue Number', fontsize=10)
        jj += 1
        fig.autofmt_xdate()
        np.savetxt("new_nucleotide-mg-total.%d.%d.txt" % (i,j), weight, fmt='%s')
    ii += 1
plt.savefig('new_nucleotide-mg-contact_at_umb_XXX.png')
print "Nucleotide-MG-contact plots done!"
