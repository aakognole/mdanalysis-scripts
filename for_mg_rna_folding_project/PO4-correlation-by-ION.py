import MDAnalysis as mda
import MDAnalysis.analysis.rdf as rdf
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sys import argv

first = 15.0
last = 40.1
for umb in np.arange(first, last, 1.0):
    dcds = []
    print umb
    dcds.append(('./at_umb_%s/1/solution.1.merged.dcd' % (umb)))
    dcds.append(('./at_umb_%s/2/solution.2.merged.dcd' % (umb)))
    #print dcds
    u = mda.Universe('./at_umb_%s/prep_system/system.nowater.pdb' % (umb), dcds)
    skip = 5
    print "    Working on %d frames..." % (u.trajectory.n_frames/skip)
    ions = ( "MG", "POT")
    for ion in ions:
        g1 = u.select_atoms("name O1P or name O2P")
        x = np.arange(1, g1.atoms.n_residues+2, 1)
        y = np.arange(1, g1.atoms.n_residues+2, 1)
        Z1 = np.zeros((g1.atoms.n_residues+2,g1.atoms.n_residues+2))
        Z2 = np.zeros((g1.atoms.n_residues+2,g1.atoms.n_residues+2))
        Z3 = np.zeros((g1.atoms.n_residues+2,g1.atoms.n_residues+2))
        z1max, mres1, mres2 = Z1.max(), 0, 0
        z2max, mres1, mres2 = Z2.max(), 0, 0
        z3max, mres1, mres2 = Z3.max(), 0, 0
        g2 = u.select_atoms("name %s" % (ion))
        print "   ", g1.atoms.n_residues, "Phosphates,", g2.atoms.n_residues, "%ss" % (ion)
        if g2.atoms.n_residues < 1:
            print "    No %s ions in the system" % (ion)
            continue
        for ts in u.trajectory[0::skip]:
            for i in g2.atoms.ids:
                pomg = u.select_atoms("(around 4.6 (bynum %d)) and (name O1P or name O2P)" % (i))
                contacts = np.unique(pomg.atoms.resids)
                if len(contacts) > 1:
                    j = 0
                    while j < len(contacts):
                        res2 = contacts[j]
                        for res1 in contacts:
                            if (res1-res2) > 1:
                                Z1[(res1),(res2)] = Z1[(res1),(res2)] + 1
                                if Z1[(res1),(res2)] > z1max:
                                    z1max = Z1[(res1),(res2)]
                                    #print res1, res2, z1max
                        j = j + 1
                if len(contacts) > 1:
                    j = 0
                    while j < len(contacts):
                        res2 = contacts[j]
                        for res1 in contacts:
                            if (res1-res2) > 1:
                                Z2[(res1),(res2)] = Z2[(res1),(res2)] + 1
                                if Z2[(res1),(res2)] > z2max:
                                    z2max = Z2[(res1),(res2)]
                                    #print res1, res2, z2max
                        j = j + 1
                if len(contacts) > 1:
                    j = 0
                    while j < len(contacts):
                        res2 = contacts[j]
                        for res1 in contacts:
                            if (res1-res2) > 1:
                                Z3[(res1),(res2)] = Z3[(res1),(res2)] + 1
                                if Z3[(res1),(res2)] > z3max:
                                    z3max = Z3[(res1),(res2)]
                                    #print res1, res2, z3max
                        j = j + 1
        Z1 = Z1*100/g2.atoms.n_residues/(u.trajectory.n_frames/skip)
        Z2 = Z2*100/g2.atoms.n_residues/(u.trajectory.n_frames/skip)
        Z3 = Z3*100/g2.atoms.n_residues/(u.trajectory.n_frames/skip)
        # For type-1 contacts
        WHERE = np.array(np.where(Z1 > ((Z1.max())/4))).transpose()
        k = 0
        while k < len(WHERE[:,0]):
            print "%3s %3s %10.5f" % (WHERE[k,0],WHERE[k,1],Z[WHERE[k,0],WHERE[k,1]])
            k += 1
        #
        plt.figure()
        plt.title('Phosphates correlated by %s at %s' % (ion, umb))
        plt.imshow(Z1, cmap='coolwarm', interpolation='nearest', origin='lower')
        plt.xlabel('Residue Number')
        plt.xticks(np.arange(1,53,5))
        plt.ylabel('Residue Number')
        plt.yticks(np.arange(1,53,5))
        plt.colorbar()
        plt.savefig('around_%s_at_umb_%s_type-1.png' % (ion, umb), dpi=300)
        plt.close()
        #
        plt.figure()
        plt.title('Phosphates correlated by %s at %s' % (ion, umb))
        plt.imshow(Z2, cmap='coolwarm', interpolation='nearest', origin='lower')
        plt.xlabel('Residue Number')
        plt.xticks(np.arange(1,53,5))
        plt.ylabel('Residue Number')
        plt.yticks(np.arange(1,53,5))
        plt.colorbar()
        plt.savefig('around_%s_at_umb_%s_type-2.png' % (ion, umb), dpi=300)
        plt.close()
        #
        plt.figure()
        plt.title('Phosphates correlated by %s at %s' % (ion, umb))
        plt.imshow(Z3, cmap='coolwarm', interpolation='nearest', origin='lower')
        plt.xlabel('Residue Number')
        plt.xticks(np.arange(1,53,5))
        plt.ylabel('Residue Number')
        plt.yticks(np.arange(1,53,5))
        plt.colorbar()
        plt.savefig('around_%s_at_umb_%s_type-3.png' % (ion, umb), dpi=300)
        plt.close()

exit()

