import MDAnalysis as mda
import MDAnalysis.analysis.rdf as rdf
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

pdb = ("./all.mg.positions.pdb")
u = mda.Universe(pdb)

#g1 = u.select_atoms("segid RNAA and (name C3' or name C4' or name C5' or name O3' or name O5' or name P)")
g1 = u.select_atoms("name O1P or name O2P")
#print g1.atoms.n_atoms
g2 = u.select_atoms("name MG")
#print g2.atoms.n_atoms
#print u.trajectory.n_frames

omg = rdf.InterRDF(g1, g2, nbins=100, range=(0.0,10.0))
omg.run()

plt.figure(figsize=(12,12))
plt.subplot(211)
plt.title("Phosphate Oxygen - MG : Radial Distribution Function", fontsize=18)
plt.ylabel("G(r)", fontsize=18)
plt.xlabel("O-Mg distance in Angstroms", fontsize=12)
plt.grid()
plt.plot(omg.bins, omg.rdf)
plt.subplot(212)
plt.ylabel("Count", fontsize=18)
plt.xlabel("O-Mg distance in Angstroms", fontsize=12)
plt.grid()
plt.plot(omg.bins, omg.count, 'r-')
plt.savefig('O-MG-rdf.png')
print "Completed ploting... Phosphate Oxygen - MG : Radial Distribution Function"
