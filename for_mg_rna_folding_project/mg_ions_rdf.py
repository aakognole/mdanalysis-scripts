import MDAnalysis as mda
import MDAnalysis.analysis.rdf as rdf
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

pdb = ("./all.mg.positions.pdb")
u = mda.Universe(pdb)

g1 = u.select_atoms("name O1P or name O2P")
g2 = u.select_atoms("name MG")

omg = rdf.InterRDF(g1, g2, nbins=100, range=(0.0,10.0))
omg.run()

plt.figure()
plt.subplot(211)
plt.title("Phosphate Oxygen - MG : Radial Distribution Function")
plt.ylabel("G(r)")
plt.xlabel("O-Mg distance in Angstroms")
plt.grid()
plt.plot(omg.bins, omg.rdf)
plt.subplot(212)
plt.ylabel("Count")
plt.xlabel("O-Mg distance in Angstroms")
plt.grid()
plt.plot(omg.bins, omg.count, 'r-')
plt.savefig('O-MG-rdf.png', dpi=300)
print "Completed ploting... Phosphate Oxygen - MG : Radial Distribution Function"
