import MDAnalysis as mda
import MDAnalysis.analysis.rdf as rdf
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

dcds = []
for i in np.arange(15.0,40.0,2):
    print i
    dcds.append(('./at_umb_%s/1/solution.1.merged.dcd' % (i)))
    dcds.append(('./at_umb_%s/2/solution.2.merged.dcd' % (i)))
#print dcds
u = mda.Universe('./at_umb_15.0/prep_system/system.nowater.pdb', dcds)

g1 = u.select_atoms("name O1P or name O2P or name O6 or name O4 or name O2")
g11 = u.select_atoms("name O1P or name O2P")
g12 = u.select_atoms("nucleic and (name O6 or name O4 or name O2)")

g2 = u.select_atoms("name MG or name POT")
g21 = u.select_atoms("name MG")
g22 = u.select_atoms("name POT")

#g11g2 = rdf.InterRDF(g11, g2, nbins=100, range=(0.0,10.0))
#g11g2.run()
g11g21 = rdf.InterRDF(g11, g21, nbins=100, range=(0.0,10.0))
g11g21.run()
g11g22 = rdf.InterRDF(g11, g22, nbins=100, range=(0.0,10.0))
g11g22.run()

#g12g2 = rdf.InterRDF(g11, g2, nbins=100, range=(0.0,10.0))
#g12g2.run()
g12g21 = rdf.InterRDF(g12, g21, nbins=100, range=(0.0,10.0))
g12g21.run()
g12g22 = rdf.InterRDF(g12, g22, nbins=100, range=(0.0,10.0))
g12g22.run()

plt.title("MG-RDF bb vs base")
plt.plot(g11g21.bins, g11g21.rdf, 'b-', label="backbone")
plt.plot(g12g21.bins, g12g21.rdf, 'r-', label="base")
plt.legend()
plt.savefig("MG-RDF_bb-vs-base.png", dpi=600)
plt.close()

plt.title("POT-RDF bb vs base")
plt.plot(g11g22.bins, g11g22.rdf, 'b-', label="backbone")
plt.plot(g12g22.bins, g12g22.rdf, 'r-', label="base")
plt.legend()
plt.savefig("POT-RDF_bb-vs-base.png", dpi=600)
plt.close()

