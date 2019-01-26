import MDAnalysis
import numpy as np
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import MDAnalysis.analysis.rdf as rdf
import MDAnalysis.analysis.distances as dist

pdb = ('solution.<runn>.prod.<curr>.pdb')
dcd = ('solution.<runn>.prod.<curr>.dcd')

u = MDAnalysis.Universe(pdb, dcd)
nowater = u.select_atoms("all and not resname SOL")
nowater.write('solution.<runn>.prod.<curr>.nowater.pdb')
with MDAnalysis.Writer("temp.nowater.dcd", nowater.n_atoms) as W:
    for ts in u.trajectory:
        W.write(nowater)

nwpdb = ('solution.<runn>.prod.<curr>.nowater.pdb')
nw = MDAnalysis.Universe(nwpdb, 'temp.nowater.dcd')
ref = MDAnalysis.Universe(nwpdb)
alignment = align.AlignTraj(nw, ref, select="nucleic and not name H*", filename='solution.<runn>.prod.<curr>.nowater.dcd')
alignment.run()

a = MDAnalysis.Universe(nwpdb, nwpdb, 'solution.<runn>.prod.<curr>.nowater.dcd')
rna = a.select_atoms("nucleic and not name H*")
refrna = ref.select_atoms("nucleic and not name H*")
RMSD = []
for ts in a.trajectory:
    R = rmsd(rna.positions, refrna.positions)
    RMSD.append((ts.frame, R))
RMSD = np.array(RMSD)
np.savetxt('rmsd.<runn>.<curr>.txt', RMSD, fmt='%s')        
