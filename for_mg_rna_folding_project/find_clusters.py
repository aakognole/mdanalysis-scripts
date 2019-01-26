import MDAnalysis as mda
import MDAnalysis.analysis.rdf as rdf
import MDAnalysis.analysis.distances as dist
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

pdb = ("./mg.distribution.pdb")
u = mda.Universe(pdb)
mgs = u.select_atoms("resname MG")

cutoff = 2.5
contacts = []
N = []
cnt = 0
for i in mgs.atoms.ids:
    nlist = []
    A = mgs.select_atoms('bynum %d' %(i))
    print A.atoms.ids
    for j in mgs.atoms.ids:
        B = mgs.select_atoms('bynum %d' %(j))
        atob = dist.dist(A,B,offset=0)
        if cutoff > atob[2] > 0:
            nlist.append((j))
    contacts.append((cnt, i, len(nlist)))
    nlist = np.array(nlist)
    N.append((nlist))
    cnt += 1
NEIGHBORS = np.array(N)
contact = np.array(contacts)
cont = contact[np.argsort(contact[:,2])]
rcont = cont[::-1]
np.savetxt('max_contact_positions.txt', rcont, fmt='%s')
np.savetxt('neighbor-lists.txt', NEIGHBORS, fmt='%s')
# Filter contacts into distinct clusters
groups = []
centers = []
mlist = []
ik,il = 0,0
for k in rcont[:,0]:
    item = rcont[ik,1]
    imem = rcont[ik,2]
    members = NEIGHBORS[k]
    if not item in mlist:
        if imem > 0:
            sole_members = []
            for l in members:
                member = int(l)
                if not member in mlist:
                    mlist.append((member))
                    sole_members.append((member))
            sole_members = np.array(sole_members)
            groups.append((sole_members))
            centers.append((il, item, len(sole_members)))
            il += 1
    ik += 1
centers = np.array(centers)
groups = np.array(groups)
sort_centers = centers[np.argsort(centers[:,2])]
rcenters = sort_centers[::-1]
# Find the significant cluster centers
mlist = []
cluster = []
ik = 0
for k in rcenters[:,0]:
    item = rcenters[ik,1]
    imem = rcenters[ik,2]
    members = groups[k]
    if not item in mlist:
        # Criteria for significant clusters
        if imem > rcenters[0,2]/10:
            for l in members:
                member = int(l)
                if not member in mlist:
                    mlist.append((member))
            cluster.append((item, imem))
    ik = ik + 1
cluster = np.array(cluster)
np.savetxt('cluster_positions.txt', cluster, fmt='%s')
# Writing the coordinates of cluster centers
rna = u.select_atoms("all and not resname MG")
rna.write("./clusters/rna.pdb")
im = 0
for m in cluster[:,0]:
    item = m
    imem = cluster[im,1]
    centers = u.select_atoms("bynum %d" % (item))
    if imem > 99:
        imem = 99
    centers.tempfactors = imem
    centers.write("./clusters/mg_%d_size_%d.pdb" % (item, imem))
    im = im + 1

