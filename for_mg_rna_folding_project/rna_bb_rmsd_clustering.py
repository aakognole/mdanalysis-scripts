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

def rmsd_clustering(pdb, trj, jobname, cutoff, skip, align_to):
    import MDAnalysis as mda
    import numpy as np
    import matplotlib.pyplot as plt
    from MDAnalysis.analysis import align
    from MDAnalysis.analysis.rms import rmsd
    u1 = mda.Universe(pdb, pdb, trj)
    u2 = mda.Universe(pdb, pdb, trj)
    ref = mda.Universe(pdb)
    align.alignto(u1, ref, select=align_to, weights="mass")
    align.alignto(u2, ref, select=align_to, weights="mass")
    R = np.zeros((u1.trajectory.n_frames, u2.trajectory.n_frames))
    N = []
    C = []
    cnt = 0
    for ts1 in u1.trajectory[0::skip]:
        ref1 = u1.select_atoms(align_to)
        nlist = []
        for ts2 in u2.trajectory[0::skip]:
            ref2 = u2.select_atoms(align_to)
            r = rmsd(ref2.positions, ref1.positions) #, superposition=True)
            time = (ts2.frame)/100.0
            R[ts1.frame,ts2.frame] = r
            if r < cutoff:
                nlist.append((ts2.frame))
        nlist = np.array(nlist)
        N.append((nlist))
        C.append((cnt, ts1.frame, len(nlist)))
        cnt += 1
    RMSD = np.array(R)
    NEIGHBORS = np.array(N)
    COUNT = np.array(C)
    #np.savetxt('rmsd-clustering.%s.txt' % (jobname), RMSD, fmt='%5.3f')
    np.savetxt('neighbor-lists.%s.txt' % (jobname), NEIGHBORS, fmt='%s')
    sort_count = COUNT[np.argsort(COUNT[:,2])]
    rcount = sort_count[::-1]
    np.savetxt('max-neighbors.%s.txt' % (jobname), rcount, fmt='%s')
    # Filter the neighbor lists into clusters
    groups = []
    centers = []
    mlist = []
    ik,il = 0,0
    for k in rcount[:,0]:
        frame = rcount[ik,1]
        imem = rcount[ik,2]
        members = NEIGHBORS[k]
        if not frame in mlist:
            if imem > 0: #rcount[0,2]/5:
                sole_members = []                
                for l in members:
                    member = int(l)
                    if not member in mlist:
                        mlist.append((member))
                        sole_members.append((member))
                sole_members = np.array(sole_members)
                coverage = int(float(len(mlist))/len(rcount[:,0])*100)
                groups.append((sole_members))
                centers.append((il,frame,len(sole_members)))
                il += 1
            else:
                pass
        ik += 1
    centers = np.array(centers)
    groups = np.array(groups)
    sort_centers = centers[np.argsort(centers[:,2])]
    rcenters = sort_centers[::-1]
    # Finding the significant cluster centers
    mlist = []
    cluster = []
    ik,coverage = 0,0
    for k in rcenters[:,0]:
        frame = rcenters[ik,1]
        imem = rcenters[ik,2]
        members = groups[k]
        if not frame in mlist:
            # Criteria for significant clusters
            if coverage < 95  or imem > rcenters[0,2]/10:
                for l in members:
                    member = int(l)
                    if not member in mlist:
                        mlist.append((member))
                coverage = int(float(len(mlist))/u1.trajectory.n_frames*100*skip)
                cluster.append((k, frame, imem, coverage))
            else:
                pass
        ik = ik + 1
    cluster = np.array(cluster)
    np.savetxt('list-of-cluster-centers.%s.txt' % (jobname), cluster, fmt='%s')
    ncluster = len(cluster[:,0])
    print "--->", ncluster, "significant clusters found in '%s' for given cutoff of %d" % (jobname,cutoff)
    print "[(count)","(frame)","(#members)","(cumulative coverage %)]"
    print cluster
    # Writing the coordinates of cluster centers 
    full = u1.select_atoms("all")
    with mda.Writer("cluster-centers.%s.pdb" % (jobname), full.n_atoms) as W:
        ik,temp = 0,0
        for frame in cluster[:,1]:
            temp = int(cluster[ik,3]) - temp
            u1.trajectory[frame]
            full.tempfactors = temp
            W.write(full)
            temp = int(cluster[ik,3])
            ik += 1

cutoff = 2.5
skip = 5
align_to = "nucleic and (name P or name O3' or name C3' or name C4' or name C5' or name O5')"

for i in range(1,runs+1):
    print "for run =",i
    jobname = ("at_umb_XXX.run.%d" % (i))
    trj = ('./%d/solution.%d.merged.dcd' % (i,i))
    pdb = ('./prep_system/system.nowater.pdb')
    rmsd_clustering(pdb, trj, jobname, cutoff, skip, align_to)
