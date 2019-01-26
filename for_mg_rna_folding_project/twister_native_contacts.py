import MDAnalysis as mda
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#from sys import argv
#import glob

#dcds = glob.glob("./openm*/step5*dcd")
#print len(dcds)
first = 15.0
last = 40.1

def dist(v1,v2):
    d = np.sqrt((np.square(v1[0]-v2[0])+np.square(v1[1]-v2[1])+np.square(v1[2]-v2[2])))
    return d

WT1 = []
WT2 = []
for umb in np.arange(first, last, 1.0):
    dcds = []
    print umb
    dcds.append(('./at_umb_%s/1/solution.1.merged.dcd' % (umb)))
    dcds.append(('./at_umb_%s/2/solution.2.merged.dcd' % (umb)))
    u = mda.Universe('./at_umb_%s/prep_system/system.nowater.pdb' % (umb), dcds)
    skip = 1
    print "Working on %d frames..." % (u.trajectory.n_frames/skip)
    po4s = u.select_atoms("nucleic and (name P O1P O2P)")
    bases = u.select_atoms("nucleic and not (name P O1P O2P H5T O5' C5' H5' H5'' C4' H4' O4' C1' H1' C2' H2'' O2' H2' C3' H3' O3' H3T H*)")
    # Secondary and tertiary contacts
    dt0 = [10.0]
    t01 = [24]
    t02 = [29]
    T0 = np.zeros((u.trajectory.n_frames/skip, len(t01)+1))
    dt1 = [6.5, 6.5, 6.5, 8.0]
    t11 = [25, 26, 27, 28]
    t12 = [49, 48, 47, 46]
    T1 = np.zeros((u.trajectory.n_frames/skip, len(t11)+1))
    dt2 = [6.5, 6.5]
    t21 = [13, 14]
    t22 = [31, 30]
    T2 = np.zeros((u.trajectory.n_frames/skip, len(t21)+1))
    dp1 = [14.0, 6.5, 6.5, 6.5, 7.5]
    p11 = [1, 2, 3, 4, 5]
    p12 = [54, 53, 52, 51, 50]
    P1 = np.zeros((u.trajectory.n_frames/skip, len(p11)+1))
    dp2 = [9.0, 6.5, 6.5, 8.5, 7.5]
    p21 = [8, 9, 10, 11, 12]
    p22 = [45, 44, 43, 42, 41]
    P2 = np.zeros((u.trajectory.n_frames/skip, len(p21)+1))
    dp4 = [6.5, 6.5, 6.5, 8.5]
    p41 = [20, 21, 22, 23]
    p42 = [35, 34, 33, 32]
    P4 = np.zeros((u.trajectory.n_frames/skip, len(p41)+1))
    wt_T0, wt_T1, wt_T2, wt_P1, wt_P2, wt_P4 = 0, 0, 0, 0, 0, 0
    # Coordination sites
    dpp1 = [9.5]
    pp11 = [20]
    pp12 = [30]
    PP1 = np.zeros((u.trajectory.n_frames/skip, len(pp11)+1))
    dpp2 = [10.5, 13.0]
    pp21 = [32, 32]
    pp22 = [41, 42]
    PP2 = np.zeros((u.trajectory.n_frames/skip, len(pp21)+1))
    dpp3 = [12.0]
    pp31 = [9]
    pp32 = [28]
    PP3 = np.zeros((u.trajectory.n_frames/skip, len(pp31)+1))
    dpp4 = [9.0]
    pp41 = [6]
    pp42 = [25]
    PP4 = np.zeros((u.trajectory.n_frames/skip, len(pp41)+1))
    dpp5 = [14.0, 14.0]
    pp51 = [4, 5]
    pp52 = [47, 46]
    PP5 = np.zeros((u.trajectory.n_frames/skip, len(pp51)+1))
    wt_PP1, wt_PP2, wt_PP3, wt_PP4, wt_PP5 = 0, 0, 0, 0, 0
    for ts in u.trajectory[0::skip]:
        #print ts.frame/skip
        j = 0
        T0[ts.frame/skip,j] = ts.frame/skip
        while j < len(t01):
            res1 = bases.select_atoms("resid %d" % (t01[j]))
            res2 = bases.select_atoms("resid %d" % (t02[j]))
            T0[ts.frame/skip,j+1] = dist(res1.center_of_mass(),res2.center_of_mass())
            if T0[ts.frame/skip,j+1] < dt0[j]:
                wt_T0 = wt_T0 + (1/float(len(t01)))
            j += 1
        j = 0
        T1[ts.frame/skip,j] = ts.frame/skip
        while j < len(t11):
            res1 = bases.select_atoms("resid %d" % (t11[j]))
            res2 = bases.select_atoms("resid %d" % (t12[j]))
            T1[ts.frame/skip,j+1] = dist(res1.center_of_mass(),res2.center_of_mass())
            if T1[ts.frame/skip,j+1] < dt1[j]:
                wt_T1 = wt_T1 + (1/float(len(t11)))
            j += 1
        j = 0
        T2[ts.frame/skip,j] = ts.frame/skip
        while j < len(t21):
            res1 = bases.select_atoms("resid %d" % (t21[j]))
            res2 = bases.select_atoms("resid %d" % (t22[j]))
            T2[ts.frame/skip,j+1] = dist(res1.center_of_mass(),res2.center_of_mass())
            if T2[ts.frame/skip,j+1] < dt2[j]:
                wt_T2 = wt_T2 + (1/float(len(t21)))
            j += 1
        j = 0
        P1[ts.frame/skip,j] = ts.frame/skip
        while j < len(p11):
            res1 = bases.select_atoms("resid %d" % (p11[j]))
            res2 = bases.select_atoms("resid %d" % (p12[j]))
            P1[ts.frame/skip,j+1] = dist(res1.center_of_mass(),res2.center_of_mass())
            if P1[ts.frame/skip,j+1] < dp1[j]:
                wt_P1 = wt_P1 + (1/float(len(p11)))
            j += 1
        j = 0
        P2[ts.frame/skip,j] = ts.frame/skip
        while j < len(p21):
            res1 = bases.select_atoms("resid %d" % (p21[j]))
            res2 = bases.select_atoms("resid %d" % (p22[j]))
            P2[ts.frame/skip,j+1] = dist(res1.center_of_mass(),res2.center_of_mass())
            if P2[ts.frame/skip,j+1] < dp2[j]:
                wt_P2 = wt_P2 + (1/float(len(p21)))
            j += 1
        j = 0
        P4[ts.frame/skip,j] = ts.frame/skip
        while j < len(p41):
            res1 = bases.select_atoms("resid %d" % (p41[j]))
            res2 = bases.select_atoms("resid %d" % (p42[j]))
            P4[ts.frame/skip,j+1] = dist(res1.center_of_mass(),res2.center_of_mass())
            if P4[ts.frame/skip,j+1] < dp4[j]:
                wt_P4 = wt_P4 + (1/float(len(p41)))
            j += 1
        j = 0
        PP1[ts.frame/skip,j] = ts.frame/skip
        while j < len(pp11):
            res1 = po4s.select_atoms("resid %d" % (pp11[j]))
            res2 = po4s.select_atoms("resid %d" % (pp12[j]))
            PP1[ts.frame/skip,j+1] = dist(res1.center_of_mass(),res2.center_of_mass())
            if PP1[ts.frame/skip,j+1] < dpp1[j]:
                wt_PP1 = wt_PP1 + (1/float(len(pp11)))
            j += 1
        j = 0
        PP2[ts.frame/skip,j] = ts.frame/skip
        while j < len(pp21):
            res1 = po4s.select_atoms("resid %d" % (pp21[j]))
            res2 = po4s.select_atoms("resid %d" % (pp22[j]))
            PP2[ts.frame/skip,j+1] = dist(res1.center_of_mass(),res2.center_of_mass())
            if PP2[ts.frame/skip,j+1] < dpp2[j]:
                wt_PP2 = wt_PP2 + (1/float(len(pp21)))
            j += 1
        j = 0
        PP3[ts.frame/skip,j] = ts.frame/skip
        while j < len(pp31):
            res1 = po4s.select_atoms("resid %d" % (pp31[j]))
            res2 = po4s.select_atoms("resid %d" % (pp32[j]))
            PP3[ts.frame/skip,j+1] = dist(res1.center_of_mass(),res2.center_of_mass())
            if PP3[ts.frame/skip,j+1] < dpp3[j]:
                wt_PP3 = wt_PP3 + (1/float(len(pp31)))
            j += 1
        j = 0
        PP4[ts.frame/skip,j] = ts.frame/skip
        while j < len(pp41):
            res1 = po4s.select_atoms("resid %d" % (pp41[j]))
            res2 = po4s.select_atoms("resid %d" % (pp42[j]))
            PP4[ts.frame/skip,j+1] = dist(res1.center_of_mass(),res2.center_of_mass())
            if PP4[ts.frame/skip,j+1] < dpp4[j]:
                wt_PP4 = wt_PP4 + (1/float(len(pp41)))
            j += 1
        j = 0
        PP5[ts.frame/skip,j] = ts.frame/skip
        while j < len(pp51):
            res1 = po4s.select_atoms("resid %d" % (pp51[j]))
            res2 = po4s.select_atoms("resid %d" % (pp52[j]))
            PP5[ts.frame/skip,j+1] = dist(res1.center_of_mass(),res2.center_of_mass())
            if PP5[ts.frame/skip,j+1] < dpp5[j]:
                wt_PP5 = wt_PP5 + (1/float(len(pp51)))
            j += 1
    print "T0 weight =", (wt_T0/float(u.trajectory.n_frames/skip)), "T1 weight =", (wt_T1/float(u.trajectory.n_frames/skip)), "T2 weight =", (wt_T2/float(u.trajectory.n_frames/skip))
    print "P1 weight =", (wt_P1/float(u.trajectory.n_frames/skip)), "P2 weight =", (wt_P2/float(u.trajectory.n_frames/skip)), "P4 weight =", (wt_P4/float(u.trajectory.n_frames/skip))
    print "PP1 weight =", (wt_PP1/float(u.trajectory.n_frames/skip)), "PP2 weight =", (wt_PP2/float(u.trajectory.n_frames/skip)), "PP3 weight =", (wt_PP3/float(u.trajectory.n_frames/skip)), "PP4 weight =", (wt_PP4/float(u.trajectory.n_frames/skip)), "PP5 weight =", (wt_PP5/float(u.trajectory.n_frames/skip))
    if umb == 15.0:
        print "T0 ---"
        print "mean :", T0[:,1].mean()
        print "std :", T0[:,1].std()
        print "cut-off :", T0[:,1].mean()+T0[:,1].std()*3
        print "T1 ---"
        print "mean :", T1[:,1].mean(),T1[:,2].mean(),T1[:,3].mean(),T1[:,4].mean()
        print "std :", T1[:,1].std(),T1[:,2].std(),T1[:,3].std(),T1[:,4].std()
        print "cut-off :", T1[:,1].mean()+T1[:,1].std()*3,T1[:,2].mean()+T1[:,2].std()*3,T1[:,3].mean()+T1[:,3].std()*3,T1[:,4].mean()+T1[:,4].std()*3
        print "T2 ---"
        print "mean :", T2[:,1].mean(),T2[:,2].mean()
        print "std :", T2[:,1].std(),T2[:,2].std()
        print "cut-off :", T2[:,1].mean()+T2[:,1].std()*3,T2[:,2].mean()+T2[:,2].std()*3
        print "P1 ---"
        print "mean :", P1[:,1].mean(),P1[:,2].mean(),P1[:,3].mean(),P1[:,4].mean(),P1[:,5].mean()
        print "std :", P1[:,1].std(),P1[:,2].std(),P1[:,3].std(),P1[:,4].std(),P1[:,5].std()
        print "cut-off :", P1[:,1].mean()+P1[:,1].std()*3,P1[:,2].mean()+P1[:,2].std()*3,P1[:,3].mean()+P1[:,3].std()*3,P1[:,4].mean()+P1[:,4].std()*3,P1[:,5].mean()+P1[:,5].std()*3
        print "P2 ---"
        print "mean :", P2[:,1].mean(),P2[:,2].mean(),P2[:,3].mean(),P2[:,4].mean(),P2[:,5].mean()
        print "std :", P2[:,1].std(),P2[:,2].std(),P2[:,3].std(),P2[:,4].std(),P2[:,5].std()
        print "cut-off :", P2[:,1].mean()+P2[:,1].std()*3,P2[:,2].mean()+P2[:,2].std()*3,P2[:,3].mean()+P2[:,3].std()*3,P2[:,4].mean()+P2[:,4].std()*3,P2[:,5].mean()+P2[:,5].std()*3
        print "P4 ---"
        print "mean :", P4[:,1].mean(),P4[:,2].mean(),P4[:,3].mean(),P4[:,4].mean()
        print "std :", P4[:,1].std(),P4[:,2].std(),P4[:,3].std(),P4[:,4].std()
        print "cut-off :", P4[:,1].mean()+P4[:,1].std()*3,P4[:,2].mean()+P4[:,2].std()*3,P4[:,3].mean()+P4[:,3].std()*3,P4[:,4].mean()+P4[:,4].std()*3
        print "PP1 ---"
        print "mean :", PP1[:,1].mean()
        print "std :", PP1[:,1].std()
        print "cut-off :", PP1[:,1].mean()+PP1[:,1].std()*3
        print "PP2 ---"
        print "mean :", PP2[:,1].mean(),PP2[:,2].mean()
        print "std :", PP2[:,1].std(),PP2[:,2].std()
        print "cut-off :", PP2[:,1].mean()+PP2[:,1].std()*3,PP2[:,2].mean()+PP2[:,2].std()*3
        print "PP3 ---"
        print "mean :", PP3[:,1].mean()
        print "std :", PP3[:,1].std()
        print "cut-off :", PP3[:,1].mean()+PP3[:,1].std()*3
        print "PP4 ---"
        print "mean :", PP4[:,1].mean()
        print "std :", PP4[:,1].std()
        print "cut-off :", PP4[:,1].mean()+PP4[:,1].std()*3
        print "PP5 ---"
        print "mean :", PP5[:,1].mean(),PP5[:,2].mean()
        print "std :", PP5[:,1].std(),PP5[:,2].std()
        print "cut-off :", PP5[:,1].mean()+PP5[:,1].std()*3,PP5[:,2].mean()+PP5[:,2].std()*3

    WT1.append((umb, (wt_T0/float(u.trajectory.n_frames/skip)), (wt_T1/float(u.trajectory.n_frames/skip)), (wt_T2/float(u.trajectory.n_frames/skip)), (wt_P1/float(u.trajectory.n_frames/skip)), (wt_P2/float(u.trajectory.n_frames/skip)), (wt_P4/float(u.trajectory.n_frames/skip))))
    WT2.append((umb, (wt_PP1/float(u.trajectory.n_frames/skip)), (wt_PP2/float(u.trajectory.n_frames/skip)), (wt_PP3/float(u.trajectory.n_frames/skip)), (wt_PP4/float(u.trajectory.n_frames/skip)), (wt_PP5/float(u.trajectory.n_frames/skip))))
    umb += 1.0

WT1 = np.array(WT1)
WT2 = np.array(WT2)

np.savetxt('native-contact-weights.txt', WT1)
plt.plot(WT1[:,0], WT1[:,1], 'xkcd:blue', label="T0")
plt.plot(WT1[:,0], WT1[:,2], 'xkcd:red', label="T1")
plt.plot(WT1[:,0], WT1[:,3], 'xkcd:green', label="T2")
plt.ylabel('fraction of contact')
plt.ylim(0,1)
plt.xlabel('Reaction coordinate')
plt.legend()
plt.savefig('T0_T1_T2-weights.png', dpi=300)
plt.close()
plt.plot(WT1[:,0], WT1[:,4], 'xkcd:blue', label="P1")
plt.plot(WT1[:,0], WT1[:,5], 'xkcd:red', label="P2")
plt.plot(WT1[:,0], WT1[:,6], 'xkcd:green', label="P4")
plt.ylabel('fraction of contact')
plt.ylim(0,1)
plt.xlabel('Reaction coordinate')
plt.legend(fontsize=10)
plt.savefig('P1_P2_P4-weights.png', dpi=300)
plt.close()

np.savetxt('po4-contact-weights.txt', WT2)
plt.plot(WT2[:,0], WT2[:,1], 'xkcd:blue', label="PP1")
plt.plot(WT2[:,0], WT2[:,2], 'xkcd:red', label="PP2")
plt.plot(WT2[:,0], WT2[:,3], 'xkcd:green', label="PP3")
plt.plot(WT2[:,0], WT2[:,4], 'xkcd:pink', label="PP4")
plt.plot(WT2[:,0], WT2[:,5], 'xkcd:black', label="PP5")
plt.ylabel('fraction of contact')
plt.ylim(0,1)
plt.xlabel('Reaction coordinate')
plt.legend(fontsize=10)
plt.savefig('PP1_PP2_PP3_PP4_PP5-weights.png', dpi=300)
plt.close()

exit()
