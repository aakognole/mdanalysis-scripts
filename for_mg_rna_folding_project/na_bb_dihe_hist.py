import MDAnalysis as mda
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import MDAnalysis.analysis.rdf as rdf

def nucleic_acid_backbone_dihedrals_histogram(u):
    print "returning histogram as matplotlib object..."
    u = u
    output = filename
    bb = u.select_atoms("nucleic and (name C4' or name C3' or name O3' or name P or name O5' or name C5')")
    D = []
    j=0
    while j <= len(bb.ids)-4:
        first = u.select_atoms("bynum %d" % (bb.ids[j]+1))
        if first.names[0] == "C4'":
            epsilon = [bb.ids[j], bb.ids[j+1], bb.ids[j+2], bb.ids[j+3]]
            zeta = [bb.ids[j+1], bb.ids[j+2], bb.ids[j+3], bb.ids[j+4]]
            alpha = [bb.ids[j+2], bb.ids[j+3], bb.ids[j+4], bb.ids[j+5]]
            beta = [bb.ids[j+3], bb.ids[j+4], bb.ids[j+5], bb.ids[j+6]]
            gamma = [bb.ids[j+4], bb.ids[j+5], bb.ids[j+6], bb.ids[j+7]]
            for ts in u.trajectory:
                alp = mda.core.topologyobjects.Dihedral(alpha, u, type=None)
                bet = mda.core.topologyobjects.Dihedral(beta, u, type=None)
                gam = mda.core.topologyobjects.Dihedral(gamma, u, type=None)
                eps = mda.core.topologyobjects.Dihedral(epsilon, u, type=None)
                zet = mda.core.topologyobjects.Dihedral(zeta, u, type=None)    
                a = alp.value()
                b = bet.value()
                g = gam.value()
                e = eps.value()
                z = zet.value()
                if a < 0:
                    a = 360 + a
                if b < 0:
                    b = 360 + b
                if g < 0:
                    g = 360 + g
                if e < 0:
                    e = 360 + e
                if z < 0:
                    z = 360 + z
                D.append((a, b, g, e, z))
        j += 1
    D = np.array(D)
    histogram = plt.figure(figsize=(9,15))
    plt.subplots_adjust(hspace=0.5)
    plt.subplot(511)
    plt.title('alpha')
    plt.xlim(0,360)
    plt.hist(D[:,0], bins=180, range=(0,360))
    plt.subplot(512)
    plt.title('beta')
    plt.xlim(0,360)
    plt.hist(D[:,1], bins=180, range=(0,360))
    plt.subplot(513)
    plt.title('gamma')
    plt.xlim(0,360)
    plt.hist(D[:,2], bins=180, range=(0,360))
    plt.subplot(514)
    plt.title('epsilon')
    plt.xlim(0,360)
    plt.hist(D[:,3], bins=180, range=(0,360))
    plt.subplot(515)
    plt.title('zeta')
    plt.xlim(0,360)
    plt.hist(D[:,4], bins=180, range=(0,360))
    return histogram
