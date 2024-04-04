import moments
import numpy as np

def pulse_migration(num_pops, pop0, pop1, f):
    """
    Pulse migration from pop0 into pop1 with proportion f.
    """
    # create the matrix that would build a new population via admixture
    A = moments.LD.Matrices.admix_ld(num_pops, pop0, pop1, f)
    # the last population created replaces pop1
    mom_from = moments.LD.Util.moment_names(num_pops + 1)[0]
    mom_to = moments.LD.Util.moment_names(num_pops)[0]
    P = np.zeros((len(mom_to), len(mom_from)))
    for i, m in enumerate(mom_to):
        l = m.split("_")
        for k in range(1, len(l)):
            if l[k] == str(pop1):
                l[k] = str(num_pops)
        m_from = "_".join(l)
        m_from = moments.LD.Util.map_moment(m_from)
        j = mom_from.index(m_from)
        P[i, j] = 1
    return P.dot(A)
