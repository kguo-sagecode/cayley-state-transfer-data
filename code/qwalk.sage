import scipy
from scipy import linalg
import numpy


def partn( rls, eps):
    dc={rls[0]:[0]}
    keys=[rls[0]]
    for i in [1..len(rls)-1]:
        x = rls[i]
        key_exists = False
        for j in [0..len(keys)-1]:
            if abs( x-keys[j] ) < eps:
                key_exists = True
                dc[keys[j]] += [i]
                break
        if not key_exists:
            dc[x] = [i]
            keys += [x]
    return dc

def spcdec( A):
    eps = 1e-6
    n = A.nrows()
    la, v = linalg.eigh(A)
    dc = partn( la, eps)
    evals = dc.keys()
    Es = []
    for ev in dc.keys():
        E  = Matrix(CC,n,n)
        for ind in dc[ev]:
            col = numpy.asarray(v[:,ind]).flatten()
            E += Matrix(CC, numpy.outer(col, col))
        Es += [(ev, E)]
    return Es

def partn(lst, eps = 0):
    dc = {lst[0]:[0]}
    keys = [lst[0]]
    for it in [1..len(lst) - 1]:
        new_key = lst[it]
        key_exists = False
        for key in keys:
            if abs(new_key - key) <= eps:
                key_exists = True
                dc[key].append(it)
                break
        if not key_exists:
                dc[new_key] = [it]
                keys.append(new_key)
    return dc

def spd_herm(A):
    eps = 1e-6
    n = A.nrows()
    evs, evecs = linalg.eigh(A.change_ring(CC))
    dc = partn(evs, eps)
    Es = []
    for ev in dc.keys():
        emtx = Matrix(CC, evecs[:, dc[ev]])
        eproj = emtx * emtx.conjugate_transpose()
        Es.append((ev, eproj, emtx))
    return Es

def spd(M):
    A = 1/2*(M + M.conjugate_transpose())
    B = i/2*(M.conjugate_transpose() - M)
    Es_A = spd_herm(A.change_ring(CC))
    Es_M = []
    for it in Es_A:
        C = it[2].conjugate_transpose()*B*it[2]
        Es_C = spd_herm(C.change_ring(CC))
        for j in Es_C:
            Es_M.append((it[0] + i*j[0], it[2]*j[1]*it[2].conjugate_transpose()))
    return Es_M


def entropy( pv):
    return -sum([ vv* log(vv) for vv in pv if vv != 0])

def matrix_entropy( M):
    return sum([entropy(vv) for vv in M])

def transmat( sp):
    return lambda t: sum([exp(I*t*it[0])*it[1] for it in sp])

def transmat_col( sp, a):
    return lambda t: sum([exp(I*t*it[0])*it[1][a] for it in sp])

def transmat_entry( sp, a, b):
    return lambda t: sum([exp(I*t*it[0])*it[1][a][b] for it in sp])

def stateprob( sp):
    return lambda t: sum([exp(I*t*it[0])*it[1] for it in sp]).apply_map(norm)

def stateprob_row( sp, a):
    return lambda t: list(map( norm, sum([exp(I*t*it[0])*it[1][a] for it in sp])))

def stateprob_entry( sp, a,b):
    return lambda t: abs(sum([exp(I*t*it[0])*it[1][a][b] for it in sp]))^2

def entropy_row( sp, a):
    nn = sp[0][1].nrows()
    return lambda t: entropy( list(map( norm, sum([exp(I*t*it[0])*it[1][a] for it in sp]))))/log(float(nn))

def avg( sp):
    return sum([it[1].elementwise_product(it[1]) for it in sp])

def avgcx( sp):
    return sum([Matrix(it[1]).elementwise_product(Matrix(it[1]).conjugate()) for it in sp])

def abs_sum( sp):
    return sum( it[1].apply_map(abs) for it in sp)

def parplot(sp,a,b,t0,t1, pltpts=None):
    if pltpts==None: pltpts = max(1000,10*t1)
    t=var('t')
    return parametric_plot( (transmat_entry(sp,a,b)(t).real(), transmat_entry(sp,a,b)(t).imag()), (t,t0,t1),\
        plot_points=pltpts)


def temp(A):
    eps = 1e-6
    n = A.nrows()
    _,evs,_ = linalg.svd(A.change_ring(CC))
    L = list(map(rd,evs))
    return L.count(0)

def rd(x):
    y = CC(x)
    return round(real(y),6) + i*round(imag(y),6)

def rd_mat(A):
    return A.apply_map(rd)

def coeffmtx(X):
	n = X.num_verts()
	cm = Matrix(n,n)
	polys = []
	for it in range(n):
		now = Graph(X.am())
		now.delete_vertex(it)
		curr = now.am().charpoly()
		cm.set_row(it,list(curr))
		polys.append(curr)
	return [cm,polys]

def hasSimpleEigs(m):
	now = m.fcp()
	simp = true
	for it in now:
		simp = simp and (it[1]==1)
	return simp

def hasIntegerEigs(m):
	now = m.fcp()
	allint = true
	for it in now:
		allint = allint and (it[0].degree() == 1)
	return allint

def parallelVxs(X,u,v):
	Y=copy(X)
	Y.delete_vertices([u,v])
	g = gcd(X.am().charpoly(), Y.am().charpoly())
	f = X.am().charpoly()/g
	now = list(f.factor())
	degs = list(set([it[1] for it in now]))
	return [(len(degs) ==1 and degs[0] == 1), X.am().fcp(), Y.am().fcp()]

def parallelVxsL(X,u,v):
	Y=copy(X)
	Y.delete_vertices([u,v])
	g = gcd(X.laplacian_matrix().charpoly(), Y.laplacian_matrix().charpoly())
	f = X.laplacian_matrix().charpoly()/g
	now = list(f.factor())
	degs = list(set([it[1] for it in now]))
	return [(len(degs) ==1 and degs[0] == 1), X.laplacian_matrix().fcp(), Y.laplacian_matrix().fcp()]

def periodicVx(X,v):
	Y=copy(X)
	Y.delete_vertices([v])
	g = gcd(X.am().charpoly(), Y.am().charpoly())
	f = ((X.am().charpoly())^2 - (Y.am().charpoly())^2)/g^2
	now = list(f.factor())
	degs = set([it[0].degree() for it in now])
	return (degs == set([1,2]))

def periodicVxL(X,v):
	Y=copy(X)
	Y.delete_vertices([v])
	g = gcd(X.laplacian_matrix().charpoly(), Y.laplacian_matrix().charpoly())
	f = ((X.laplacian_matrix().charpoly())^2 - (Y.laplacian_matrix().charpoly())^2)/g^2
	now = list(f.factor())
	degs = set([it[0].degree() for it in now])
	return (degs == set([1,2]))

def findPeriodicVxs(X):
	vxs = X.vertices(sort=True)
	fd = []
	for it in vxs:
		if periodicVxL(X,it):
			fd.append(it)
	return fd
