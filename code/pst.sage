def _strongly_cospectral_pair(sp, u, v, eps=1e-6):
    r"""u and v are strongly cospectral iff, for every spectral idempotent E_r,
    (E_r)_{u,u} = (E_r)_{v,v} and |(E_r)_{u,v}| = (E_r)_{u,u}."""
    for (theta, E) in sp:
        Eu  = real(CC(E[u, u]))
        Ev  = real(CC(E[v, v]))
        Euv = real(CC(E[u, v]))
        if abs(Eu - Ev) > eps:
            return False
        if abs(abs(Euv) - Eu) > eps:
            return False
    return True


def hasPerfectStateTransferPair(M, u, v, eps=1e-6, _pre=None):
    r"""Determine whether there is perfect state transfer between distinct
    vertices u and v with respect to the symmetric matrix M.

    Returns False, or (True, tau0) where tau0 is the minimum positive time
    at which PST occurs.
    """
    if _pre is None:
        _pre = _peakST_precompute(M)
    if not _strongly_cospectral_pair(_pre['sp'], u, v, eps=eps):
        return False
    return hasPeakStateTransferPair(M, u, v, eps=eps, _pre=_pre)


def hasPerfectStateTransfer(G, M=None, vt=False):
    r"""Check whether graph G admits perfect state transfer between any pair
    of distinct vertices, with respect to M (default: adjacency matrix).

    Returns False, or (True, tau0) -- the smallest positive PST time across
    all pairs that admit it.

    If vt=True, treat G as vertex-transitive: only test pairs (0, v) and
    return on the first hit.  By transitivity, PST anywhere implies PST
    from vertex 0 to some v with the same time tau, so this is sound.
    """
    if M is None:
        M = G.adjacency_matrix()
    M = Matrix(QQ, M)
    n = M.nrows()
    pre = _peakST_precompute(M)
    if vt:
        for v in range(1, n):
            if not _strongly_cospectral_pair(pre['sp'], 0, v):
                continue
            res = hasPeakStateTransferPair(M, 0, v, _pre=pre)
            if res is False:
                continue
            return (True, res[1])
        return False
    best = None
    for u in range(n):
        for v in range(u + 1, n):
            if not _strongly_cospectral_pair(pre['sp'], u, v):
                continue
            res = hasPeakStateTransferPair(M, u, v, _pre=pre)
            if res is False:
                continue
            _, tau = res
            if best is None or tau < best:
                best = tau
    if best is None:
        return False
    return (True, best)


def perfectStateTransferPairs(G, M=None):
    r"""Return a list of triples (u, v, tau0) for every pair u < v of G that
    admits perfect state transfer with respect to M (default: adjacency)."""
    if M is None:
        M = G.adjacency_matrix()
    M = Matrix(QQ, M)
    n = M.nrows()
    pre = _peakST_precompute(M)
    out = []
    for u in range(n):
        for v in range(u + 1, n):
            if not _strongly_cospectral_pair(pre['sp'], u, v):
                continue
            res = hasPeakStateTransferPair(M, u, v, _pre=pre)
            if res is not False:
                out.append((u, v, res[1]))
    return out
