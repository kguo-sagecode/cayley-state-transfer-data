def _peakST_precompute(M):
    r"""
    Precompute per-matrix data used by hasPeakStateTransferPair:
      - numerical spectral decomposition (for signs of idempotent entries);
      - factorization of the characteristic polynomial over QQ (for exact
        eigenvalues as algebraic integers);
      - a matching from each numerical eigenvalue to its exact AA root and
        factor index.
    """
    M = Matrix(QQ, M)
    sp = spcdec(Matrix(CC, M))
    chi = M.charpoly()
    factored = list(chi.factor())
    factor_data = [(f, f.roots(AA, multiplicities=False)) for (f, _m) in factored]
    numerical_exact = []
    for (theta_num, E) in sp:
        th = float(real(CC(theta_num)))
        matched = None
        for fi, (f, rts) in enumerate(factor_data):
            for r in rts:
                if abs(float(r) - th) < 1e-4:
                    matched = (fi, f, r)
                    break
            if matched:
                break
        numerical_exact.append((th, E, matched))
    return {'M': M, 'sp': sp, 'factor_data': factor_data,
            'numerical_exact': numerical_exact}


def hasPeakStateTransferPair(M, u, v, eps=1e-6, _pre=None):
    r"""
    Check peak state transfer between distinct vertices u and v with respect
    to the symmetric matrix M (rational entries, e.g. adjacency or Laplacian).

    Implements the spectral characterization of arXiv:2505.11986.

    Returns False if there is no peak state transfer; otherwise (True, tau0)
    where tau0 = pi/(g*sqrt(D)) is the minimum positive time of peak state
    transfer.
    """
    if _pre is None:
        _pre = _peakST_precompute(M)
    numerical_exact = _pre['numerical_exact']
    factor_data     = _pre['factor_data']

    support = []
    for (th, E, matched) in numerical_exact:
        entry = real(CC(E[v, u]))
        if abs(entry) <= eps:
            continue
        if matched is None:
            return False
        fi, f, r = matched
        support.append((fi, r, +1 if entry > 0 else -1, th))

    if not support:
        return False
    if not any(s[2] == +1 for s in support):
        return False
    if not any(s[2] == -1 for s in support):
        return False

    included = set(s[0] for s in support)
    for fi in included:
        f, rts = factor_data[fi]
        for r in rts:
            rn = float(r)
            if not any(abs(sr[3] - rn) < 1e-4 for sr in support):
                return False

    for fi in included:
        f, _ = factor_data[fi]
        if f.degree() > 2:
            return False

    D = None
    for fi in included:
        f, _ = factor_data[fi]
        if f.degree() == 2:
            cs = f.coefficients(sparse=False)
            p_coef = ZZ(cs[1]); q_coef = ZZ(cs[0])
            disc = p_coef*p_coef - 4*q_coef
            if disc <= 0:
                return False
            D_this = squarefree_part(disc)
            if D is None:
                D = ZZ(D_this)
            elif ZZ(D_this) != D:
                return False
    if D is None:
        D = ZZ(1)

    if D == 1:
        a_common = ZZ(0)
    else:
        a_common = None
        for fi in included:
            f, _ = factor_data[fi]
            cs = f.coefficients(sparse=False)
            if f.degree() == 2:
                a_this = -ZZ(cs[1])
            else:
                a_this = 2 * (-ZZ(cs[0]))
            if a_common is None:
                a_common = a_this
            elif a_this != a_common:
                return False

    sqrtD_AA = AA(D).sqrt()
    support_b = []
    for (fi, r_AA, sign, _) in support:
        val = (2*r_AA - a_common) / sqrtD_AA
        try:
            val_Q = QQ(val)
        except (TypeError, ValueError):
            return False
        if val_Q.denominator() != 1:
            return False
        support_b.append((sign, ZZ(val_Q)))

    a_par = int(a_common) % 2
    for (_, b) in support_b:
        if (b - a_par) % 2 != 0:
            return False

    s_idx = next(i for i, (sg, _) in enumerate(support_b) if sg == 1)
    b_s = support_b[s_idx][1]

    c_values = [(sg, ZZ((b_s - b) / 2)) for (sg, b) in support_b]

    nonzero_cs = [abs(c) for (_, c) in c_values if c != 0]
    if not nonzero_cs:
        return False
    g = ZZ(nonzero_cs[0])
    for c in nonzero_cs[1:]:
        g = gcd(g, c)

    for (sg, c) in c_values:
        if c % g != 0:
            return False
        quot = c // g
        if sg == +1 and quot % 2 != 0:
            return False
        if sg == -1 and quot % 2 == 0:
            return False

    tau0 = pi / (g * sqrt(D))
    return (True, tau0)


def peakStateTransferAmount(M, u, v, tau, _pre=None):
    r"""
    Numerically evaluate |U(tau)_{v,u}|^2 where U(tau) = exp(i*tau*M) and
    M is symmetric.  Uses the cached spectral decomposition pre['sp']
    (a list of (theta, E_theta) pairs of float eigenvalues and idempotent
    matrices).  Returns a Python float in [0, 1].
    """
    if _pre is None:
        _pre = _peakST_precompute(M)
    tau_f = float(tau)
    amp = complex(0)
    for (theta_num, E) in _pre['sp']:
        th = float(real(CC(theta_num)))
        amp += complex(CC(E[v, u])) * complex(cos(th * tau_f), sin(th * tau_f))
    return float(amp.real * amp.real + amp.imag * amp.imag)


def hasPeakStateTransfer(G, M=None, vt=False):
    r"""
    Check whether graph G admits peak state transfer between any pair of
    distinct vertices.  By default the adjacency matrix is used; pass M
    explicitly to use a different symmetric matrix.

    Returns False if no pair admits peak state transfer, otherwise
    (True, tau0, amount) where tau0 is the smallest positive peak-ST time
    across all qualifying pairs and `amount` is |U(tau0)_{v,u}|^2 for that
    minimizing pair.

    If vt=True, treat G as vertex-transitive: only test pairs (0, v).
    Different (Aut)_0-orbits on V \ {0} can have different amplitudes, so
    we still scan all v and return the (tau, amount) with maximum amount.
    """
    if M is None:
        M = G.adjacency_matrix()
    M = Matrix(QQ, M)
    n = M.nrows()
    pre = _peakST_precompute(M)
    if vt:
        best = None
        for v in range(1, n):
            res = hasPeakStateTransferPair(M, 0, v, _pre=pre)
            if res is False:
                continue
            _, tau = res
            amount = peakStateTransferAmount(M, 0, v, tau, _pre=pre)
            if best is None or amount > best[0]:
                best = (amount, tau)
        if best is None:
            return False
        amount, tau = best
        return (True, tau, amount)
    best = None
    best_pair = None
    for u in range(n):
        for v in range(u + 1, n):
            res = hasPeakStateTransferPair(M, u, v, _pre=pre)
            if res is False:
                continue
            _, tau = res
            if best is None or tau < best:
                best = tau
                best_pair = (u, v)
    if best is None:
        return False
    u, v = best_pair
    amount = peakStateTransferAmount(M, u, v, best, _pre=pre)
    return (True, best, amount)


def peakStateTransferPairs(G, M=None):
    r"""
    Return a list of (u, v, tau0, amount) for every pair u < v of vertices
    of G that admits peak state transfer with respect to M (default:
    adjacency).  `amount` = |U(tau0)_{v,u}|^2.
    """
    if M is None:
        M = G.adjacency_matrix()
    M = Matrix(QQ, M)
    n = M.nrows()
    pre = _peakST_precompute(M)
    out = []
    for u in range(n):
        for v in range(u + 1, n):
            res = hasPeakStateTransferPair(M, u, v, _pre=pre)
            if res is False:
                continue
            _, tau = res
            amount = peakStateTransferAmount(M, u, v, tau, _pre=pre)
            out.append((u, v, tau, amount))
    return out
