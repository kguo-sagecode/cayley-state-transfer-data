import itertools


def _group_data(o, i):
    """Elements of SmallGroup(o,i) (identity first), its multiplication table
    on indices, and the classes {g, g^-1} of non-identity elements."""
    H = libgap.SmallGroup(o, i)
    els = list(H.AsList())
    one = H.Identity()
    els.sort(key=lambda g: g != one)
    idx = {g: j for j, g in enumerate(els)}
    mult = [[idx[a * b] for b in els] for a in els]
    inv = [idx[g.Inverse()] for g in els]
    classes, seen = [], set([0])
    for j in range(1, o):
        if j in seen:
            continue
        cls = sorted(set([j, inv[j]]))
        seen.update(cls)
        classes.append(cls)
    return H, els, mult, classes


def cayley_graphs(o, i, connected=True, sizes=None):
    """{canonical graph6 : connection set} for the Cayley graphs of
    SmallGroup(o, i), one entry per isomorphism class.

    The connection set is returned as a list of GAP elements.  With
    connected=True (default) only connection sets that generate the group are
    kept.  sizes, if given, restricts the valency |S|."""
    H, els, mult, classes = _group_data(o, i)
    out = {}
    for r in range(1, len(classes) + 1):
        for combo in itertools.combinations(classes, r):
            S = [s for cls in combo for s in cls]
            if sizes is not None and len(S) not in sizes:
                continue
            X = Graph([range(o), [(a, mult[a][s]) for a in range(o) for s in S]],
                      format='vertices_and_edges')
            if connected and not X.is_connected():
                continue
            g6 = X.canonical_label().graph6_string()
            if g6 not in out:
                out[g6] = [els[s] for s in S]
    return out


def cayley_graphs_of_order(n, connected=True):
    """{canonical graph6 : [i, ...]}: every Cayley graph on n vertices up to
    isomorphism, with the indices i such that it is a Cayley graph of
    SmallGroup(n, i)."""
    out = {}
    for i in range(1, int(libgap.NumberSmallGroups(n)) + 1):
        for g6 in cayley_graphs(n, i, connected=connected):
            out.setdefault(g6, []).append(i)
    return out


def is_cayley_graph_of(X, o, i):
    """True iff the graph X is a Cayley graph of SmallGroup(o, i).
    (Exhaustive over connection sets of size equal to the valency of X.)"""
    if X.order() != o or not X.is_regular():
        return False
    target = X.canonical_label().graph6_string()
    forms = cayley_graphs(o, i, connected=X.is_connected(),
                          sizes=[X.degree(X.vertices()[0])])
    return target in forms
