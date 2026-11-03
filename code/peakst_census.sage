load("qwalk.sage")
load("peakST.sage")
load("cayley_graphs.sage")

import os, sys, tarfile, gzip, time

NMIN, NMAX = int(sys.argv[1]), int(sys.argv[2])
VT_DIR, OUT_DIR = sys.argv[3], sys.argv[4]

HEADER = ["graph6", "n", "k", "charpoly", "groups", "partner", "tau",
          "amount", "amount_approx", "pst"]


def vt_graphs(n):
    """Yield (graph6, k) for the connected vertex-transitive graphs on n
    vertices."""
    if n <= 9:
        for k in range(1, n):
            if (n * k) % 2 or (k == 1 and n > 2):
                continue
            for G in graphs.nauty_geng("%d -c -d%d -D%d" % (n, k, k)):
                if G.is_vertex_transitive():
                    yield G.graph6_string(), k
        return
    with tarfile.open(os.path.join(VT_DIR, "alltrans%d.tar" % n)) as tf:
        members = sorted((m for m in tf.getmembers() if m.name.endswith(".gz")),
                         key=lambda m: m.name)
        for m in members:
            k = int(os.path.basename(m.name).split("k=")[1].split(".")[0])
            if k == 0:
                continue
            for raw in gzip.GzipFile(fileobj=tf.extractfile(m)):
                g6 = raw.strip().decode("ascii")
                if g6 and Graph(g6).is_connected():
                    yield g6, k


def exact_idempotents(A):
    """{theta: E_theta} for the eigenvalues theta of A of degree <= 2 over Q,
    computed exactly (theta in QQ or in QuadraticField(D))."""
    out = {}
    for f, _m in A.charpoly().factor():
        if f.degree() == 1:
            F = QQ
        elif f.degree() == 2:
            b, c = f[1], f[0]
            F = QuadraticField(squarefree_part(ZZ(b*b - 4*c)), 'r')
        else:
            continue
        AF = A.change_ring(F)
        for theta in f.roots(F, multiplicities=False):
            B = (AF - theta).right_kernel_matrix()
            out[theta] = B.transpose() * (B * B.transpose()).inverse() * B
    return out


def common_field(xs):
    """The quadratic field among the parents of xs, or QQ."""
    for x in xs:
        if x.parent() is not QQ:
            return x.parent()
    return QQ


def fmt_quadratic(x):
    """a + b*sqrt(D) as a string Sage can read back."""
    F = x.parent()
    if F is QQ or x in QQ:
        return str(QQ(x))
    a, b = x.list()
    D = ZZ(F.gen()**2)
    if a == 0:
        return "%s*sqrt(%d)" % (b, D)
    return "%s %s %s*sqrt(%d)" % (a, "-" if b < 0 else "+", abs(b), D)


def peak_rows(g6, n, k, groups):
    G = Graph(g6)
    A = G.adjacency_matrix().change_ring(QQ)
    pre = _peakST_precompute(A)
    hits = []
    for v in range(1, n):
        res = hasPeakStateTransferPair(A, 0, v, _pre=pre)
        if res is not False:
            hits.append((v, res[1]))
    if not hits:
        return []
    E = exact_idempotents(A)
    cp = str(A.charpoly().factor())
    rows = []
    for v, tau in hits:
        support = [Et[0, v] for Et in E.values() if Et[0, v] != 0]
        F = common_field(support)
        s = sum((F(x) if F(x) >= 0 else -F(x)) for x in support)
        amount = s * s
        approx = peakStateTransferAmount(A, 0, v, tau, _pre=pre)
        if abs(float(amount) - approx) > 1e-8:
            raise RuntimeError("amount mismatch %s v=%d: %s vs %r"
                               % (g6, v, amount, approx))
        rows.append([g6, n, k, cp, groups, v, tau, fmt_quadratic(amount),
                     "%.10f" % approx, "yes" if amount == 1 else "no"])
    return rows


def group_labels(n):
    """{canonical g6 : '(n,i):Struct:ab;...'} for the connected Cayley graphs
    on n vertices."""
    names = {}
    for i in range(1, int(libgap.NumberSmallGroups(n)) + 1):
        H = libgap.SmallGroup(n, i)
        names[i] = "(%d,%d):%s:%s" % (n, i,
                                      str(H.StructureDescription()).replace(" ", ""),
                                      "ab" if H.IsAbelian() else "nonab")
    return {g6: ";".join(names[i] for i in idx)
            for g6, idx in cayley_graphs_of_order(n).items()}


summary = [["n", "vt_graphs", "cayley_graphs", "peakst_graphs",
            "peakst_cayley", "peakst_noncayley", "pst_graphs"]]
for n in range(NMIN, NMAX + 1):
    t0 = time.time()
    labels = group_labels(n)
    n_vt = n_cay = 0
    found = []
    seen_cayley = set()
    for g6, k in vt_graphs(n):
        n_vt += 1
        canon = Graph(g6).canonical_label().graph6_string()
        groups = labels.get(canon, "-")
        if groups != "-":
            n_cay += 1
            seen_cayley.add(canon)
        rows = peak_rows(g6, n, k, groups)
        if rows:
            found.append((g6, rows))
    assert seen_cayley == set(labels), "n=%d: Cayley graphs missing from VT list" % n

    if found:
        d = os.path.join(OUT_DIR, "Order%02d" % n)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "Ord%d-g6.txt" % n), "w") as f:
            for g6, _ in found:
                f.write(g6 + "\n")
        with open(os.path.join(d, "Ord%d-peakST.tsv" % n), "w") as f:
            f.write("\t".join(HEADER) + "\n")
            for _, rows in found:
                for r in rows:
                    f.write("\t".join(str(x) for x in r) + "\n")
    n_peak = len(found)
    n_peak_cay = sum(1 for _, rows in found if rows[0][4] != "-")
    n_pst = sum(1 for _, rows in found if any(r[-1] == "yes" for r in rows))
    summary.append([n, n_vt, n_cay, n_peak, n_peak_cay, n_peak - n_peak_cay, n_pst])
    print("n=%2d: %4d VT, %4d Cayley, %4d peak ST (%d Cayley), %d PST  [%.0fs]"
          % (n, n_vt, n_cay, n_peak, n_peak_cay, n_pst, time.time() - t0))
    sys.stdout.flush()

with open(os.path.join(OUT_DIR, "summary.tsv"), "w") as f:
    for r in summary:
        f.write("\t".join(str(x) for x in r) + "\n")
