load("cayley_graphs.sage")

import os, sys, ast

DATA = sys.argv[1]
NMIN = int(sys.argv[2]) if len(sys.argv) > 2 else 6
NMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 23

counts = {}
with open(os.path.join(DATA, "NrOfGraphs.txt")) as f:
    for line in f:
        if ":" in line and line.split(":")[0].strip().isdigit():
            a, b = line.split(":")
            counts[int(a)] = int(b)


def canon(g6):
    return Graph(g6).canonical_label().graph6_string()


ok = True
for n in range(NMIN, NMAX + 1):
    d = os.path.join(DATA, "Order%02d" % n)
    with open(os.path.join(d, "Ord%d-g6.txt" % n)) as f:
        g6s = [l.strip() for l in f if l.strip()]
    ref = cayley_graphs_of_order(n)
    theirs = [canon(s) for s in g6s]
    problems = []
    if len(set(theirs)) != len(theirs):
        problems.append("repeated graphs")
    if set(theirs) != set(ref):
        problems.append("graph sets differ (%d listed, %d computed)"
                        % (len(set(theirs)), len(ref)))
    if counts.get(n) != len(g6s):
        problems.append("NrOfGraphs.txt says %s, file has %d"
                        % (counts.get(n), len(g6s)))
    gg = os.path.join(d, "Ord%d-GraphGroups.txt" % n)
    if os.path.exists(gg):
        with open(gg) as f:
            pairs = [ast.literal_eval(l.strip()) for l in f if l.strip()]
        if [p[0] for p in pairs] != g6s:
            problems.append("GraphGroups.txt and g6.txt list different graphs")
        for s, idx in pairs:
            if sorted(idx) != sorted(ref.get(canon(s), [])):
                problems.append("groups of %s: listed %s, computed %s"
                                % (s, sorted(idx), sorted(ref.get(canon(s), []))))
    elif int(libgap.NumberSmallGroups(n)) != 1:
        problems.append("no GraphGroups.txt, but %d groups of this order"
                        % libgap.NumberSmallGroups(n))
    print("n=%2d: %5d graphs  %s" % (n, len(g6s), "; ".join(problems) or "OK"))
    sys.stdout.flush()
    ok = ok and not problems

print("ALL OK" if ok else "PROBLEMS FOUND")
