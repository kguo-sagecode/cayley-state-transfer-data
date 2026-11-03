For n (currently <30), the folder "Order n" contains two text files.

Case 1: There are more than one group of order n.
The two files are "Ordn-g6.txt" and "Ordn-GraphGroups.txt". The former contains the graph6- strings of all non-isomorphic Cayley graphs on n vertices, one on each line. The latter contains the same graphs, but here, each line is a list [G6,Gr], where G6 is the graph6-string and Gr is a list of indices of groups in the gap group library. So if i is in the list Gr, then G6 is a Cayley graph for the group gap.SmallGroup(n,i). The graph on line k should be the same in both files.

Case 2: There is exactly one group of order n.
This group is the cyclic group. The two files are "Ordn-g6.txt" and "Ordn-GraphSet.txt". The former is as in Case 1. The latter again has the same graphs, but here, each line is a list, [G6,S] where G6 is the graph6-string and S is a connection set for the the group C_n = <f>, that yields the corresponding graph. It can be read into Sage with the following code:

f = open('Ord%s-GraphSet.txt'%n,'r')
strings = f.readlines()
f.close()

G = AbelianGroup([n])
(f,) = G.gens()

s = []
for st in strings:
    el = sage_eval(st.strip(), locals={'f': f})
    s.append(el)
