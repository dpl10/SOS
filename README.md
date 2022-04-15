# SOS
SOS: Sequestering Orthologous Subclades

# Algorithm description

1. Trees are read from a tree file (TNT or Newick formats).

2. Then they are unrooted.

3. All internal nodes are tested for the ortholog condition. Therefore, for a each node:

	a. All descendant taxa conform to a set of mutually exclusive clades.

	b. Node has two or more descendant taxa.

4. The most inclusive orthologous sets are tagged.
