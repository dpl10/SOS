# SOS
SOS: Sequestering Orthologous Subclades

# Algorithm description

**Input:** A phylogenetic tree stored in TNT or Newick formats. Polytomies are treated as such. 

**Output:** A tsv table in which rows are tree terminal labels and columns are encoded orthologous clades.

1. The tree is unrooted and stored in memory as an adjacency matrix.

2. Internal nodes are tested for the ortholog condition. This condition assumes that a valid ortholog node:

	a. Has two or more descendant taxa.

	b. All descendant taxa conform to a set of mutually exclusive clades.

3. The unrooted tree is traversed postorder, but traversal is stopped to retrieve the most inclusive orthologous clades.

4. An output matrix is initialized, and a column is encoded for each orthologous clade from step 3.
