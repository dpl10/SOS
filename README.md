# SOS
SOS: Sequestering Orthologous Subclades

# Algorithm description

**Input:** A phylogenetic tree stored in Newick or TNT formats. Polytomies are treated as such. 

**Output:** A tsv table in which rows are tree terminal labels and columns are encoded orthologous clades.

1. The tree is unrooted and stored in memory as an adjacency matrix.

2. Internal nodes are tested for the ortholog condition. This condition assumes that a valid ortholog node:

	a. Has two or more descendant taxa.

	b. All descendant taxa conform to a set of mutually exclusive clades.

3. The unrooted tree is traversed postorder, but traversal is stopped to retrieve the most inclusive orthologous clades.

4. An output matrix is initialized, and a column is encoded for each orthologous clade from step 3.

# Requirements

- Python 3 interpreter.
- Numpy.
- Scipy.
- Pytest (for testing).

# Installation


# Usage

`sos.py -t <tree_file> -m <#> [-v] > output_tsv_file`

__where__:

| Option | Explanation |
|:---|:---|
| `-t` | Input tree file in Newick or TNT format. The file should contain a single tree. |
| `-m` | Minimum taxa per orthologous set (default = 4) |
| `-v` | Activates verbose mode. |

By default the script does not produce any output under non-informative circunstances, 
such as the tree having a single species, a single representative per species or 
not containing any orthologous set. In verbose mode, it will output a non-informative 
matrix for each one of these cases.

# License



# Contact

Nelson R. Salinas  
nrsalinas@gmail.com  

Damon Little  
dlittle@nybg.org

New York Botanical Garden
