# SOS
SOS: Sequestering Orthologous Subclades.

A Python utility to encode sets of ortholog sequences in a phylogenetic tree, assuming a simple ortholog membership rule.

Polytomies are treated as such.

#### Procedure description

**Input:** A text file containing a single phylogenetic tree, following the Newick or TNT standard. 

**Output:** A tsv table in which rows are tree terminal labels and columns are encoded orthologous clades.

1. The tree is unrooted and stored in memory as an sparse adjacency matrix.

2. All internal nodes are tested for the ortholog condition. This condition assumes that a valid ortholog node:

	a. Has two or more descendant taxa.

	b. All descendant taxa conform to a set of mutually exclusive clades.

3. The unrooted tree is traversed postorderly, and the most inclusive orthologous clades are retrieved.

4. An output matrix is initialized, and a column is encoded for each orthologous clade from step 3.

# Requirements

- Python 3 interpreter.
- Numpy.
- Scipy.
- Pytest (for testing).

# Installation

This software is distributed through [Github](https://github.com/dpl10/SOS): you can either download the compressed zip file of the repository or clone it:

`git clone https://github.com/dpl10/SOS`


# Usage

`sos.py -t <tree_file> [-m <#>] [-v] > <output_file>`

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

Copyright (c) 2022 Nelson R. Salinas & Damon Little.

SOS is available under the MIT License. See LICENSE for more information.

# Contact

Nelson R. Salinas  
nrsalinas@gmail.com  

Damon Little  
dlittle@nybg.org

New York Botanical Garden
