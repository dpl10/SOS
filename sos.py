#!/usr/bin/env python3

import sys
import os
from phylo import Tree

tree_file = None
min_tax = 4
verbose_mode = False
debug_mode = False

myhelp = """
SOS: Sequestering Orthologous Subclades

A Python script to encode orthologous clades in a phylogenetic tree.

USAGE:	sos.py -t <tree_file> -m <#> [-v] > output_tsv_file

WHERE:
	-t	Input tree file.

	-m	Minimum taxa per orthologous set (default = 4).

	-v	Verbose mode: outputs a matrix if there are a single orthologous set (single
		column of ones) or if there are none (single column of zeros).

CITATION: Salinas, Sondervan, Tessler and Little 2022. SOS: Sequestering Orthologous 
Subclades. https://github.com/dpl10/SOS.

"""

for ia in range(len(sys.argv)):

	if sys.argv[ia] == "-t" and len(sys.argv) > (ia + 1):

		if os.path.exists(sys.argv[ia+1]):
			tree_file = sys.argv[ia+1]

		else:
			print(f"Could not open {sys.argv[ia+1]}.")
			exit()

	elif sys.argv[ia] == "-m" and len(sys.argv) > (ia + 1):
		
		try:
			min_tax = int(sys.argv[ia+1])

		except ValueError:
			print("Minimum taxa per set (-m) should be an integer.")
			exit()

		except:
			raise

		if min_tax <= 0:
			print("Minimum taxa per set (-m) should be an integer greater than zero.")
		

	elif sys.argv[ia] == "-v":
		verbose_mode = True

	elif sys.argv[ia] == "-d":
		debug_mode = True


if tree_file and min_tax > 0:
	print(f"{tree_file=}")
	print(f"{min_tax=}")
	print(f"{verbose_mode=}")

	tr = Tree(tree_file, debug=debug_mode)
	print(tr.tsv_table(min_tax, verbose=verbose_mode, debug=debug_mode))


else:
	print(myhelp)


exit()