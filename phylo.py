import re
import warnings
from functools import reduce
from itertools import combinations
import numpy as np
from scipy.sparse import csr_matrix, find
from time import time

class Tree:
	"""Simple class to manipulate unrooted phylogenetic trees."""
	def __init__(self, tree_file: str):
		"""Instantiates Tree class from a Newick tree file. Trees could have 
		polytomies, branch lengths are discarded."""
		self.list = []
		self.lengths = {}
		self.node_count = 0
		self.labels = {}
		self.taxa = {}
		self.adj_table = None
		self.results = None
		self.tree_file = tree_file

		#########################################
		# Test on tree files with superflous parentheses

		with open(tree_file , "r") as fh:
			for line in fh:
				line = line.strip()

				if re.search(r"#nexus", line, re.IGNORECASE):
					raise ValueError("Input file should conform to either TNT or Newick format.")

				if line.startswith("("):

					if re.search(',', line): # Newick
						line = re.sub(r'\)', ',)', line)
						line = re.sub(r'\),\(', ') (', line)
						line = re.sub(r',', ' ', line)
						line = re.sub(r'\):', ')=', line)
						line = re.sub(r':', '=', line)
						line = re.sub(r'\s+', ' ', line)

					line = re.sub(r"\s+\)", ")", line)
					node_pointer = -1
					label = ""
					in_br_len = False
					br_len = ""

					for char in line:
					
						if in_br_len:
							if char == " " or char == ")":
								self.lengths[node_pointer] = float(br_len)
								br_len = ""
								in_br_len = False

							elif char == "[" or char == "]":
								continue

							else:
								br_len += char

						if not in_br_len:

							if char == '(':
								pa = node_pointer
								node_pointer = self.node_count
								self.node_count += 1
								self.list.append([pa, node_pointer])
								
							elif char == ')':
								if len(label) > 0:
									self.labels[node_pointer] = label
									label = ""
								node_pointer = self.get_parent(node_pointer)
								
							elif char == " ":
								if len(label) > 0:
									self.labels[node_pointer] = label
									label = ""

								pa = self.get_parent(node_pointer)
								node_pointer = self.node_count
								self.node_count += 1
								self.list.append([pa, node_pointer])

							elif char == "=":
								if len(label) > 0:
									self.labels[node_pointer] = label
									label = ""
								in_br_len = True

							elif char == ';':
								break

							else:
								label += char
		
		for node in self.labels:
			self.taxa[node] = self.labels[node].split('#')[0]

		#"""
		for pa,n in self.list:
			if n in self.labels:
				print(pa, n, self.labels[n])
			else:
				print(pa, n)
		#"""
		
		root_edges_idx = []
		root_descendants = []

		for idx, edge in enumerate(self.list):
			if edge[0] == -1:
				root_descendants.append(edge[1])
				root_edges_idx.append(idx)
		
		for i,d in combinations(root_descendants, 2):
			self.list.append([i, d])

		self.list = [x for i,x in enumerate(self.list) if not i in root_edges_idx]


		"""
		Result table. Index meaning:
		First dimension = main node
		Second dimension = excluded node
		
		Values
		0 = Not computed yet
		1 = Negative
		2 = Positive
		"""
		rows = [x[0] for x in self.list] + [x[1] for x in self.list]
		cols = [x[1] for x in self.list] + [x[0] for x in self.list]
		vals = [1 for x in rows]
		vzeros = [0 for x in rows]
		self.adj_table = csr_matrix((vals, (rows, cols)), 
			shape=(self.node_count, self.node_count), dtype=np.int8)
		self.results = csr_matrix((vzeros, (rows, cols)), 
			shape=(self.node_count, self.node_count), dtype=np.int8)
		#self.adj_table = self.adj_table.tolil()
		
		self.edge_coors_x, self.edge_coors_y = np.where(self.adj_table.toarray() > 0)
		"""
		# Coordinates of internal edges in adjacency table
		ba = np.isin(self.edge_coors_x, list(self.taxa.keys()))
		bb = np.isin(self.edge_coors_y, list(self.taxa.keys()))
		self.edge_coors_x = self.edge_coors_x[~(ba | bb)]
		self.edge_coors_y = self.edge_coors_y[~(ba | bb)]
		"""


	def get_parent(self, node: int) -> int:
		""""To be used exclusively with the preliminary linked list during class
		instantiation."""
		for no, des in self.list:
			if des == node:
				return no


	def leaves_from_node(self, node : int, excluded: int) -> list:
		"""
		Returns the list of leaves (ints) that are descendant of a given node. 
		It is required to select a neighbor (`excluded`) of the latter node to 
		orientate (root) the operation.
		"""
		th = self.adj_table[node].toarray().flatten()
		children = np.where(th == 1)[0]
		children = children[children != excluded]
		leaves, no_leaves = [], []
		
		for x in children:
			leaves.append(x) if x in self.taxa else no_leaves.append(x)
		
		for child in no_leaves:
			leaves += self.leaves_from_node(child, node)
		
		return leaves


	def names_struc_from_node(self, node : int, excluded: int) -> list:
		"""
		Returns a two-dimensional list of descentdant leaves from a node. The 
		first dimension of the output indicates their ancestorship relative to 
		the inmediate descendants of the input node. It is required to select a 
		neighbor (`excluded`) of the latter node to orientate (root) the operation.
		"""


		th = self.adj_table[node].toarray().flatten()
		children = np.where(th == 1)[0]
		children = children[children != excluded]
		leaves = []
		no_leaves = []
		struc = [] 
		
		for x in children:
			leaves.append(x) if x in self.taxa else no_leaves.append(x)

		if len(leaves) > 0:
			struc.append(leaves)

		for child in no_leaves:
			struc.append(self.leaves_from_node(child, node))
		
		return struc
		# ===>>  Reorder nodes
		# ===>>  Create dict or Bloom filter (?) to store times a sp is inserted in the tree  
		#  [[24, 25], [27, 29, 30, 31]]
		# [[24], [25], [27, 29, 30, 31]]


	def orthology_test(self, target_node: int, excluded_node: int) -> bool:
		"""
		Test orthology condition on a node of an unrooted tree. The test excludes 
		one set of descendants of the node, to be set by `excluded_node` argument.
		Returns boolean. Results are stored in a sparse matrix (`self.result`) for 
		backreference. 
		"""
		
		pass_test = True
		#print(f"{target_node=}, {excluded_node=}")
		if self.results[target_node, excluded_node] == 0:
			
			#r = self.adj_table_[target_node]
			r = self.adj_table[target_node].toarray().flatten() # <<===
			icr = np.where(r == 1)[0]
			icr = icr[icr != excluded_node]
			icr = icr.tolist()
			#internal = [x for x in icr if not x in self.labels]
			#print(f"{icr=}")

			#if len(internal) > 0:
			if len(icr) > 0:
				falsies = []
				for i in icr:
					if self.orthology_test(i, target_node) == False:
						falsies.append(i)

				if len(falsies) > 0:
					pass_test = False

			if pass_test:				

				# Init name_struct to the final shape  <<===
				thleaves = self.names_struc_from_node(target_node, excluded_node)
				#print(f"{thleaves=}")

				# if target node is leaf (then thleaves == []), let's assume it is ortholog set
				if len(thleaves) > 0: 

					name_struc = []

					for brleaves in thleaves:
						name_struc.append(set([self.taxa[x] for x in brleaves]))


					#print(name_struc)

					superset = reduce(lambda x,y: x | y, name_struc)
					#print(f"{superset=}")
					if len(superset) == 1:
						pass_test = True

					else:
						for i,d in combinations(name_struc, 2):
							if len(i & d) > 0:
								pass_test = False
								break
			#print(f"{pass_test=}")
			if pass_test:
				self.results[target_node, excluded_node] = 2
			else:
				self.results[target_node, excluded_node] = 1


		elif self.results[target_node, excluded_node] == 1:
			pass_test = False

		elif self.results[target_node, excluded_node] == 2:
			pass_test = True

		return pass_test


	def ortholog_encoder(self, min_taxa : int, verbose : bool) -> list:
		"""
		Encodes orthologous groups found in the input tree into a matrix of shape
		leaves x orthologous group. Matrix is a two-dimensional list. Orthologous 
		sets are automatically encoded as ones.
		"""

		encoded_edges = []
		encoding = []
		labels = [x for x in self.labels]
		leaves = [x for x in self.taxa]		
		uniq_tax = set(self.taxa.values())
		get_table = True
		
		if len(uniq_tax) < min_taxa:
			
			warnings.warn(f"Tree in file {self.tree_file} contains less than {min_taxa} species (threshold set by the user).",
				stacklevel=2)

			if not verbose: get_table = False

		if get_table:
		
			if len(uniq_tax) == len(self.taxa):	# tree is perfect		

				warnings.warn(f"Tree in file {self.tree_file} is non-problematic (all terminals are different species).",
					stacklevel=2)
				encoding = [[1 for x in self.labels]]

			else:

				for i,d in zip(self.edge_coors_x, self.edge_coors_y):
					if not i in self.taxa:
						_ = self.orthology_test(i,d)
					if not d in self.taxa:
						_ = self.orthology_test(d,i)
				
				# Get starting nodes for traversal
				inits = {}
				for i,d in combinations(leaves, 2):
					pa0 = self.adj_table[i].toarray().flatten() # <<====
					pa1 = self.adj_table[d].toarray().flatten()
					pa0 = np.where(pa0 == 1)[0][0]
					pa1 = np.where(pa1 == 1)[0][0]

					if pa0 == pa1:
						#print(f"{inits=}")
						if not pa0 in inits:
							inits[pa0] = {i:0, d:0}

				# Find orthologous clades
				print(f"{inits=}")
				for start in inits:
					print(f"{start=}")
					prev_node = None
					curr_node = start
					#curr_excluded = list(inits[start].keys())
					curr_excluded = list(inits[curr_node].keys())
					still = True
					
					while still:
						print(f"{curr_node=}, {prev_node=}")

						neighs = self.get_neighbors(curr_node, curr_excluded)

						if neighs.shape[0] == 0: # probably in special case
							
							if curr_node in self.taxa:
								
								thnames = self.names_struc_from_node(prev_node, curr_node)
								thnames = reduce(lambda x,y: x+y, thnames)
								thnames = {self.taxa[x] for x in thnames}
								print(f"{thnames=}")

								if self.taxa[curr_node] in thnames: # maybe a single spp
								
									if len(thnames) == 1: # definitively a single spp

										warnings.warn(f"Tree in file {self.tree_file} contains a single species.", stacklevel=2)
										encoding = [[1 for x in self.labels]]

										return encoding


									else: # ===>> What could possibly be this case???
										pass

								else: # Simple set of non-problematic taxa, less species than leaves
								
									warnings.warn(f"Tree in file {self.tree_file} is non-problematic.", stacklevel=2)
									encoding = [[1 for x in self.labels]]
								
									return encoding

						#curr_excluded = [curr_node] + labels
						curr_excluded.append(curr_node)
						cands = []
						print(f"{neighs=}")

						for nei in neighs:
							test = self.results[curr_node, nei]
							#print(f"self.results[{curr_node}, {nei}] : {self.results[curr_node, nei]}")

							if test == 2:
								cands.append(nei)
							else:
								curr_excluded.append(nei)
						print(f"{cands=}")
						if len(cands) == 0:
							still = False
							#print("1")
					
						elif len(cands) == 1:
							prev_node = curr_node
							curr_excluded.append(curr_node)
							curr_node = cands[0]
							#print("2")
						else:
							choosen = None
							choosen_size = 0
					
							for c in cands:
								th = len(self.leaves_from_node(curr_node, c))
								#print(f"cand: {c}, leaves: {th}")
								if choosen_size < th:
									choosen = c
									choosen_size = th
					
							prev_node = curr_node
							curr_excluded += [c for c in cands if c != choosen]
							curr_node = choosen

					#print(f"Encoded edge: {curr_node} - {prev_node}")
					#print(f"Excluded: {curr_excluded}")

					# encode char
					if prev_node is None or (prev_node, curr_node) in encoded_edges:
						continue

					else:
						thchar = [0 for x in self.labels]
						#print(f"{curr_node=}, {prev_node=}")
						thleaves = self.leaves_from_node(prev_node, curr_node)
						thtaxa = {self.taxa[x] for x in thleaves}
						#print(f"{thtaxa=}")
						if len(thtaxa) >= min_taxa:
							#*******************************
							# Possible solutions to mini trees with low taxa:
							# - Store node richness in results np.ndarray
							# - Simply count the number of species at the beggining and 
							#   print a single zero column output
							#
							#*******************************
							for l in thleaves:
								thchar[labels.index(l)] = 1
							if len(set(thchar)) >= 2: # only append informative chars
								encoding.append(thchar) 
								encoded_edges.append((prev_node, curr_node))
			
		return encoding


	def tsv_table(self, min_spp : int = 3, verbose : bool = False):

		bffr = ''
		encoding = self.ortholog_encoder(min_taxa=min_spp, verbose=verbose)

		if len(encoding) > 0:
	
			charnames = [f"state_{i+1}" for i,x in enumerate(encoding)]
			header = 'sequence\t' + '\t'.join(charnames) + '\n'
		
			for idx, node in enumerate(self.labels):
				bffr += f'{self.labels[node]}'
				
				for ichar in range(len(encoding)):
					bffr += f'\t{encoding[ichar][idx]}'
				
				bffr += '\n'
			
			bffr = header + bffr
		
		return bffr


	def get_neighbors(self, node: int, excluded: list = []):
		"""Retrieve list of adjacent nodes."""
		th = self.adj_table[node].toarray().flatten()
		th = np.where(th > 0)[0]
		th = th[~np.isin(th, excluded)]
		return th


if __name__ == "__main__":
	import os

	tsv = False

	if tsv:
		for d, s, f in os.walk('test_trees/'):
			for filito in f:
				if filito.endswith('.newick'):
					print(filito)
					file = os.path.join(d, filito)
					root = file.rstrip('.newick')
					thnet = Tree(file)
					res = thnet.tsv_table(3)
					with open(f'{root}.tsv', 'w') as wh:
						wh.write(res)

	else:
		#tfile = "test_trees/ygob/3162.newick" # Perfect medium tree
		#tfile = "test_trees/ygob/3162_der.newick" # Perfect tree with a single duplicated species
		#tfile = "test_trees/ygob/3162_der_der.newick" # Medium tree in which clipping a single duplicated species makes a perfect case
		#tfile = "test_trees/ygob/322.newick" # medium size tree with 3 ortholog sets 
		#tfile = "test_trees/ygob/4741.newick" # Perfect small tree
		#tfile = "test_trees/ygob/301.newick"  # medium size tree with 2 ortholog sets
		#tfile = "test_trees/ygob/4777.newick"  #  small tree of a single sp
		tfile = "test_trees/ygob/4777_.newick"  #  small tree of two sp
		#tfile = "simple.newick"
		tr = Tree(tfile)
		#print(tr.list)
		#print("\n".join([f"{x[0]}:{x[1]}" for x in tr.labels.items()]))
		#print(tr.orthology_test(29, 5))
		#print(tr.orthology_test(2, 4))
		#print(tr.orthology_test(7, 6))
		#print(tr.orthology_test(6, 7))
		print("<<====  Table  ====>>")
		print(tr.tsv_table(1, verbose=True))