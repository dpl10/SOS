import re
import numpy as np
from itertools import combinations

class Treell:

	def __init__(self, tree_file):

		self.list = []
		self.lengths = {}
		self.node_count = 0
		self.labels = {}
		self.taxa = {}
		self.adj_table = None
		self.comp_record = None
		self.results = None

		#
		# Change numpy array types to np.int8
		# Move unroot func to constructor
		#


		with open(tree_file , "r") as fh:
			for line in fh:
				line = line.strip()

				if line.startswith("("):

					if re.search(',', line): # Newick
						line = re.sub(r'\)', ',)', line)
						line = re.sub(r'\),\(', ') (', line)
						line = re.sub(r',', ' ', line)
						line = re.sub(r'\):', ')=', line)
						line = re.sub(r':', '=', line)
						line = re.sub(r'\s+', ' ', line)
	
						#print(line)

					line = re.sub(r"\s+\)", ")", line)
					#print(line)
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
								#print(f"(: {papa_pointer} to {node_pointer}")
								
							elif char == ')':
								if len(label) > 0:
									self.labels[node_pointer] = label
									label = ""
								node_pointer = self.get_parent(node_pointer)
								#print(f"): back to {node_pointer}")
								
							elif char == " ":
								if len(label) > 0:
									self.labels[node_pointer] = label
									label = ""

								pa = self.get_parent(node_pointer)
								node_pointer = self.node_count
								self.node_count += 1
								self.list.append([pa, node_pointer])
								#print(f"Space: {pa} to {node_pointer}")

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

		root_edges_idx = []
		root_descendants = []

		for idx, edge in enumerate(self.list):
			if edge[0] == -1:
				root_descendants.append(edge[1])
				root_edges_idx.append(idx)
		
		for i,d in combinations(root_descendants, 2):
			self.list.append([i, d])

		self.list = [x for i,x in enumerate(self.list) if not i in root_edges_idx]
		self.adj_table = np.zeros((self.node_count, self.node_count))

		for i,d in self.list:
			self.adj_table[i,d] = 1
			self.adj_table[d,i] = 1

		self.results = np.zeros_like(self.adj_table)
		# Result table. Index meaning:
		# First dimension = main node
		# Second dimension = excluded node
		#
		# Values
		# 0 = Not computed yet
		# 1 = Negative
		# 2 = Positive


	def get_parent(self, node: int) -> int:
		""""To be used exclusively with the linked list."""
		for no, des in self.list:
			if des == node:
				return no


	def leaves_from_node(self, node : int, excluded: int) -> list:
		th = self.adj_table[node]
		children = np.where(th == 1)[0]
		children = children[children != excluded]
		leaves, no_leaves = [], []
		
		for x in children:
			leaves.append(x) if x in self.taxa else no_leaves.append(x)
		
		for child in no_leaves:
			leaves += self.leaves_from_node(child, node)
		
		return leaves


	def orthology_test(self, target_node: int, excluded_node: int) -> bool:
		"""
		Test orthology condition on a node of an unrooted tree. The test excludes 
		one set of descendants of the node, to be set by `excluded_node` argument.
		Returns boolean. Updates table at `self.result`. 
		"""
		
		pass_test = True
		
		if self.results[target_node, excluded_node] == 0:
			
			r = self.adj_table[target_node]
			icr = np.where(r == 1)[0]
			icr = icr[icr != excluded_node]
			name_origin = {}
			a_son_failed = False
			
			for child in icr:
				a_son_failed = False
				if child in self.taxa:
					if not self.taxa[child] in name_origin:
						name_origin[self.taxa[child]] = 1
					else:
						name_origin[self.taxa[child]] += 1

				else:
					thtest = self.orthology_test(child, target_node)

					if thtest:
						thleaves = self.leaves_from_node(child, target_node)
						thnames = [self.taxa[x] for x in thleaves]
						
						for tn in thnames:
							if tn in name_origin:
								name_origin[tn] += 1
							else:
								name_origin[tn] = 1

					else:
						a_son_failed = True
						break

			#print(f"{name_origin=}\n")

			if a_son_failed:
				pass_test = False

			else:
				if len(name_origin) == 1:
					pass_test = True

				elif len(name_origin) > 1:
					for tn in name_origin:
						if name_origin[tn] > 1:
							pass_test = False
							break

				else: # redundant?
					pass_test = False

			if pass_test:
				self.results[target_node, excluded_node] = 2
			else:
				self.results[target_node, excluded_node] = 1


		elif self.results[target_node, excluded_node] == 1:
			pass_test = False

		elif self.results[target_node, excluded_node] == 2:
			pass_test = True

		return pass_test


	def ortholog_encoder(self):
		encoded_edges = []
		encoding = []
		labels = [x for x in self.labels]
		internal = np.copy(self.adj_table)
		leaves = [x for x in self.taxa]
		internal[leaves] = 0
		internal[:,leaves] = 0
		int_coors = np.where(internal > 0)
		for i,d in zip(int_coors[0], int_coors[1]):
			_ = self.orthology_test(i,d)
			_ = self.orthology_test(d,i)
		
		# Get starting nodes for traversal
		inits = {}
		for i,d in combinations(leaves, 2):
			pa0 = self.adj_table[i]
			pa1 = self.adj_table[d]
			pa0 = np.where(pa0 == 1)[0][0]
			pa1 = np.where(pa1 == 1)[0][0]

			if pa0 == pa1:
				if not pa0 in inits:
					inits[pa0] = {i:0, d:0}
				else:
					inits[pa0][i] = 0
					inits[pa0][d] = 0

		# Find orthologous clades
		# print(inits)
		for start in inits:
			prev_node = None
			curr_node = start
			curr_excluded = list(inits[start].keys())
			still = True
			
			while still:
				neighs = self.get_neighbors(curr_node, curr_excluded)
				curr_excluded = [curr_node] + labels
				cands = []
			
				for nei in neighs:
					test = self.results[curr_node, nei]
			
					if test == 2:
						cands.append(nei)
					else:
						curr_excluded.append(nei)
			
				if len(cands) == 0:
					still = False
			
				elif len(cands) == 1:
					prev_node = curr_node
					curr_excluded.append(curr_node)
					curr_node = cands[0]
			
				else:
					choosen = None
					choosen_size = 0
			
					for c in cands:
						th = len(self.leaves_from_node(curr_node, c))
						if choosen_size < th:
							choosen = c
							choosen_size = th
			
					prev_node = curr_node
					curr_excluded += [c for c in cands if c != choosen]
					curr_node = choosen

			# encode char
			if (curr_node, prev_node) in encoded_edges:
				continue
			else:
				thchar = [0 for x in self.labels]
				thleaves = self.leaves_from_node(curr_node, prev_node)
				for l in thleaves:
					thchar[labels.index(l)] = 1
				if len(set(thchar)) >= 2: # only append informative chars
					encoding.append(thchar) 
					encoded_edges.append((curr_node, prev_node))
		
		return encoding


	def tsv_table(self):
		#####################################################
		# Should we add char names at the top of the table?
		#
		bffr = ''
		encoding = self.ortholog_encoder()
		for idx, node in enumerate(self.labels):
			bffr += f'{self.labels[node]}'
			for ichar in range(len(encoding)):
				bffr += f'\t{encoding[ichar][idx]}'
			bffr += '\n'
		return bffr


	def get_neighbors(self, node, excluded=[]):
		th = self.adj_table[node]
		th = np.where(th > 0)[0]
		th = th[~np.isin(th, excluded)]
		return th


if __name__ == "__main__":

	import os

	for d, s, f in os.walk('test_trees/'):
		for filito in f:
			if filito.endswith('.newick'):
				file = os.path.join(d, filito)
				root = file.rstrip('.newick')
				thnet = Treell(file)
				res = thnet.tsv_table()
				with open(f'{root}.tsv', 'w') as wh:
					wh.write(res)
