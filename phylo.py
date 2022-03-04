import re
import numpy as np
from itertools import combinations

class Treell:

	def __init__(self, tnt_file):

		self.list = []
		self.lengths = {}
		self.node_count = 0
		self.labels = {}
		self.taxa = {}
		self.adj_table = None

		with open(tnt_file , "r") as fh:
			for line in fh:
				line = line.strip()
				if line.startswith("("):
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

						

	def get_parent(self, node):
		for no, des in self.list:
			if des == node:
				return no


	def unroot(self):
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

		return None


	def orthology_test(self, target_node, excluded_node):
		pass_test = True
		r = self.adj_table[target_node]
		icr = np.where(r == 1)[0]
		icr = icr[icr != excluded_node]
		print(icr)
		names = []
		name_origin = {}

		for child in icr:
			if child in self.taxa:
				names.append(self.taxa[child])
				name_origin[self.taxa[child]] = 1

			else:
				thnames = self.orthology_test(child, target_node)
				if len(thnames) == 0:
					pass_test = False
					name_origin = {}
					break
				else:
					names += thnames
					for tn in thnames:
						if tn in name_origin:
							name_origin[tn] += 1
						else:
							name_origin[tn] = 1

		print(names)
		print(name_origin)

		if len(name_origin) == 1:
			pass_test = True
		elif len(name_origin) > 1:
			for tn in name_origin:
				if name_origin[tn] > 1:
					pass_test = False
					break
		else:
			pass_test = False

		if pass_test:
			names = list(set(names))
			return names

		else:
			return []

		"""
		for descendants that are leaves:
			get names

		for descendants that are not leaves:
			call function on descendants, target_node as excluded_node
			get their name lists

		if there is at least two taxon names in descendants AND one of them is repeated in different descendant lists:

			fail, return an empty list

		else:
			pass, return list of unique names
		
		"""

	def ortholog_finder(self):
		pass
		"""
		for all internal edges:
			for each node in the edge:
				if node pass orthology property:
					append node to ortholog list

		"""

	# Ortholog identification iin the tree conducted using the maximum inclusion,
	# but not inclusion and multiples leaves for the same taxa only if they are
	# monophyletic. 
	# Evaluation different techniques may be done by measuring phylogenetic noise
	# in the final tree 



if __name__ == "__main__":

	tntfile = "toy.tree"

	al = Treell(tntfile)

	al.unroot()

	for pair in al.list:
		print(pair)
	for la in al.labels:
		print(la, al.labels[la])
	for a in al.taxa:
		print(a, al.taxa[a])
	for le in al.lengths:
		print(le, al.lengths[le])

	print(al.adj_table)

	print(al.orthology_test(1, 3))




				


	