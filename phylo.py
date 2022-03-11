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
		self.comp_record = None

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

		self.results = np.zeros_like(self.adj_table)

		return None

	def leaves_from_node(self, node : int, excluded: int) -> list:
		th = self.adj_table[node]
		children = np.where(th == 1)[0]
		children = children[children != excluded]
		leaves, no_leaves = [], []
		for x in children:
			leaves.append(x) if x in self.taxa else no_leaves.append(x)
		for child in no_leaves:
			leaves += self.leaves_from_node(child)
		return leaves

	def orthology_test_(self, target_node: int, excluded_node: int) -> bool:
		
		pass_test = True
		
		if self.results[target_node, excluded_node] == 0:

			r = self.adj_table[target_node]
			icr = np.where(r == 1)[0]
			icr = icr[icr != excluded_node]
			names = []
			name_origin = {}

			for child in icr:

				if child in self.taxa:
					names.append(self.taxa[child])

					if not self.taxa[child] in name_origin:
						name_origin[self.taxa[child]] = 1
					else:
						name_origin[self.taxa[child]] += 1

				else:
					thtest = self.orthology_test(child, target_node)

					if thtest:
						thnames = self.leaves_from_node(child, target_node)
						names += thnames

						for tn in thnames:
							if tn in name_origin:
								name_origin[tn] += 1

							else:
								name_origin[tn] = 1
					else:
						pass_test = False
						break

			#print(f"{names=}")
			#print(f"{name_origin=}")

			if len(name_origin) == 1:
				pass_test = True

			elif len(name_origin) > 1:

				for tn in name_origin:

					if name_origin[tn] > 1:
						pass_test = False
						break

			else:
				pass_test = False

		elif self.results[target_node, excluded_node] == 1:
			pass_test = False

		elif self.results[target_node, excluded_node] == 2:
			pass_test = True

		return pass_test

	def orthology_test(self, target_node: int, excluded_node: int) -> list:
		"""
		Test orthology condition on a node of an unrooted tree. The test excludes 
		one set of descendants of the node, to be set by `excluded_node` argument.
		Returns an empty list if fails, or a list of the terminal taxa if passes. 
		"""

		pass_test = True
		r = self.adj_table[target_node]
		icr = np.where(r == 1)[0]
		icr = icr[icr != excluded_node]
		names = []
		name_origin = {}

		for child in icr:

			if child in self.taxa:
				names.append(self.taxa[child])

				if not self.taxa[child] in name_origin:
					name_origin[self.taxa[child]] = 1
				else:
					name_origin[self.taxa[child]] += 1

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

		#print(f"{names=}")
		#print(f"{name_origin=}")

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


	def ortholog_finder(self):
		internal_ = np.zeros_like(self.adj_table)
		leaves = [x for x in range(self.adj_table.shape[0]) if x in self.taxa]
		internal_[leaves] = 0
		internal_[:,leaves] = 0
		int_coors = np.where(internal_ > 0)

		print(internal_)

		for i,d in zip(int_coors[0], int_coors[1]):
			testi = self.orthology_test(i,d)
			testd = self.orthology_test(i,d)

			if len(testi):
				#encode descendants excluding
				pass


if __name__ == "__main__":

	tntfile = "../toy.tree"

	al = Treell(tntfile)

	al.unroot()

	for pair in al.list:
		print(pair)
	#for la in al.labels:
	#	print(la, al.labels[la])
	for a in al.taxa:
		print(a, al.taxa[a])
	for le in al.lengths:
		print(le, al.lengths[le])

	print(al.adj_table)

	print(al.orthology_test(0, 15))

	al.ortholog_finder()


				


	