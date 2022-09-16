from phylo import Tree
import re
import os
import warnings
import numpy as np

tree_dir = "test_files"
t0 = "3162.newick" # Perfect medium tree
t1 = "4777.newick"  #  Small tree of a single sp
t2 = "322.newick" # Medium size tree with 4 ortholog sets
t3 = "3162_der.newick" # Perfect tree with a single duplicated species
t4 = "3162_der_der.newick" # Medium tree in which clipping a single duplicated species makes a perfect case
t5 = "149.newick" # Medium size tree with polytomies, one at base

def get_char_stats(tsv_text):

	chars = None
	last_char_states = None
	lines = re.split(r'\n', tsv_text)

	if re.search(r'sequence(\tstate_\d+)+', lines[0]):

		bits = re.split(r'\t', lines[0])
		chars = len(bits) - 1

		last_char_states = {x[-1] for x in lines[1:] if len(x) > 0}

	return (chars, last_char_states)


def get_matrix_rep(tsv_text):

	lines = re.split(r'\n', tsv_text)
	terms = []
	term_dat = []

	if re.search(r'sequence(\tstate_\d+)+', lines[0]):
		for line in lines[1:]:
			if len(line) > 0:

				bits = re.split(r'\t', line)
				terms.append(bits[0])
				term_dat.append(bits[1:])

		term_dat = np.array(term_dat, np.int8)

	return (terms, term_dat)


def test_perfect_mid():

	tfile = os.path.join(tree_dir, t0)
	tr = Tree(tfile)
	
	with warnings.catch_warnings(record=True) as w:
		warnings.simplefilter("always")
		out = tr.tsv_table(1, verbose=True)
		chars, last_states = get_char_stats(out)

		assert chars == 1
		assert last_states == set('1')
		#print(w[0].message)
		assert len(w) == 1
		assert re.search(r"Tree in file .+\.newick is non\-problematic \(all terminals are different species\)\.", str(w[0].message))


def test_tiny_single():

	tfile = os.path.join(tree_dir, t1)
	tr = Tree(tfile)
	
	with warnings.catch_warnings(record=True) as w:
		warnings.simplefilter("always")
		out = tr.tsv_table(1, verbose=True)
		chars, last_states = get_char_stats(out)
		#print(chars, states)
		assert chars == 1
		assert last_states == set('1')
		#print(w[0].message)
		assert len(w) == 1
		assert re.search(r'Tree in file .+\.newick contains a single species.', str(w[0].message))


def test_med_four():

	tfile = os.path.join(tree_dir, t2)
	tr = Tree(tfile)

	with warnings.catch_warnings(record=True) as w:
		warnings.simplefilter("always")
		out = tr.tsv_table(1, verbose=True)
		chars, last_states = get_char_stats(out)

		assert chars == 4
		assert len(last_states) == 2
		assert len(w) == 0


def test_perfect_duplication():

	tfile = os.path.join(tree_dir, t3)
	tr = Tree(tfile)

	with warnings.catch_warnings(record=True) as w:
		warnings.simplefilter("always")
		out = tr.tsv_table(1, verbose=True)
		#print(w[0].message)
		chars, last_states = get_char_stats(out)

		assert chars == 1
		assert len(last_states) == 1
		assert len(w) == 1
		assert re.search(r'Tree in file .+\.newick is a non-problematic \(but some species have multiple terminals\).', str(w[0].message))
		

def test_almost_perfect():

	tfile = os.path.join(tree_dir, t4)
	tr = Tree(tfile)

	with warnings.catch_warnings(record=True) as w:
		warnings.simplefilter("always")
		out = tr.tsv_table(1, verbose=True)
		chars, last_states = get_char_stats(out)
		terms, mat = get_matrix_rep(out)
		idx0 = terms.index('Knag#KNAG0D02940')
		idx1 = terms.index('Knag#KNAG0D02940b')

		assert chars == 2
		assert len(last_states) == 2
		assert len(w) == 0
		assert mat.shape == (21,2)
		assert mat[idx0, 1] == 0
		assert mat[idx1, 0] == 0
		assert np.unique(mat[~idx0, 1]) == [1]
		assert np.unique(mat[~idx1, 0]) == [1]


def test_base_polytomy():

	tfile = os.path.join(tree_dir, t5)
	tr = Tree(tfile)

	with warnings.catch_warnings(record=True) as w:
		warnings.simplefilter("always")
		out = tr.tsv_table(1, verbose=True)
		print(out)
		chars, last_states = get_char_stats(out)

		assert chars == 1
		assert len(last_states) == 1
		assert len(w) == 1
		assert re.search(r'Tree in file .+\.newick is a non-problematic \(but some species have multiple terminals\).', str(w[0].message))



if __name__ == "__main__":

	test_base_polytomy()