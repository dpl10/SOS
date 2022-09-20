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
t6 = "301.newick"  # medium size tree with 4 ortholog sets
t7 = "4925.newick" # Star phylogeny of 4 terms - 2 taxa

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


def check_group_membership(term_group, term_matrix, char_matrix):
	char_number = None
	idxs = [term_matrix.index(x) for x in term_group]
	antiidxs = [x for x in range(len(term_matrix)) if not x in idxs]
	#print(f"{idxs=}")
	#print(f"{antiidxs=}")

	for col in range(char_matrix.shape[1]):
		#print(f"{col=}")

		pro = set(char_matrix[idxs, col])
		anti = set(char_matrix[antiidxs, col])
		#print(f"{pro=}")
		#print(f"{anti=}")

		if pro == {1} and anti == {0}:
			char_number = col
			break

	return char_number


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
		antiidx0 = [x for x in range(len(terms)) if x != idx0]
		idx1 = terms.index('Knag#KNAG0D02940b')
		antiidx1 = [x for x in range(len(terms)) if x != idx1]

		assert chars == 2
		assert len(last_states) == 2
		assert len(w) == 0
		assert mat.shape == (21,2)
		assert mat[idx0, 1] == 0
		assert mat[idx1, 0] == 0
		assert set(mat[antiidx0, 1]) == {1}
		assert set(mat[antiidx1, 0]) == {1}


def test_base_polytomy():

	tfile = os.path.join(tree_dir, t5)
	tr = Tree(tfile)

	with warnings.catch_warnings(record=True) as w:
		warnings.simplefilter("always")
		out = tr.tsv_table(1, verbose=True)
		#print(out)
		chars, last_states = get_char_stats(out)

		assert chars == 1
		assert len(last_states) == 1
		assert len(w) == 1
		assert re.search(r'Tree in file .+\.newick is a non-problematic \(but some species have multiple terminals\).', str(w[0].message))


def test_med_four():

	tfile = os.path.join(tree_dir, t6)
	tr = Tree(tfile)
	group_0 = ['Kafr#KAFR0B05530' , 'Tbla#TBLA0F03060', 'Tbla#TBLA0G01630']
	group_1 = ['Cgla#CAGL0L12254g', 'Kafr#KAFR0B02730', 'Knag#KNAG0G02040', 'Tdel#TDEL0F03900',
		'Zrou#ZYRO0C01628g', 'Ecym#Ecym_4307', 'Egos#AGR085W', 'Lthe#KLTH0G13662g', 
		'Lwal#Kwal_56.23791', 'Klac#KLLA0F19162g', 'Lklu#SAKL0H17006g', 'Ncas#NCAS0B05020', 
		'Ndai#NDAI0B01990', 'Suva#Suva_10.173', 'Skud#Skud_12.157', 'Scer#YLR089C', 
		'Smik#Smik_12.148', 'Tbla#TBLA0E04440', 'Tpha#TPHA0J00750', 'Vpol#Kpol_392.5']
	group_2 = ['Ncas#NCAS0B03840', 'Ndai#NDAI0J01360']
	group_3 = ['Suva#Suva_2.271', 'Skud#Skud_4.372', 'Scer#YDR111C', 'Smik#Smik_4.357', 
		'Tpha#TPHA0A01800', 'Vpol#Kpol_543.39']

	with warnings.catch_warnings(record=True) as w:
		warnings.simplefilter("always")
		out = tr.tsv_table(1, verbose=True)
		chars, last_states = get_char_stats(out)
		terms, mat = get_matrix_rep(out)
		
		char_nbs = []
		char_nbs.append(check_group_membership(group_0, terms, mat))
		char_nbs.append(check_group_membership(group_1, terms, mat))
		char_nbs.append(check_group_membership(group_2, terms, mat))
		char_nbs.append(check_group_membership(group_3, terms, mat))

		assert set(char_nbs) == {0, 1, 2, 3}
		assert chars == 4
		assert len(last_states) == 2
		assert len(w) == 0


def test_tiny_star():

	tfile = os.path.join(tree_dir, t7)
	tr = Tree(tfile)
	
	with warnings.catch_warnings(record=True) as w:
		warnings.simplefilter("always")
		out = tr.tsv_table(1, verbose=True)
		chars, last_states = get_char_stats(out)
		#print(chars, states)
		assert chars == 1
		assert last_states == set('0')
		#print(w[0].message)
		assert len(w) == 1
		assert re.search(r'Tree in file .+\.newick is uninformative \(star phylogeny\)', str(w[0].message))


if __name__ == "__main__":

	test_tiny_star()