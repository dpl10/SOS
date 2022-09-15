from phylo import Tree
import re
import os
import warnings

tree_dir = "test_files"
t0 = "3162.newick" # Perfect medium tree
t1 = "4777.newick"  #  small tree of a single sp

def get_char_stats(tsv_text):

	chars = None
	last_char_states = None
	lines = re.split(r'\n', tsv_text)

	if re.search(r'sequence(\tstate_\d+)+', lines[0]):

		bits = re.split(r'\t', lines[0])
		chars = len(bits) - 1

		last_char_states = {x[-1] for x in lines[1:] if len(x) > 0}

	return (chars, last_char_states)
	


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



if __name__ == "__main__":

	test_tiny_single()