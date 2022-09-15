from phylo import Tree
import re
import warnings

def test_perfect_mid():

	tfile = "test_trees/ygob/3162.newick" # Perfect medium tree
	tr = Tree(tfile)
	
	with warnings.catch_warnings(record=True) as w:
		warnings.simplefilter("always")
		out = tr.tsv_table(1, verbose=True)
		lines = re.split(r'\n', out)
		#print(lines)
		codes = {x[-1] for x in lines if len(x) > 0}

		assert len(re.split(r'\t', lines[0])) == 2
		assert codes == set('1')
		#print(w[0].message)
		assert len(w) == 1
		assert re.search(r"Tree in file .+\.newick is non\-problematic \(all terminals are different species\)\.", str(w[0].message))

if __name__ == "__main__":

	test_perfect_mid()