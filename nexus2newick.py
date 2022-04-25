import sys
import os

"""
Use:

	python nexus2newick.py <directory>

where:

`directory`: folder containing files in nexus format, extension `.tre`.
"""


indir = None
if os.path.exists(sys.argv[1]):
	indir = sys.argv[1]

for d, s, f in os.walk(indir):
	for file in f:
		if file.endswith('.tre'):
			print(file)
			outfile = file.rstrip('tre') + 'newick'
			with open(os.path.join(d, file), 'r') as fhandle:
				with open(os.path.join(d, outfile), 'w') as outhandle:
					for line in fhandle:
						line = line.strip()
						if line.startswith('(') and line.endswith(';'):
							outhandle.write(line)

exit()