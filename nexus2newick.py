import sys
import os

indir = None
for arg in sys.argv:
	if os.path.exists(arg):
		indir = arg

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