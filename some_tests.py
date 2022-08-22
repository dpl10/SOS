from time import time
import tracemalloc
import phylo
import os
import sys

ar = "test_trees/ygob"
count = 0
bffr = "File\tTerminals\tReading time (s)\tEncoding time (s)\tTotal time (s)\tCurrent memory (KB)\tPeak memory (KB)\n"

if os.path.exists(ar):

	for d, s, f in os.walk(ar):

		for infile in f:

			if infile.endswith(".newick"):

				print(infile)
				tracemalloc.start()
				time0 = time()
				tr = phylo.Tree(os.path.join(d, infile))
				time1 = time()
				itime = time1 - time0
				#print("\n\nInstantiation time:", itime, "seconds")
				#print(f"{len(tr.labels)=}")
				#l0 = tr.ortholog_encoder()
				l1 = tr.tsv_table(1, verbose=True, debug=False)
				print(l1)
				time2 = time()
				etime = time2 - time1
				#print("\n\nEncoding time:", etime, "seconds")
				ttime = time() - time0
				#print("\n\nTotal time:", ttime, "seconds")
				#print("\n\nMemory use:")
				mem = tracemalloc.get_traced_memory()
				#print(f"Current: {mem[0] / 1000} KB, peak: {mem[1] / 1000} KB")
				tracemalloc.stop()
				bffr += f"{infile}\t{len(tr.labels)}\t{itime}\t{etime}\t{ttime}\t{mem[0] / 1000}\t{mem[1] / 1000}\n"
				count += 1
				#if count > 10: break
				print("#" * 80)
				
else:
	print(f"{ar} is not a valid folder!")

with open("some_stats_ygob.tsv", "w") as fh:
	fh.write(bffr)

exit()
