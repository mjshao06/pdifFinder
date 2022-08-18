from Bio import SeqIO


infile = "E:/Python/pdifFinder/Already/pS30-1/inputFile.fasta"
for rec in SeqIO.parse(infile, "fasta"):
    seq = str(rec.seq)[5404:9139]
    pdifmoduleseq = "E:/Python/pdifFinder/Already/pS30-1/pdifmoduleseq.txt"
    with open(pdifmoduleseq, 'a') as w:
        w.write(seq + '\n')
    print(rec.id)
    print(repr(rec.seq))
    print(len(rec))
