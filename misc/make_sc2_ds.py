from Bio import SeqIO
import pandas as pd

allele_matrix = pd.read_csv("/home/amorris/BioInf/Ardal/data/usher_barcodes.csv", index_col=0)
seqrecs = SeqIO.parse("/home/amorris/BioInf/Ardal/data/NC_045512_Hu-1.fasta", "fasta")

refseq = str(next(seqrecs).seq)  # More efficient way to get the first sequence

with open("./data/usher_seqs.fasta", "w") as fout:
    for sample_id, row in allele_matrix.iterrows():
        fout.write(f">{sample_id}\n")
        mod_seq = list(refseq)

        alleles = row[row == 1].index.tolist()

        for a in alleles:
            ref = a[0]
            p = int(a[1:-1]) - 1
            alt = a[-1]

            if refseq[p] == ref:
                mod_seq[p] = alt
            else:
                print(f"Error: Reference base mismatch at position {p+1} for allele {a}. Expected {ref}, found {refseq[p]}. Skipping this allele.")

        fout.write("".join(mod_seq) + "\n")
