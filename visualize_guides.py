import matplotlib.pyplot as plt
from Bio import Entrez, SeqIO
from Bio.SeqUtils import gc_fraction

Entrez.email="skyng867@gmail.com"
handle=Entrez.efetch(db="nucleotide", id="NM_000518", rettype="fasta", retmode="text")
record=SeqIO.read(handle, "fasta")
handle.close()

seq=record.seq.upper()

positions=[]
gc_values=[]

for i in range(20, len(seq)-3):
    pam=seq[i:i+3]
    if pam[1]=="G" and pam[2]=="G":
        guide=seq[i-20:i]
        gc=gc_fraction(guide)*100
        if 40<=gc<=70:
            positions.append(i)
            gc_values.append(gc)

plt.figure(figsize=(12, 5))
plt.scatter(positions, gc_values, color="steelblue", alpha=0.6)
plt.axhline(y=60, color="red", linestyle="--", label="Ideal GC 60%")
plt.xlabel("Position in HBB gene")
plt.ylabel("GC content (%)")
plt.title("CRISPR guide RNA candidates -- HBB gene")
plt.legend()
plt.tight_layout()
plt.savefig("hbb_guides.png")
plt.show()

print(f"Plot saved as hbb_guides.png")