from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

Entrez.email="skyng867@gmail.com"

genes = [
    ("APP",   "NM_000484"),
    ("PSEN1", "NM_000021"),
    ("PSEN2", "NM_000447"),
]

for gene_name, accession in genes:
    print(f"\n{'='*40}")
    print(f"Gene: {gene_name}   ({accession})")

    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    record=SeqIO.read(handle, "fasta")
    handle.close()
    seq=record.seq.upper()
    print(f"Length: {len(seq)} bases")

    guides=[]
    for i in range(20,len(seq)-3):
        pam=seq[i:i+3]
        if pam[1]=="G" and pam[2]=="G":
            guide=seq[i-20:i]
            gc=gc_fraction(guide)*100
            if 40<=gc<=70:
                guides.append((i,guide,gc))

    print(f"Guide RNAs passing GC filter: {len(guides)}")
    print(f"Best guide: {guides[0][1]}  GC={guides[0][2]:.1f}%")