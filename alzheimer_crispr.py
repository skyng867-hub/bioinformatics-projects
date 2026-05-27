"""
Alzheimer's Disease - CRISPR Guide RNA Scanner
Genes: APP, PSEN1, Psen2
Pipline: Fetch from NCBI, find PAM sites, filter by GC%
"""
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

Entrez.email="skyng867@gmail.com"

genes = [
    ("APP",   "NM_000484"),
    ("PSEN1", "NM_000021"),
    ("PSEN2", "NM_000447"),
]

def find_guides(seq, min_gc=40, max_gc=70):
    guides=[]
    for i in range(20, len(seq)-3):
        pam=seq[i:i+3]
        if pam[1]=="G" and pam[2]=="G":
            guide=seq[i-20:i]
            gc=gc_fraction(guide)*100
            if min_gc<=gc<=max_gc:
                guides.append((i,str(guide), round(gc, 1)))
    return guides

print("Alzzheimer's disease - CRISPR guide RNA analysis")
print("="*40)

for gene_name, accession in genes:
    handle=Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    record=SeqIO.read(handle, "fasta")
    handle.close()

    seq=record.seq.upper()
    guides=find_guides(seq)

    guides.sort(key=lambda g: abs(g[2]-60))
    print(f"\nGene:{gene_name} ({accession})")
    print(f"Sequence length : {len(seq)} bases")
    print(f"Candidate guides: {len(guides)}")
    print(f"Top 3 guides:")
    for pos, guide, gc in guides[:3]:
        print(f"Posistion {pos:4d} | {guide} | GC={gc}%")
    
print("\n" + "="*40)

