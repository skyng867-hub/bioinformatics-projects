from Bio import Entrez, SeqIO
from Bio.Seq import Seq

Entrez.email = "skyng867@gmail.com"

handle = Entrez.efetch(db="nucleotide", id="NM_000518", rettype="fasta", retmode ="text")
record = SeqIO.read(handle, "fasta")
handle.close()

seq = record.seq.upper()
print(f"Fetched: {record.id}")
print(f"Length: {len(seq)}")

guides = []

for i in range(20, len(seq) - 3):
    pam = seq [i:i+3]
    if pam == "NGG" or (pam[1] == "G" and pam[2] == "G"):
        guide = seq[i-20:i]
        guides.append((i,guide))

print(f"\nFound {len(guides)} potential guide RNAs on forward strand")
print("\nFirst 5 candidates.")
for pos, guide in guides[:5]:
    print(f"Posistion {pos} : {guide}")

rev_seq = seq.reverse_complement()
rev_guides = []

for i in range(20, len(rev_seq) - 3):
    pam = rev_seq[i:i+3]
    if pam[1] == "G" and pam[2] == "G":
        guide = rev_seq[i-20:i]
        rev_guides.append((i,guide))

print (f"Found {len(rev_guides)} potential guides RNAs on reverse strand")
print(f"Total candidates: {len(guides) + len(rev_guides)}")

from Bio.SeqUtils import gc_fraction

print("\nFiltering by GC content(40-70%)")
good_guides = []

for pos, guide in guides:
    gc = gc_fraction(guide) * 100
    if 40 <= gc <= 70:
        good_guides.append((pos,guide,gc))
    
print(f"Guides passing GC filter: {len(good_guides)}")
print("\nTop 5 filtered guides:")

for pos, guide, gc in good_guides[:5]:
    print(f"Posistion {pos}: {guide} GC={gc:.1f}%")

from Bio.Blast import NCBIWWW, NCBIXML

best_guide = good_guides[0][1]
print(f"\nBLASTING guide: {best_guide}")
print("Searching against human genome...(around 30-60 sec)")

result_handle = NCBIWWW.qblast(
    "blastn",
    "nt",
    str(best_guide),
    entrez_query="Homo sapiens [organism]"
)

blast_records = NCBIXML.read(result_handle)

print(f"\nTop BLAST hits:")
for alignment in blast_records.alignments[:5]:
    for hsp in alignment.hsps[:1]:
        print(f"Match: {alignment.title[:60]}")
        print(f"Score: {hsp.score} Gaps: {hsp.gaps} Identites: {hsp.identities}/{hsp.align_length}")
        print()