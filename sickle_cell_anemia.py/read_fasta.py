from Bio import SeqIO

normal = SeqIO.read("hbb_normal.fasta", "fasta")
sickle = SeqIO.read("hbb_sickle.fasta", "fasta")

print(f"ID:{normal.id}")
print(f"Description: {normal.description}")
print(f"Sequence: {normal.seq}")
print(f"Length: {len(normal.seq)}")

print(f"\nTranslated normal: {normal.seq.translate()}")
print(f"\nTranslated sickle: {sickle.seq.translate()}")