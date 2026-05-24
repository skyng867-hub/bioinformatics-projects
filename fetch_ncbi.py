from Bio import Entrez, SeqIO

Entrez.email = "skyng867@gmail.com"

handle = Entrez.efetch(db="nucleotide", id="NM_000518", rettype="fasta", retmode="text")

record = SeqIO.read(handle, "fasta")
handle.close()

print(f"ID: {record.id}")
print(f"Length: {len(record.seq)}")
print(f"First 27 bases: {record.seq[:27]}")

start = record.seq.upper().find("ATGGTGCAT")
print(f"Coding sequence starts at posistion {start}")
print(f"First codon: {record.seq[start:start+3]}")

coding_seq = record.seq.upper()[start:]
protein = coding_seq.translate(to_stop=True)

print(f"Protien length: {len(protein)} amino acids")
print(f"Protein: {protein}")

normal_coding = record.seq.upper()[start:]
normal_str = str(normal_coding)
sickle_str = normal_str[:19] + "T" + normal_str[20:]

from Bio.Seq import Seq

normal_protein = Seq(normal_str).translate(to_stop=True)
sickle_protien = Seq(sickle_str).translate(to_stop=True)

print(f"\nNormal protein position 6-8: {normal_protein[5:8]}")
print(f"Sickle protein position 6-8: {sickle_protien[5:8]}")
print(f"Same length: {len(normal_protein) == len(sickle_protien)}")

print(f"Normal str posistion 13-18: {normal_str[13:18]}")
print(f"Sickle str posistion 13-18: {sickle_str[13:18]}")
print(f"Are they different: {normal_str != sickle_str}")

print(f"\n=== Sickle Cell Mutation Summary ===")
print(f"Gene: HBB (beta-globin)")
print(f"Position: codon 7")
print(f"DNA change: GAG → GTG")
print(f"Protein change: Glutamic acid (E) → Valine (V)")
print(f"Normal protein:  {normal_protein}")
print(f"Sickle protein:  {sickle_protien}")