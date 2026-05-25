from Bio.Seq import Seq

normal_dna = Seq("ATGGTGCACCTGACTCCTGAGGAGAAG")
sickle_dna = Seq("ATGGTGCACCTGACTCCTGTGGAGAAG")

print("DNA sequences")
print(f"Normal: {normal_dna}")
print(f"Sickle: {sickle_dna}")

print("\nWhere is the mutation?")

for i, (normal_base, sickle_base) in enumerate(zip(normal_dna, sickle_dna)):
    if normal_base != sickle_base:
        print(f"Position {i + 1}: Normal = {normal_base} Sickle = {sickle_base}")

normal_protein = normal_dna.translate()
sickle_protein = sickle_dna.translate()

print("\n Amino acid sequence")
print(f"Normal: {normal_protein}")
print(f"Sickle: {sickle_protein}")

for i, (n_aa, s_aa) in enumerate(zip(normal_protein, sickle_protein)):
    if n_aa != s_aa:
        print(f"Amino acid #{i + 1}: Normal = {n_aa} Sickle = {s_aa}")
        print("E = Glutamic acid (charged, water loving)")
        print("V = Valine        (greasy,water-loving)")
        print("V = Valine        (greasy, water-hating → clumps)")                                          



def find_mutations(seq1, seq2, label1 = "Normal", label2 = "Mutant"):
    mutations = []

    for i, (a, b) in enumerate(zip(seq1, seq2)):
        if a != b:
            mutations.append(i + 1)
            print(f"Position {i + 1}: {label1} = {a} {label2}={b}")
            print(f"\nTotal differences found: {len(mutations)}")
            return mutations


print("DNA differences")
dna_diffs = find_mutations(normal_dna, sickle_dna)

print("\nProtein differences")
aa_diffs = find_mutations(normal_protein, sickle_protein, label1="Glu", label2="Val")
