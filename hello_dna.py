from Bio.Seq import Seq

dna = Seq("ATCGATCGTAGC")
print("My first DNA sequence:", dna)
print("Complement:", dna.complement())
print("Reverse complement:", dna.reverse_complement())