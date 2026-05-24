# ==============================
# 1. VARIABLES — no types, just assign
# ==============================
gene_name = "HBB"           # string
length = 27                 # integer
gc_content = 0.52           # float
is_mutant = True            # boolean

print(gene_name, length, gc_content, is_mutant)


# ==============================
# 2. STRINGS — DNA is just a string
# ==============================
seq = "ATGGTGCAC"

print(seq[0])        # A        — first character (index 0)
print(seq[-1])       # C        — last character
print(seq[0:3])      # ATG      — first codon (index 0 up to BUT NOT including 3)
print(seq[3:6])      # GTG      — second codon
print(len(seq))      # 9        — length

# Slice every 3 bases (codons)
print(seq[0:3])      # ATG
print(seq[3:6])      # GTG
print(seq[6:9])      # CAC


# ==============================
# 3. f-STRINGS — how you print cleanly
# ==============================
position = 19
normal = "A"
mutant = "T"

print(f"Position {position}: Normal={normal}  Mutant={mutant}")
# anything inside {} gets evaluated — that's it


# ==============================
# 4. FOR LOOP — the most important thing in bioinformatics
# ==============================
seq = "ATCG"

# Loop through every character
for base in seq:
    print(base)         # prints A, then T, then C, then G

# Loop with position number
for i, base in enumerate(seq):
    print(i, base)      # 0 A / 1 T / 2 C / 3 G

# enumerate() returns TWO things each step: (index, value)
# you catch them with two variable names: i, base


# ==============================
# 5. ZIP — walk two sequences together
# ==============================
normal = "GAGGAG"
sickle = "GTGGAG"

for n, s in zip(normal, sickle):
    print(n, s)         # G G / A T / G G / G G / A A / G G

# zip() pairs them up like a zipper: (G,G) (A,T) (G,G)...
# combine with enumerate to also get position:
for i, (n, s) in enumerate(zip(normal, sickle)):
    print(i, n, s)      # 0 G G / 1 A T / 2 G G ...
# note the extra () around (n, s) — required when unpacking inside enumerate


# ==============================
# 6. IF / ELIF / ELSE
# ==============================
base = "T"

if base == "A":
    print("Adenine")
elif base == "T":
    print("Thymine")    # this one runs
elif base == "G":
    print("Guanine")
else:
    print("Cytosine")

# != means "not equal"
# == means "is equal" (NOT assignment, that's just =)


# ==============================
# 7. LISTS — store many sequences, results, positions
# ==============================
mutations = []          # empty list

mutations.append(5)     # add to it
mutations.append(19)
mutations.append(47)

print(mutations)        # [5, 19, 47]
print(mutations[0])     # 5
print(len(mutations))   # 3

# Loop through a list
for pos in mutations:
    print(f"Mutation at position {pos}")


# ==============================
# 8. LIST COMPREHENSION — one-line loop (used constantly)
# ==============================
seq = "ATGGTGCAC"

# Old way — loop and append
gc_bases = []
for base in seq:
    if base == "G" or base == "C":
        gc_bases.append(base)

# New way — same thing in one line
gc_bases = [base for base in seq if base == "G" or base == "C"]

print(gc_bases)         # ['G', 'G', 'C', 'C']
print(len(gc_bases))    # 4

# pattern: [what_to_keep  for  variable  in  sequence  if  condition]


# ==============================
# 9. FUNCTIONS — wrap reusable logic
# ==============================
def gc_content(sequence):
    gc = [base for base in sequence if base == "G" or base == "C"]
    return len(gc) / len(sequence)

print(gc_content("ATGGTGCAC"))   # 0.555...
print(gc_content("AAATTTAAA"))   # 0.0

# def name(inputs):
#     do stuff
#     return result