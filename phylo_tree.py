# Import necessary libraries
from Bio import AlignIO, Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# Step 1: Load the input sequences from a .txt file
def load_sequences_from_txt(txt_file):
    sequences = []
    with open(txt_file, 'r') as file:
        for i, line in enumerate(file):
            line = line.strip()
            if line:
                # Create SeqRecord for each sequence
                seq_record = SeqRecord(Seq(line), id=f"sequence_{i+1}")
                sequences.append(seq_record)
    return sequences

txt_file = "input.txt"
sequences = load_sequences_from_txt(txt_file)

# Step 2: Create a Multiple Sequence Alignment object from the sequences
alignment = MultipleSeqAlignment(sequences)

# Step 3: Calculate the Distance Matrix
calculator = DistanceCalculator('identity')  # You can use 'identity' or other models
distance_matrix = calculator.get_distance(alignment)

# Step 4: Construct the Phylogenetic Tree using Neighbor-Joining method
constructor = DistanceTreeConstructor()
phylo_tree = constructor.nj(distance_matrix)  # Use nj() for Neighbor-Joining or upgma() for UPGMA

# Step 5: Visualize the Tree
Phylo.draw(phylo_tree)
Phylo.write(phylo_tree, "output_tree.xml", "phyloxml")  # Save the tree in phyloXML format
