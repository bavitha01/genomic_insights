import matplotlib.pyplot as plt
from sequence_analysis import fetch_sequence
from sequence_analysis import calculate_gc_content
# Example: Plot GC content of multiple sequences
def plot_gc_content(sequence_ids):
    gc_values = []
    for seq_id in sequence_ids:
        record = fetch_sequence(seq_id)
        gc_values.append(calculate_gc_content(record))

    plt.bar(sequence_ids, gc_values)
    plt.xlabel("Sequence ID")
    plt.ylabel("GC Content (%)")
    plt.title("GC Content of Sequences")
    plt.show()

# Example: Plot GC content for BRCA1, TP53, and EGFR
plot_gc_content(["NM_007294", "NM_000546", "NM_005228"])