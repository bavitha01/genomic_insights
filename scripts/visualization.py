import matplotlib.pyplot as plt
import os
from sequence_analysis import fetch_sequence, calculate_gc_content

# Function to plot GC content of multiple sequences
def plot_gc_content(sequence_ids):
    gc_values = []
    for seq_id in sequence_ids:
        record = fetch_sequence(seq_id)
        if record:
            gc_content = calculate_gc_content(record)
            gc_values.append(gc_content)
        else:
            print(f"Failed to fetch sequence: {seq_id}")

    # Plot the GC content
    plt.bar(sequence_ids, gc_values)
    plt.xlabel("Sequence ID")
    plt.ylabel("GC Content (%)")
    plt.title("GC Content of Sequences")

    # Save the plot to the results directory
    output_path = os.path.join("results", "gc_content_plot.png")
    plt.savefig(output_path)
    print(f"Plot saved to: {output_path}")
    plt.show()

# Example: Plot GC content for BRCA1, TP53, and EGFR
plot_gc_content(["NM_007294", "NM_000546", "NM_005228"])

