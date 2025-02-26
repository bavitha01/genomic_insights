import pandas as pd
from sequence_analysis import fetch_sequence
from sequence_analysis import calculate_gc_content
import os

# Function to save GC content results to a CSV file
def save_gc_content_results(sequence_ids):
    data = []
    for seq_id in sequence_ids:
        record = fetch_sequence(seq_id)
        if record:
            gc_content = calculate_gc_content(record)
            data.append({"Sequence ID": seq_id, "GC Content (%)": gc_content})
        else:
            print(f"Failed to fetch sequence: {seq_id}")

    # Convert the data to a DataFrame
    df = pd.DataFrame(data)

    # Save the results to a CSV file in the results directory
    output_path = os.path.join("results", "gc_content_results.csv")
    df.to_csv(output_path, index=False)
    print(f"Results saved to: {output_path}")

# Example: Save GC content results for BRCA1, TP53, and EGFR
save_gc_content_results(["NM_007294", "NM_000546", "NM_005228"])