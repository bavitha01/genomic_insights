from sequence_analysis import fetch_sequence
def generate_report(sequence_ids):
    report = []
    for seq_id in sequence_ids:
        record = fetch_sequence(seq_id)
        if record:
            gc_content = (record.seq.count("G") + record.seq.count("C")) / len(record.seq) * 100
            report.append(f"Sequence ID: {seq_id}, GC Content: {gc_content:.2f}%")
        else:
            report.append(f"Failed to fetch sequence: {seq_id}")

    # Save the report to a text file in the results directory
    output_path = "results/analysis_report.txt"
    with open(output_path, "w") as f:
        f.write("\n".join(report))
    print(f"Report saved to: {output_path}")

# Example: Generate a report for BRCA1, TP53, and EGFR
generate_report(["NM_007294", "NM_000546", "NM_005228"])