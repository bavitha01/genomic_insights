from Bio import Entrez, SeqIO

# Set your email for NCBI access
Entrez.email = "bavithareddy08@gmail.com"

# Fetch a sequence from NCBI
def fetch_sequence(sequence_id):
    try:
        handle = Entrez.efetch(db="nucleotide", id=sequence_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return record
    except Exception as e:
        print(f"Error fetching sequence: {e}")
        return None

# Calculate GC content
def calculate_gc_content(record):
    if record is None:
        print("Error: Record is None.")
        return None
    try:
        gc_content = (record.seq.count("G") + record.seq.count("C")) / len(record.seq) * 100
        return gc_content
    except Exception as e:
        print(f"Error calculating GC content: {e}")
        return None