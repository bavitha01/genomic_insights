from Bio import Entrez, SeqIO
from Bio.SeqUtils import gc_fraction  # Correct import

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
        gc_content = gc_fraction(record.seq) * 100
        return gc_content
    except Exception as e:
        print(f"Error calculating GC content: {e}")
        return None

# Example: Fetch BRCA1 gene sequence
brca1_record = fetch_sequence("NM_007294")
if brca1_record:
    print(f"Sequence ID: {brca1_record.id}")
    print(f"Sequence Description: {brca1_record.description}")
    print(f"Sequence Length: {len(brca1_record.seq)}")
    gc_content = calculate_gc_content(brca1_record)
    if gc_content is not None:
        print(f"GC Content: {gc_content}%")
    else:
        print("Failed to calculate GC content.")
else:
    print("Failed to fetch sequence.")