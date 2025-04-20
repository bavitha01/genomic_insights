from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import os
import re

# Set email for NCBI access - this is required for GenBank API access
Entrez.email = "bavithareddy08@gmail.com"  # Update with your real email for production use
Entrez.tool = "GenomicInsightsWebApp"
Entrez.api_key = os.environ.get("f8d4168273ea9ff4fe422025e6ec3fd65a08")  # Add your API key if you have one for higher rate limits

def process_fasta_file(file_path):
    """
    Process a FASTA file and return a dictionary of sequences
    """
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = record.seq
    
    if not sequences:
        raise ValueError("No valid sequences found in the FASTA file")
    
    return sequences

def fetch_sequence_from_genbank(sequence_id):
    """
    Fetch a sequence from GenBank by its ID
    """
    try:
        # First, check if the ID exists
        search_handle = Entrez.esearch(db="nucleotide", term=sequence_id, retmax=1)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        if int(search_results["Count"]) == 0:
            print(f"No results found for ID: {sequence_id}")
            return None
            
        # If ID exists, fetch the sequence
        handle = Entrez.efetch(db="nucleotide", id=sequence_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return record
    except Exception as e:
        print(f"Error fetching sequence {sequence_id}: {e}")
        return None

def analyze_sequence(sequence):
    """
    Analyze a DNA sequence and return various metrics
    """
    if isinstance(sequence, str):
        sequence = Seq(sequence)
    
    # Calculate GC content manually
    gc_count = sequence.count("G") + sequence.count("C")
    gc_content = (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0
    
    return {
        "gc_content": round(gc_content, 2),
        "length": len(sequence),
    }

def get_transcription(sequence):
    """
    Return the RNA transcription of a DNA sequence
    """
    if isinstance(sequence, str):
        sequence = Seq(sequence)
    
    # DNA to RNA transcription (replace T with U)
    return str(sequence.transcribe())

def get_translation(sequence):
    """
    Return the protein translation of a DNA sequence
    """
    if isinstance(sequence, str):
        sequence = Seq(sequence)
    
    # DNA to protein translation
    try:
        return str(sequence.translate())
    except Exception as e:
        print(f"Error translating sequence: {e}")
        return "Translation error"

def validate_sequence(sequence):
    """
    Validate if a sequence is a valid DNA sequence
    """
    if isinstance(sequence, Seq):
        sequence = str(sequence)
    
    # Check if sequence only contains valid DNA characters (A, T, G, C, N)
    pattern = re.compile(r'^[ATGCN]+$', re.IGNORECASE)
    return bool(pattern.match(sequence)) 
