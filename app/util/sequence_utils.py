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
    
    # Validate sequence
    if not validate_sequence(str(sequence)):
        raise ValueError("Invalid DNA sequence")
    
    # Calculate GC content
    gc_count = sequence.count("G") + sequence.count("C")
    gc_content = (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0
    
    return {
        "gc_content": round(gc_content, 2),
        "length": len(sequence)
    }

def get_transcription(sequence):
    """
    Return the RNA transcription of a DNA sequence
    """
    if isinstance(sequence, str):
        sequence = Seq(sequence)
    
    try:
        # DNA to RNA transcription
        return str(sequence.transcribe())
    except Exception as e:
        print(f"Error transcribing sequence: {e}")
        return None

def get_translation(sequence):
    """
    Return the protein translation of a DNA sequence
    """
    if isinstance(sequence, str):
        sequence = Seq(sequence)
    
    try:
        # DNA to protein translation
        protein = str(sequence.translate())
        return protein if protein else None
    except Exception as e:
        print(f"Error translating sequence: {e}")
        return None

def validate_sequence(sequence):
    """
    Validate if a sequence is a valid DNA sequence
    """
    if isinstance(sequence, Seq):
        sequence = str(sequence)
    
    # Check if sequence only contains valid DNA characters (A, T, G, C, N)
    pattern = re.compile(r'^[ATGCN]+$', re.IGNORECASE)
    return bool(pattern.match(sequence))

def compare_sequences(sequences_dict):
    """
    Compare multiple sequences and return comparative analysis
    
    Args:
        sequences_dict (dict): Dictionary of sequence_id -> sequence
        
    Returns:
        dict: Comparative analysis results including GC content
    """
    if not sequences_dict or len(sequences_dict) < 2:
        raise ValueError("At least two sequences are required for comparison")
    
    # Analyze each sequence
    comparison_results = {
        "sequences": [],
        "gc_content": {
            "labels": [],
            "values": []
        }
    }
    
    # Process each sequence
    for seq_id, sequence in sequences_dict.items():
        if isinstance(sequence, str):
            sequence = Seq(sequence)
            
        # Calculate metrics
        analysis = analyze_sequence(sequence)
        
        # Add to comparison results
        comparison_results["sequences"].append({
            "seq_id": seq_id,
            "length": analysis["length"],
            "gc_content": analysis["gc_content"]
        })
        
        # Add data for visualization
        comparison_results["gc_content"]["labels"].append(seq_id)
        comparison_results["gc_content"]["values"].append(analysis["gc_content"])
    
    return comparison_results 