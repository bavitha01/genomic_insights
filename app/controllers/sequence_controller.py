"""
Controller for sequence analysis operations.
This module handles the business logic for sequence analysis.
"""

from app.models.sequence_model import Sequence
from app.util.sequence_utils import (
    fetch_sequence_from_genbank,
    process_fasta_file,
    analyze_sequence,
    get_transcription,
    get_translation,
    validate_sequence
)

class SequenceController:
    """Controller for sequence analysis operations."""
    
    @staticmethod
    def analyze_from_fasta(file_path):
        """
        Analyze sequences from a FASTA file.
        
        Args:
            file_path (str): Path to the FASTA file
            
        Returns:
            list: List of analyzed sequence dictionaries
            
        Raises:
            ValueError: If no valid sequences found or file processing error
        """
        # Process the file
        sequences_dict = process_fasta_file(file_path)
        results = []
        
        for seq_id, sequence in sequences_dict.items():
            # Create a sequence object
            seq_obj = Sequence(seq_id, str(sequence))
            
            # Validate the sequence
            if not validate_sequence(sequence):
                results.append({
                    "seq_id": seq_id,
                    "error": "Invalid DNA sequence - contains non-nucleotide characters"
                })
                continue
            
            # Analyze the sequence
            analysis = analyze_sequence(sequence)
            seq_obj.gc_content = analysis["gc_content"]
            seq_obj.length = analysis["length"]
            seq_obj.transcription = get_transcription(sequence)
            seq_obj.translation = get_translation(sequence)
            
            # Add to results
            results.append(seq_obj.to_dict())
        
        return results
    
    @staticmethod
    def analyze_from_genbank(sequence_id):
        """
        Fetch and analyze a sequence from GenBank.
        
        Args:
            sequence_id (str): GenBank sequence ID
            
        Returns:
            dict: Analyzed sequence dictionary or None if not found
            
        Raises:
            Exception: If there's an error during fetching or analysis
        """
        record = fetch_sequence_from_genbank(sequence_id)
        
        if not record:
            return None
        
        # Create a sequence object from the record
        seq_obj = Sequence.from_genbank_record(record)
        
        # Analyze the sequence
        sequence = str(record.seq)
        analysis = analyze_sequence(sequence)
        seq_obj.gc_content = analysis["gc_content"]
        seq_obj.transcription = get_transcription(sequence)
        seq_obj.translation = get_translation(sequence)
        
        return seq_obj.to_dict() 