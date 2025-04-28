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
    validate_sequence,
    compare_sequences
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
        
        # Get the sequence string
        sequence = str(record.seq)
        
        # Analyze the sequence
        analysis = analyze_sequence(sequence)
        
        # Update sequence object with analysis results
        seq_obj.sequence = sequence
        seq_obj.gc_content = analysis["gc_content"]
        seq_obj.length = analysis["length"]
        seq_obj.transcription = get_transcription(sequence)
        seq_obj.translation = get_translation(sequence)
        
        return seq_obj.to_dict()
    
    @staticmethod
    def compare_sequences(sequence_ids=None, file_path=None):
        """
        Compare multiple sequences by their IDs or from a FASTA file.
        
        Args:
            sequence_ids (list, optional): List of GenBank sequence IDs to compare
            file_path (str, optional): Path to FASTA file with sequences to compare
            
        Returns:
            dict: Comparison results
            
        Raises:
            ValueError: If no valid sequences found or not enough sequences for comparison
        """
        sequences_dict = {}
        
        # Get sequences from GenBank if IDs provided
        if sequence_ids and len(sequence_ids) >= 2:
            for seq_id in sequence_ids:
                record = fetch_sequence_from_genbank(seq_id)
                if record:
                    sequences_dict[record.id] = str(record.seq)
        
        # Get sequences from FASTA file if provided
        elif file_path:
            sequences_dict = process_fasta_file(file_path)
        
        else:
            raise ValueError("Either sequence IDs or a FASTA file path must be provided")
        
        # Make sure we have at least 2 sequences
        if len(sequences_dict) < 2:
            raise ValueError("At least two valid sequences are required for comparison")
        
        # Perform the comparison
        return compare_sequences(sequences_dict) 