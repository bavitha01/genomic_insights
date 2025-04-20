"""
Models for the Genomic Insights application.
This file contains data model classes for the application.
"""

class Sequence:
    """
    Represents a DNA sequence with analysis results.
    """
    def __init__(self, seq_id, sequence, description=None):
        """
        Initialize a sequence object.
        
        Args:
            seq_id (str): The identifier for the sequence
            sequence (str): The nucleotide sequence
            description (str, optional): Description of the sequence
        """
        self.seq_id = seq_id
        self.sequence = sequence
        self.description = description
        self.gc_content = None
        self.length = len(sequence) if sequence else 0
        self.transcription = None
        self.translation = None
        self.organism = None
        self.taxonomy = None
    
    def to_dict(self):
        """
        Convert the sequence object to a dictionary for JSON serialization.
        
        Returns:
            dict: Dictionary representation of the sequence
        """
        return {
            "seq_id": self.seq_id,
            "sequence": self.sequence,
            "description": self.description,
            "gc_content": self.gc_content,
            "length": self.length,
            "transcription": self.transcription,
            "translation": self.translation,
            "organism": self.organism,
            "taxonomy": self.taxonomy
        }
    
    @classmethod
    def from_genbank_record(cls, record):
        """
        Create a Sequence object from a GenBank record.
        
        Args:
            record: A Biopython GenBank record
            
        Returns:
            Sequence: A Sequence object populated with data from the record
        """
        seq_obj = cls(record.id, str(record.seq), record.description)
        
        # Extract additional information if available
        if hasattr(record, 'annotations'):
            if 'organism' in record.annotations:
                seq_obj.organism = record.annotations['organism']
            if 'taxonomy' in record.annotations:
                seq_obj.taxonomy = record.annotations['taxonomy']
        
        return seq_obj 