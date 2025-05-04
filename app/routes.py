from flask import Blueprint, request, jsonify, render_template, current_app
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import os
from werkzeug.utils import secure_filename
import io
from app import app
bp = Blueprint('main', __name__)

# Configure Entrez
Entrez.email = "your_email@example.com"  # Replace with your email

@bp.route('/options')
def options():
    """Render the options page"""
    return render_template('options.html')

@bp.route('/analyze')
def analyze_page():
    """Render the analyze page"""
    return render_template('analyze.html')

@bp.route('/api/analyze', methods=['POST'])
def analyze_sequence():
    """Analyze a DNA sequence"""
    try:
        if request.files and 'file' in request.files:
            file = request.files['file']
            if file.filename == '':
                return jsonify({"error": "No file selected"}), 400
            
            if file:
                # Read the FASTA file
                fasta_content = file.read().decode('utf-8')
                sequences = list(SeqIO.parse(io.StringIO(fasta_content), 'fasta'))
                if not sequences:
                    return jsonify({"error": "No valid sequences found in file"}), 400
                
                sequence = str(sequences[0].seq)
        else:
            data = request.get_json()
            if not data:
                return jsonify({"error": "No data provided"}), 400
                
            sequence_id = data.get('sequence_id')
            if sequence_id:
                # Fetch sequence from NCBI
                try:
                    handle = Entrez.efetch(db="nucleotide", id=sequence_id, rettype="fasta", retmode="text")
                    sequence_record = SeqIO.read(handle, "fasta")
                    handle.close()
                    sequence = str(sequence_record.seq)
                except Exception as e:
                    return jsonify({"error": f"Failed to fetch sequence: {str(e)}"}), 400
            else:
                sequence = data.get('sequence', '').upper()

        if not sequence:
            return jsonify({"error": "No sequence provided"}), 400

        # Create a Seq object
        seq = Seq(sequence)
        
        # Perform analysis
        rna = seq.transcribe()
        protein = seq.translate()
        gc_content = gc_fraction(seq) * 100  # Convert to percentage
        
        return jsonify({
            "dna_sequence": str(seq),
            "rna_sequence": str(rna),
            "protein_sequence": str(protein),
            "length": len(seq),
            "gc_content": gc_content
        })
    except Exception as e:
        return jsonify({"error": str(e)}), 400

@bp.route('/api/compare', methods=['POST'])
def compare_sequences():
    sequences_data = []
    
    try:
        if request.files and 'file' in request.files:
            # Handle FASTA file upload
            file = request.files['file']
            if not file.filename:
                return jsonify({'error': 'No file selected'}), 400
                
            # Read sequences from FASTA file
            fasta_content = file.read().decode('utf-8')
            sequences = parse_fasta(fasta_content)
            
            if len(sequences) < 2:
                return jsonify({'error': 'FASTA file must contain at least 2 sequences'}), 400
                
            for seq_id, sequence in sequences.items():
                gc_content = calculate_gc_content(sequence)
                sequences_data.append({
                    'seq_id': seq_id,
                    'length': len(sequence),
                    'gc_content': gc_content
                })
                
        elif request.is_json:
            # Handle sequence IDs
            data = request.get_json()
            sequence_ids = data.get('sequence_ids', [])
            
            if len(sequence_ids) < 2:
                return jsonify({'error': 'Please provide at least 2 sequence IDs'}), 400
                
            for seq_id in sequence_ids:
                sequence = fetch_sequence(seq_id)  # Implement this function to fetch from NCBI
                if sequence:
                    gc_content = calculate_gc_content(sequence)
                    sequences_data.append({
                        'seq_id': seq_id,
                        'length': len(sequence),
                        'gc_content': gc_content
                    })
                else:
                    return jsonify({'error': f'Could not fetch sequence for ID: {seq_id}'}), 400
        else:
            return jsonify({'error': 'Invalid request format'}), 400
            
        if not sequences_data:
            return jsonify({'error': 'No valid sequences found'}), 400
            
        # Prepare data for GC content chart
        gc_content_data = {
            'labels': [seq['seq_id'] for seq in sequences_data],
            'values': [seq['gc_content'] for seq in sequences_data]
        }
        
        return jsonify({
            'sequences': sequences_data,
            'gc_content': gc_content_data
        })
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500

def parse_fasta(fasta_content):
    """Parse FASTA content into a dictionary of sequences."""
    sequences = {}
    current_id = None
    current_seq = []
    
    for line in fasta_content.split('\n'):
        line = line.strip()
        if not line:
            continue
            
        if line.startswith('>'):
            if current_id:
                sequences[current_id] = ''.join(current_seq)
            current_id = line[1:].split()[0]  # Get the first word after '>'
            current_seq = []
        else:
            current_seq.append(line)
            
    if current_id:
        sequences[current_id] = ''.join(current_seq)
        
    return sequences

def calculate_gc_content(sequence):
    """Calculate GC content percentage of a sequence."""
    if not sequence:
        return 0
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    total = len(sequence)
    return round((gc_count / total) * 100, 2) if total > 0 else 0

@bp.route('/main/analyze')
def main_analyze_page():
    """Render the analyze page - main route"""
    return render_template('analyze.html')

@bp.route('/')
def index():
    """Render the home page"""
    return render_template('index.html')

@bp.route('/compare')
def compare():
    """Render the compare sequences page"""
    return render_template('compare.html')

def fetch_sequence(seq_id):
    """Fetch a sequence from NCBI using its ID."""
    try:
        # Set your email for NCBI's records
        Entrez.email = "your_email@example.com"  # Replace with your email
        
        # Fetch the sequence
        handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
        record = SeqIO.read(io.StringIO(handle.read()), "fasta")
        handle.close()
        
        return str(record.seq)
    except Exception as e:
        print(f"Error fetching sequence {seq_id}: {str(e)}")
        return None 