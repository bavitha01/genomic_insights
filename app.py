"""
Genomic Insights Web Application
Main entry point for the Flask application.
"""

from flask import Flask, render_template, request, jsonify
import os
import traceback
from app.controllers.sequence_controller import SequenceController

# Initialize Flask app
app = Flask(__name__, template_folder='app/templates', static_folder='app/static')

# Configure upload folder
UPLOAD_FOLDER = 'uploads'
if not os.path.exists(UPLOAD_FOLDER):
    os.makedirs(UPLOAD_FOLDER)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

@app.route('/')
def index():
    """Render the main page"""
    return render_template('index.html')

@app.route('/upload', methods=['POST'])
def upload_file():
    """Handle FASTA file upload"""
    if 'fastaFile' not in request.files:
        return jsonify({"error": "No file part"}), 400
    
    file = request.files['fastaFile']
    if file.filename == '':
        return jsonify({"error": "No selected file"}), 400
    
    if file:
        # Save the file
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], file.filename)
        file.save(filepath)
        
        # Process the file using the controller
        try:
            results = SequenceController.analyze_from_fasta(filepath)
            return jsonify({"results": results})
        except Exception as e:
            app.logger.error(f"Error processing file: {str(e)}")
            app.logger.error(traceback.format_exc())
            return jsonify({"error": str(e)}), 500

@app.route('/search', methods=['POST'])
def search_sequence():
    """Search for a sequence by ID in GenBank"""
    data = request.json
    seq_id = data.get('seq_id')
    
    if not seq_id:
        return jsonify({"error": "No sequence ID provided"}), 400
    
    try:
        # Use controller to fetch and analyze the sequence
        result = SequenceController.analyze_from_genbank(seq_id)
        
        if result:
            return jsonify({"result": result})
        else:
            return jsonify({"error": f"Sequence not found for ID: {seq_id}"}), 404
    except Exception as e:
        app.logger.error(f"Error searching sequence: {str(e)}")
        app.logger.error(traceback.format_exc())
        return jsonify({"error": str(e)}), 500

@app.route('/healthcheck')
def healthcheck():
    """Simple health check endpoint for testing the application"""
    return jsonify({"status": "ok"})

if __name__ == '__main__':
    # Try to use a different port if 5000 is in use (common on macOS)
    try:
        app.run(debug=True, port=5000)
    except OSError:
        print("Port 5000 is already in use, trying port 5001...")
        app.run(debug=True, port=5001) 