"""
Genomic Insights Web Application
Main entry point for the Flask application.
"""

from flask import Flask, render_template, request, jsonify, send_from_directory, make_response
from flask_cors import CORS
import os
import logging
from werkzeug.utils import secure_filename
from app.routes import bp as main_bp

# Initialize Flask app
app = Flask(__name__, 
            template_folder='app/templates',
            static_folder='app/static')

# Enable CORS with specific settings
CORS(app)

# Register blueprints
app.register_blueprint(main_bp)

@app.route('/')
def root():
    """Redirect root to main index"""
    return render_template('index.html')

@app.route('/', defaults={'path': ''})
@app.route('/<path:path>')
def catch_all(path):
    """Handle all routes including OPTIONS requests"""
    if request.method == 'OPTIONS':
        response = make_response()
        response.headers.add('Access-Control-Allow-Origin', '*')
        response.headers.add('Access-Control-Allow-Headers', 'Content-Type,Authorization')
        response.headers.add('Access-Control-Allow-Methods', 'GET,PUT,POST,DELETE,OPTIONS')
        return response
    try:
        return render_template(f'{path}.html' if path else 'index.html')
    except:
        return render_template('index.html')

# Configure upload folder
UPLOAD_FOLDER = 'uploads'
if not os.path.exists(UPLOAD_FOLDER):
    os.makedirs(UPLOAD_FOLDER)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@app.route('/api/upload', methods=['POST'])
def upload_file():
    """Handle file uploads"""
    try:
        if 'file' not in request.files:
            return jsonify({'error': 'No file part'}), 400
        
        file = request.files['file']
        if file.filename == '':
            return jsonify({'error': 'No selected file'}), 400
        
        if file:
            filename = secure_filename(file.filename)
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            file.save(filepath)
            return jsonify({'message': 'File uploaded successfully', 'filename': filename})
    except Exception as e:
        logger.error(f'Error uploading file: {str(e)}')
        return jsonify({'error': 'An error occurred during file upload'}), 500

@app.route('/api/search', methods=['POST'])
def search_sequence():
    """Search for a sequence in uploaded files"""
    try:
        data = request.get_json()
        if not data or 'sequence' not in data:
            return jsonify({'error': 'No sequence provided'}), 400
        
        sequence = data['sequence'].upper()
        results = []
        
        # Search through uploaded files
        for filename in os.listdir(app.config['UPLOAD_FOLDER']):
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            with open(filepath, 'r') as f:
                content = f.read()
                if sequence in content:
                    results.append({
                        'filename': filename,
                        'matches': content.count(sequence)
                    })
        
        return jsonify({'results': results})
    except Exception as e:
        logger.error(f'Error searching sequence: {str(e)}')
        return jsonify({'error': 'An error occurred during sequence search'}), 500

@app.route('/health')
def health_check():
    """Health check endpoint"""
    return jsonify({'status': 'healthy'})

if __name__ == '__main__':
    app.run(debug=True, port=5004) 