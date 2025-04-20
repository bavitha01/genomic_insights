# Genomic Insights

A bioinformatics web application designed to analyze and visualize genomic data. This project integrates with the GenBank API to allow users to search for sequences by ID and provides tools for analyzing DNA sequences, including GC content calculation, transcription, and translation.

![Genomic Insights](https://img.shields.io/badge/Genomic-Insights-blue)

## Features

- **FASTA File Upload**: Upload and analyze multiple DNA sequences in FASTA format
- **GenBank API Integration**: Search for sequences by their ID from NCBI's GenBank database
- **Sequence Analysis**: Calculate GC content, sequence length, and other metrics
- **Transcription & Translation**: Convert DNA to RNA and protein sequences
- **Interactive Visualization**: User-friendly interface with color-coded sequences
- **Downloadable Results**: Export analysis results in FASTA format

## Installation

1. **Clone the repository**:
   ```bash
   git clone https://github.com/bavitha01/genomic_insights
   cd genomic_insights
   ```

2. **Create and activate a virtual environment**:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Mac/Linux
   venv\Scripts\activate     # On Windows
   ```

3. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

1. **Start the application**:
   ```bash
   python app.py
   ```

2. **Access the web interface**:
   Open your browser and navigate to `http://127.0.0.1:5000` (or `http://127.0.0.1:5001` if port 5000 is in use)

3. **Use the application**:
   - Upload a FASTA file using the upload form
   - Search for a sequence by GenBank ID (e.g., NM_007294 for BRCA1)
   - View analysis results for GC content, transcription, and translation

## API Access

The application integrates with the NCBI GenBank API. To improve API access reliability:

1. Update the email in `app/util/sequence_utils.py`:
   ```python
   Entrez.email = "bavithareddy08@gmail.com"  # Replace with your actual email
   ```

2. For frequent use, consider [obtaining an NCBI API key](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/) and adding it to the code:
   ```python
   import os
   Entrez.api_key = os.environ.get("NCBI_API_KEY")
   ```

## Project Structure

The project follows an MVC (Model-View-Controller) architecture:

```
genomic_insights/
├── app/                       # Application package
│   ├── controllers/           # Controllers (business logic)
│   │   └── sequence_controller.py
│   ├── models/                # Data models
│   │   └── sequence_model.py
│   ├── static/                # Static assets
│   │   ├── css/               # CSS styles
│   │   ├── img/               # Images
│   │   └── js/                # JavaScript files
│   ├── templates/             # HTML templates (views)
│   │   └── index.html
│   └── util/                  # Utility functions
│       └── sequence_utils.py  # Sequence analysis utilities
├── uploads/                   # Directory for uploaded files
├── sample_data/               # Sample data for testing
├── app.py                     # Main application entry point
├── requirements.txt           # Project dependencies
└── README.md                  # Project documentation
```

## Example Data

To test the application, you can use:

1. **Sample FASTA file** in the `sample_data` directory
2. **GenBank IDs** like:
   - NM_007294 - BRCA1 (Breast Cancer 1)
   - NM_000546 - TP53 (Tumor Protein p53)
   - NM_005228 - EGFR (Epidermal Growth Factor Receptor)

## Technologies Used

- **Backend**: Python, Flask, Biopython
- **Frontend**: HTML, CSS, JavaScript, Bootstrap 5
- **Architecture**: MVC (Model-View-Controller)
- **API**: NCBI GenBank API (Entrez E-utilities)

## Future Enhancements

- Multiple sequence alignment
- Phylogenetic tree generation
- Secondary structure prediction
- Interactive sequence visualization tools
- Batch processing of multiple GenBank IDs

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgements

- [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/) for providing access to genetic sequence data
- [Biopython](https://biopython.org/) for providing bioinformatics tools in Python
- [Flask](https://flask.palletsprojects.com/) for the web framework

