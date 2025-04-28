document.addEventListener('DOMContentLoaded', function() {
    // DOM Elements
    const sequenceInput = document.getElementById('sequenceInput');
    const fileInput = document.getElementById('fileInput');
    const dropZone = document.getElementById('dropZone');
    const analyzeButton = document.getElementById('analyzeSequence');
    const resultsContainer = document.getElementById('resultsContainer');
    const loadingOverlay = document.querySelector('.loading-overlay');
    const dnaResult = document.getElementById('dnaResult');
    const rnaResult = document.getElementById('rnaResult');
    const proteinResult = document.getElementById('proteinResult');
    const downloadButton = document.getElementById('downloadResults');
    const sequenceLength = document.getElementById('sequenceLength');
    const gcContent = document.getElementById('gcContent');
    const fileTab = document.getElementById('file');

    // Add analyze button for file upload section
    const fileAnalyzeButton = document.createElement('button');
    fileAnalyzeButton.type = 'button';
    fileAnalyzeButton.className = 'btn option-btn analyze-btn mt-3';
    fileAnalyzeButton.innerHTML = '<i class="bi bi-play-circle me-2"></i>Analyze File';
    fileAnalyzeButton.style.display = 'none'; // Hide initially
    dropZone.insertAdjacentElement('afterend', fileAnalyzeButton);

    // Show/hide file analyze button when file is selected
    fileInput.addEventListener('change', function(e) {
        const file = e.target.files[0];
        if (file) {
            fileAnalyzeButton.style.display = 'block';
        } else {
            fileAnalyzeButton.style.display = 'none';
        }
    });

    // File analyze button click handler
    fileAnalyzeButton.addEventListener('click', function() {
        const file = fileInput.files[0];
        if (file) {
            handleFile(file);
        }
    });

    // Store current results
    let currentResults = null;

    // Toast notification system
    const toastContainer = document.createElement('div');
    toastContainer.className = 'toast-container position-fixed bottom-0 end-0 p-3';
    document.body.appendChild(toastContainer);

    function showToast(message, type = 'error') {
        const toast = document.createElement('div');
        toast.className = `toast align-items-center border-0 ${type === 'error' ? 'bg-danger' : 'bg-success'} text-white`;
        toast.setAttribute('role', 'alert');
        toast.innerHTML = `
            <div class="d-flex">
                <div class="toast-body">
                    <i class="bi ${type === 'error' ? 'bi-exclamation-circle' : 'bi-check-circle'} me-2"></i>
                    ${message}
                </div>
                <button type="button" class="btn-close btn-close-white me-2 m-auto" data-bs-dismiss="toast"></button>
            </div>
        `;
        toastContainer.appendChild(toast);
        const bsToast = new bootstrap.Toast(toast);
        bsToast.show();
        toast.addEventListener('hidden.bs.toast', () => toast.remove());
    }

    // Download results handler
    downloadButton.addEventListener('click', () => {
        if (!currentResults) return;

        const resultsText = formatResultsForDownload(currentResults);
        const blob = new Blob([resultsText], { type: 'text/plain' });
        const url = window.URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = 'sequence_analysis_results.txt';
        document.body.appendChild(a);
        a.click();
        window.URL.revokeObjectURL(url);
        document.body.removeChild(a);
    });

    function formatResultsForDownload(results) {
        return `Sequence Analysis Results
======================

DNA Sequence:
${results.dna_sequence || 'No DNA sequence available'}

RNA Sequence:
${results.rna_sequence || 'No RNA sequence available'}

Protein Sequence:
${results.protein_sequence || 'No protein sequence available'}

Statistics:
- Sequence Length: ${results.length || 0} bp
- GC Content: ${results.gc_content || 0}%

Color Legend:
- A (Adenine): Red
- T (Thymine): Teal
- G (Guanine): Blue
- C (Cytosine): Green
- U (Uracil): Yellow
`;
    }

    // Input validation
    function validateSequenceId(sequenceId) {
        if (!sequenceId) {
            showToast('Please enter a sequence ID');
            sequenceInput.classList.add('is-invalid');
            return false;
        }
        sequenceInput.classList.remove('is-invalid');
        return true;
    }

    // File Upload Handling with visual feedback
    dropZone.addEventListener('dragover', (e) => {
        e.preventDefault();
        dropZone.style.borderColor = '#2a5298';
        dropZone.style.backgroundColor = 'rgba(42, 82, 152, 0.05)';
        dropZone.classList.add('upload-zone-active');
    });

    dropZone.addEventListener('dragleave', (e) => {
        e.preventDefault();
        dropZone.style.borderColor = 'rgba(42, 82, 152, 0.3)';
        dropZone.style.backgroundColor = '';
        dropZone.classList.remove('upload-zone-active');
    });

    dropZone.addEventListener('drop', (e) => {
        e.preventDefault();
        dropZone.style.borderColor = 'rgba(42, 82, 152, 0.3)';
        dropZone.style.backgroundColor = '';
        dropZone.classList.remove('upload-zone-active');
        
        const file = e.dataTransfer.files[0];
        if (file) {
            handleFile(file);
        }
    });

    // File input change handler
    fileInput.addEventListener('change', (e) => {
        const file = e.target.files[0];
        if (file) {
            handleFile(file);
        }
    });

    // Sequence input handling
    sequenceInput.addEventListener('input', () => {
        sequenceInput.classList.remove('is-invalid');
    });

    // Analyze Button Click Handler with improved error handling
    analyzeButton.addEventListener('click', async () => {
        const sequenceId = sequenceInput.value.trim();
        
        if (!validateSequenceId(sequenceId)) {
            return;
        }

        showLoading();
        analyzeButton.disabled = true;

        try {
            const response = await fetch('/api/analyze', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    'Accept': 'application/json'
                },
                body: JSON.stringify({ sequence_id: sequenceId })
            });

            if (!response.ok) {
                const errorData = await response.json();
                throw new Error(errorData.error || 'Failed to analyze sequence');
            }

            const data = await response.json();
            currentResults = data;
            processResults(data);
            showToast('Sequence analysis completed successfully!', 'success');
        } catch (error) {
            console.error('Error:', error);
            showToast(error.message || 'Failed to analyze sequence');
            resetResults();
        } finally {
            hideLoading();
            analyzeButton.disabled = false;
        }
    });

    // File Handler with improved error handling
    async function handleFile(file) {
        if (!file.name.match(/\.(fasta|fa|txt)$/i)) {
            showToast('Please upload a valid FASTA file (.fasta, .fa, or .txt)');
            return;
        }

        showLoading();
        analyzeButton.disabled = true;

        const formData = new FormData();
        formData.append('file', file);

        try {
            const response = await fetch('/api/analyze', {
                method: 'POST',
                body: formData
            });

            if (!response.ok) {
                const errorData = await response.json();
                throw new Error(errorData.error || 'Failed to analyze sequence');
            }

            const data = await response.json();
            currentResults = data;
            processResults(data);
            showToast('File analysis completed successfully!', 'success');
        } catch (error) {
            console.error('Error:', error);
            showToast(error.message || 'Failed to analyze sequence');
            resetResults();
        } finally {
            hideLoading();
            analyzeButton.disabled = false;
        }
    }

    // Function to colorize nucleotides
    function colorizeSequence(sequence) {
        const colors = {
            'A': '#FF6B6B', // Red
            'T': '#4ECDC4', // Teal
            'G': '#45B7D1', // Blue
            'C': '#96CEB4', // Green
            'U': '#FFD93D'  // Yellow for RNA
        };
        
        return sequence.split('').map(nucleotide => {
            const color = colors[nucleotide.toUpperCase()] || '#212529';
            return `<span style="color: ${color}; font-weight: 500;">${nucleotide}</span>`;
        }).join('');
    }

    // Process and Display Results with animation
    function processResults(data) {
        resultsContainer.classList.remove('d-none');
        
        // Update result displays with colorized sequences
        const results = [
            { element: dnaResult, data: data.dna_sequence, fallback: 'No DNA sequence available' },
            { element: rnaResult, data: data.rna_sequence, fallback: 'No RNA sequence available' },
            { element: proteinResult, data: data.protein_sequence, fallback: 'No protein sequence available' }
        ];

        results.forEach(({ element, data, fallback }) => {
            element.style.opacity = '0';
            if (data) {
                // Only colorize DNA and RNA sequences, not protein
                if (element === proteinResult) {
                    element.innerHTML = data;
                } else {
                    element.innerHTML = colorizeSequence(data);
                }
            } else {
                element.textContent = fallback;
            }
            element.style.transition = 'opacity 0.3s ease-in-out';
            setTimeout(() => element.style.opacity = '1', 50);
        });

        // Update sequence statistics
        sequenceLength.textContent = `${data.length || 0} bp`;
        gcContent.textContent = `${data.gc_content || 0}%`;

        // Enable download button
        downloadButton.disabled = false;

        // Smooth scroll to results
        resultsContainer.scrollIntoView({ behavior: 'smooth', block: 'start' });
    }

    // Reset results container
    function resetResults() {
        resultsContainer.classList.add('d-none');
        [dnaResult, rnaResult, proteinResult].forEach(element => {
            element.textContent = '';
        });
        sequenceLength.textContent = '-';
        gcContent.textContent = '-';
        downloadButton.disabled = true;
        currentResults = null;
    }

    // Loading State Handlers with fade effect
    function showLoading() {
        loadingOverlay.style.display = 'block';
        setTimeout(() => loadingOverlay.classList.remove('d-none'), 0);
    }

    function hideLoading() {
        loadingOverlay.classList.add('d-none');
        setTimeout(() => loadingOverlay.style.display = 'none', 300);
    }
}); 