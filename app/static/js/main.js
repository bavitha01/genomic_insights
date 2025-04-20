document.addEventListener('DOMContentLoaded', function() {
    // Get form elements
    const fastaUploadForm = document.getElementById('fastaUploadForm');
    const sequenceSearchForm = document.getElementById('sequenceSearchForm');
    const loadingOverlay = document.getElementById('loadingOverlay');
    const resultsContainer = document.getElementById('resultsContainer');

    // Handle FASTA file upload
    fastaUploadForm.addEventListener('submit', function(e) {
        e.preventDefault();
        
        const fileInput = document.getElementById('fastaFile');
        if (!fileInput.files.length) {
            alert('Please select a FASTA file to upload.');
            return;
        }

        const formData = new FormData();
        formData.append('fastaFile', fileInput.files[0]);

        showLoading();

        fetch('/upload', {
            method: 'POST',
            body: formData
        })
        .then(response => {
            if (!response.ok) {
                return response.json().then(data => {
                    throw new Error(data.error || 'Error uploading file');
                });
            }
            return response.json();
        })
        .then(data => {
            hideLoading();
            if (data.results && data.results.length > 0) {
                // Check if any results have errors
                const hasErrors = data.results.some(result => result.error);
                if (hasErrors) {
                    const errorMessages = data.results
                        .filter(result => result.error)
                        .map(result => `${result.seq_id}: ${result.error}`)
                        .join('\n');
                    
                    // If all results have errors
                    if (data.results.every(result => result.error)) {
                        alert(`All sequences had errors:\n${errorMessages}`);
                        return;
                    } else {
                        // Some sequences had errors, but not all
                        console.warn(`Some sequences had errors:\n${errorMessages}`);
                        // Continue to display the valid ones
                        displayResults(data.results.filter(result => !result.error));
                    }
                } else {
                    displayResults(data.results);
                }
            } else {
                alert('No valid sequences found in the file.');
            }
        })
        .catch(error => {
            hideLoading();
            alert('Error: ' + error.message);
            console.error('Error:', error);
        });
    });

    // Handle sequence ID search
    sequenceSearchForm.addEventListener('submit', function(e) {
        e.preventDefault();
        
        const sequenceId = document.getElementById('sequenceId').value.trim();
        if (!sequenceId) {
            alert('Please enter a sequence ID.');
            return;
        }

        showLoading();

        fetch('/search', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({ seq_id: sequenceId })
        })
        .then(response => {
            if (!response.ok) {
                return response.json().then(data => {
                    throw new Error(data.error || 'Error searching for sequence');
                });
            }
            return response.json();
        })
        .then(data => {
            hideLoading();
            if (data.result) {
                // Display a single result in an array for consistency with the upload results
                displayResults([data.result]);
            } else {
                alert('No sequence found with that ID.');
            }
        })
        .catch(error => {
            hideLoading();
            alert('Error: ' + error.message);
            console.error('Error:', error);
        });
    });

    // Show loading overlay
    function showLoading() {
        loadingOverlay.classList.remove('d-none');
    }

    // Hide loading overlay
    function hideLoading() {
        loadingOverlay.classList.add('d-none');
    }

    // Display results
    function displayResults(results) {
        if (!results || results.length === 0) {
            alert('No results to display.');
            return;
        }

        // Show results container
        resultsContainer.classList.remove('d-none');

        // Format and display summary
        const summaryContent = document.getElementById('summaryContent');
        let summaryHTML = '<div class="table-responsive"><table class="table table-bordered table-results">';
        summaryHTML += '<thead><tr><th>Sequence ID</th><th>Length</th><th>GC Content</th></tr></thead><tbody>';
        
        results.forEach(result => {
            summaryHTML += `
                <tr>
                    <td>${result.seq_id}</td>
                    <td>${result.length} bp</td>
                    <td>
                        ${result.gc_content}%
                        <div class="gc-content-bar">
                            <div class="gc-content-fill" style="width: ${result.gc_content}%"></div>
                        </div>
                    </td>
                </tr>
            `;
        });
        
        summaryHTML += '</tbody></table></div>';
        
        // Add GenBank metadata if available
        if (results[0].description || results[0].organism || results[0].taxonomy) {
            summaryHTML += `<div class="mt-4 card">
                <div class="card-header bg-light">
                    <h5 class="mb-0">GenBank Metadata</h5>
                </div>
                <div class="card-body">`;
                
            if (results[0].description) {
                summaryHTML += `
                    <div class="mb-2">
                        <strong>Description:</strong>
                        <p>${results[0].description}</p>
                    </div>`;
            }
            
            if (results[0].organism) {
                summaryHTML += `
                    <div class="mb-2">
                        <strong>Organism:</strong>
                        <p>${results[0].organism}</p>
                    </div>`;
            }
            
            if (results[0].taxonomy) {
                summaryHTML += `
                    <div class="mb-2">
                        <strong>Taxonomy:</strong>
                        <p>${results[0].taxonomy.join(' > ')}</p>
                    </div>`;
            }
            
            summaryHTML += `</div></div>`;
        }
        
        summaryContent.innerHTML = summaryHTML;

        // Format and display sequence data
        displaySequenceData('sequenceContent', results, 'sequence');
        displaySequenceData('transcriptionContent', results, 'transcription');
        displaySequenceData('translationContent', results, 'translation');

        // Scroll to results
        resultsContainer.scrollIntoView({ behavior: 'smooth' });
    }

    // Display sequence data with formatting
    function displaySequenceData(elementId, results, dataType) {
        const element = document.getElementById(elementId);
        let html = '';

        results.forEach(result => {
            if (result[dataType]) {
                const sequenceData = result[dataType];
                html += `<h5>${result.seq_id}</h5>`;
                
                // Add download button for the sequence
                html += `<button class="btn btn-sm btn-outline-primary mb-2" 
                            onclick="downloadSequence('${result.seq_id}', '${dataType}', \`${sequenceData}\`)">
                            <i class="bi bi-download"></i> Download ${dataType}
                        </button>`;
                
                html += `<div class="sequence-data">${formatSequence(sequenceData, dataType)}</div>`;
                if (results.length > 1) {
                    html += '<hr>';
                }
            }
        });

        element.innerHTML = html;
    }

    // Format sequence with colored nucleotides
    function formatSequence(sequence, type) {
        let formattedSequence = '';
        const chunkSize = 70; // Characters per line
        
        for (let i = 0; i < sequence.length; i += chunkSize) {
            const chunk = sequence.substring(i, i + chunkSize);
            const positionNumber = i + 1;
            
            // Add position number at the beginning of each line
            formattedSequence += `<span class="position-number">${positionNumber.toString().padStart(8, ' ')}</span> `;
            
            // Format each nucleotide/amino acid with color
            for (let j = 0; j < chunk.length; j++) {
                const char = chunk[j].toUpperCase();
                if (type === 'sequence' || type === 'transcription') {
                    // For DNA or RNA sequences
                    const nucleotideClass = `nucleotide-${char.toLowerCase()}`;
                    formattedSequence += `<span class="${nucleotideClass}">${char}</span>`;
                } else {
                    // For protein sequences (translation)
                    formattedSequence += char;
                }
            }
            
            formattedSequence += '<br>';
        }
        
        return formattedSequence;
    }
});

// Function to download sequence data
function downloadSequence(seqId, dataType, sequence) {
    // Create file content
    const content = `>${seqId} (${dataType})\n${sequence}`;
    
    // Create a blob
    const blob = new Blob([content], { type: 'text/plain' });
    
    // Create a download link
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = `${seqId}_${dataType}.fasta`;
    
    // Trigger download
    document.body.appendChild(a);
    a.click();
    
    // Clean up
    document.body.removeChild(a);
    URL.revokeObjectURL(a.href);
} 