document.addEventListener('DOMContentLoaded', function() {
    // Get form elements
    const fastaUploadForm = document.getElementById('fastaUploadForm');
    const sequenceSearchForm = document.getElementById('sequenceSearchForm');
    const genbankCompareForm = document.getElementById('genbankCompareForm');
    const fastaCompareForm = document.getElementById('fastaCompareForm');
    const addGenbankIdBtn = document.getElementById('addGenbankIdBtn');
    const loadingOverlay = document.getElementById('loadingOverlay');
    const resultsContainer = document.getElementById('resultsContainer');
    const comparisonResultsContainer = document.getElementById('comparisonResultsContainer');
    const showComparisonBtn = document.getElementById('showComparisonBtn');
    const showRegularResultsBtn = document.getElementById('showRegularResultsBtn');

    // Chart references
    let gcContentChart = null;
    
    // Store results globally within this scope
    let lastAnalysisResults = null;
    let lastComparisonResults = null;

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

    // Add button for GenBank IDs
    addGenbankIdBtn.addEventListener('click', function() {
        const container = document.getElementById('genbankIdContainer');
        const inputDiv = document.createElement('div');
        inputDiv.className = 'genbank-id-input mb-2';
        inputDiv.innerHTML = `
            <div class="input-group">
                <input type="text" class="form-control" name="sequenceIds[]" placeholder="Enter sequence ID">
                <button type="button" class="btn btn-outline-danger remove-id-btn">
                    <i class="bi bi-x"></i>
                </button>
            </div>
        `;
        container.appendChild(inputDiv);
        
        // Add event listener for the remove button
        inputDiv.querySelector('.remove-id-btn').addEventListener('click', function() {
            container.removeChild(inputDiv);
        });
    });

    // Handle GenBank sequence comparison form
    genbankCompareForm.addEventListener('submit', function(e) {
        e.preventDefault();
        
        // Get all sequence IDs
        const inputs = this.querySelectorAll('input[name="sequenceIds[]"]');
        const sequenceIds = Array.from(inputs).map(input => input.value.trim()).filter(id => id !== '');
        
        if (sequenceIds.length < 2) {
            alert('Please enter at least two sequence IDs for comparison.');
            return;
        }
        
        showLoading();
        
        fetch('/compare', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({ sequence_ids: sequenceIds })
        })
        .then(response => {
            if (!response.ok) {
                return response.json().then(data => {
                    throw new Error(data.error || 'Error comparing sequences');
                });
            }
            return response.json();
        })
        .then(data => {
            hideLoading();
            if (data.result) {
                displayComparisonResults(data.result);
            }
        })
        .catch(error => {
            hideLoading();
            alert('Error: ' + error.message);
            console.error('Error:', error);
        });
    });
    
    // Handle FASTA file comparison form
    fastaCompareForm.addEventListener('submit', function(e) {
        e.preventDefault();
        
        const fileInput = document.getElementById('compareFile');
        if (!fileInput.files.length) {
            alert('Please select a FASTA file containing multiple sequences.');
            return;
        }
        
        showLoading();
        
        // First upload the file
        const formData = new FormData();
        formData.append('fastaFile', fileInput.files[0]);
        
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
            if (data.results && data.results.length >= 2) {
                // Now that we have results, request comparison
                return fetch('/compare', {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json'
                    },
                    body: JSON.stringify({ file_id: fileInput.files[0].name })
                });
            } else {
                throw new Error('Not enough valid sequences found in the file (minimum 2 required).');
            }
        })
        .then(response => {
            if (!response.ok) {
                return response.json().then(data => {
                    throw new Error(data.error || 'Error comparing sequences');
                });
            }
            return response.json();
        })
        .then(data => {
            hideLoading();
            if (data.result) {
                displayComparisonResults(data.result);
            }
        })
        .catch(error => {
            hideLoading();
            alert('Error: ' + error.message);
            console.error('Error:', error);
        });
    });

    // Toggle between results views
    showComparisonBtn.addEventListener('click', function() {
        if (lastComparisonResults) {
            resultsContainer.classList.add('d-none');
            comparisonResultsContainer.classList.remove('d-none');
            displayComparisonResults(lastComparisonResults);
        } else {
            alert('No comparison results available. Please perform a comparison first.');
        }
    });
    
    showRegularResultsBtn.addEventListener('click', function() {
        if (lastAnalysisResults) {
            comparisonResultsContainer.classList.add('d-none');
            resultsContainer.classList.remove('d-none');
        } else {
            alert('No analysis results available. Please analyze a sequence first.');
        }
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

        // Store results globally
        lastAnalysisResults = results;
        
        // If we have multiple sequences, show the comparison button
        if (results.length > 1) {
            showComparisonBtn.classList.remove('d-none');
            
            // Create comparison data and store it
            const comparisonData = {
                sequences: results.map(result => ({
                    seq_id: result.seq_id,
                    length: result.length,
                    gc_content: result.gc_content
                })),
                gc_content: {
                    labels: results.map(result => result.seq_id),
                    values: results.map(result => result.gc_content)
                }
            };
            lastComparisonResults = comparisonData;
        } else {
            showComparisonBtn.classList.add('d-none');
        }

        // Show results container
        resultsContainer.classList.remove('d-none');
        comparisonResultsContainer.classList.add('d-none');

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

    // Display comparison results
    function displayComparisonResults(results) {
        // Store results globally
        lastComparisonResults = results;
        
        // Show comparison results container
        comparisonResultsContainer.classList.remove('d-none');
        resultsContainer.classList.add('d-none');
        
        // Clear existing chart if it exists
        if (gcContentChart) {
            gcContentChart.destroy();
        }
        
        // Create GC content chart
        const gcContentCtx = document.getElementById('gcContentChart').getContext('2d');
        gcContentChart = new Chart(gcContentCtx, {
            type: 'bar',
            data: {
                labels: results.gc_content.labels,
                datasets: [{
                    label: 'GC Content (%)',
                    data: results.gc_content.values,
                    backgroundColor: 'rgba(75, 192, 192, 0.6)',
                    borderColor: 'rgba(75, 192, 192, 1)',
                    borderWidth: 1
                }]
            },
            options: {
                responsive: true,
                scales: {
                    y: {
                        beginAtZero: true,
                        max: 100,
                        title: {
                            display: true,
                            text: 'GC Content (%)'
                        }
                    },
                    x: {
                        title: {
                            display: true,
                            text: 'Sequence ID'
                        }
                    }
                },
                plugins: {
                    title: {
                        display: true,
                        text: 'GC Content Comparison'
                    }
                }
            }
        });
        
        // Fill the comparison details table
        const tableBody = document.getElementById('comparisonTableBody');
        tableBody.innerHTML = '';
        
        results.sequences.forEach(seq => {
            const row = document.createElement('tr');
            row.innerHTML = `
                <td>${seq.seq_id}</td>
                <td>${seq.length.toLocaleString()} bp</td>
                <td>
                    ${seq.gc_content}%
                    <div class="gc-content-bar">
                        <div class="gc-content-fill" style="width: ${seq.gc_content}%"></div>
                    </div>
                </td>
            `;
            tableBody.appendChild(row);
        });
        
        // Scroll to comparison results
        comparisonResultsContainer.scrollIntoView({ behavior: 'smooth' });
    }
});

// Function to download sequence data (needs to be in global scope)
function downloadSequence(seqId, dataType, sequence) {
    // Create file name
    const fileName = `${seqId}_${dataType}.txt`;
    
    // Create blob
    const blob = new Blob([sequence], { type: 'text/plain' });
    
    // Create download link
    const a = document.createElement('a');
    a.href = window.URL.createObjectURL(blob);
    a.download = fileName;
    
    // Trigger download
    document.body.appendChild(a);
    a.click();
    
    // Cleanup
    window.URL.revokeObjectURL(a.href);
    document.body.removeChild(a);
} 