document.addEventListener('DOMContentLoaded', function() {
    // Get UI elements
    const inputMethods = document.getElementById('inputMethods');
    const resultsTable = document.getElementById('resultsTable');
    const addIdField = document.getElementById('addIdField');
    const loadingOverlay = document.getElementById('loadingOverlay');
    const comparisonResultsContainer = document.getElementById('comparisonResultsContainer');

    // Initialize Bootstrap tabs
    const tabElements = document.querySelectorAll('[data-bs-toggle="tab"]');
    tabElements.forEach(tab => {
        tab.addEventListener('shown.bs.tab', function(event) {
            // Clear any error messages when switching tabs
            const forms = document.querySelectorAll('form');
            forms.forEach(form => form.reset());
        });
    });

    // Check for stored sequences first
    try {
        const storedData = sessionStorage.getItem('compareSequences');
        if (storedData) {
            const data = JSON.parse(storedData);
            if (data && data.sequences && data.gc_content && data.sequences.length >= 1) {
                // Hide input methods if we have stored data
                inputMethods.classList.add('d-none');
                resultsTable.classList.remove('d-none');
                
                // Display the data
                populateOverviewTable(data.sequences);
                createGCContentChart(data.gc_content);
                displayStatistics(data.sequences);
                
                // Clear session storage
                sessionStorage.removeItem('compareSequences');
                return;
            }
        }
        
        // If we get here, show the input methods
        inputMethods.classList.remove('d-none');
        resultsTable.classList.add('d-none');
        
    } catch (error) {
        console.error('Error initializing comparison:', error);
        showError('Error loading comparison data. Please try again.');
    }

    // Add ID field button handler
    if (addIdField) {
        addIdField.addEventListener('click', function() {
            const container = document.getElementById('sequenceIdInputs');
            const inputGroup = document.createElement('div');
            inputGroup.className = 'input-group mb-3';
            inputGroup.innerHTML = `
                <input type="text" class="form-control" name="sequenceIds[]" placeholder="Enter sequence ID">
                <button type="button" class="btn btn-outline-danger" onclick="this.parentElement.remove()">
                    <i class="bi bi-x"></i>
                </button>
            `;
            container.appendChild(inputGroup);
        });
    }

    // Handle sequence ID comparison form
    const idCompareForm = document.getElementById('idCompareForm');
    if (idCompareForm) {
        idCompareForm.addEventListener('submit', async function(e) {
            e.preventDefault();
            
            const inputs = this.querySelectorAll('input[name="sequenceIds[]"]');
            const sequenceIds = Array.from(inputs).map(input => input.value.trim()).filter(Boolean);
            
            if (sequenceIds.length < 2) {
                alert('Please enter at least two sequence IDs');
                return;
            }

            try {
                const response = await fetch('/api/compare', {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json',
                        'Accept': 'application/json'
                    },
                    body: JSON.stringify({ sequence_ids: sequenceIds })
                });

                if (!response.ok) {
                    const errorData = await response.json();
                    throw new Error(errorData.error || 'Failed to compare sequences');
                }

                const data = await response.json();
                displayResults(data);
            } catch (error) {
                console.error('Error:', error);
                alert(error.message || 'Failed to compare sequences');
            }
        });
    }

    // Handle FASTA file comparison form
    const fastaCompareForm = document.getElementById('fastaCompareForm');
    if (fastaCompareForm) {
        fastaCompareForm.addEventListener('submit', async function(e) {
            e.preventDefault();
            
            const fileInput = document.getElementById('fastaFile');
            const file = fileInput.files[0];
            
            if (!file) {
                alert('Please select a FASTA file');
                return;
            }

            const formData = new FormData();
            formData.append('file', file);

            try {
                const response = await fetch('/api/compare', {
                    method: 'POST',
                    body: formData
                });

                if (!response.ok) {
                    const errorData = await response.json();
                    throw new Error(errorData.error || 'Failed to compare sequences');
                }

                const data = await response.json();
                displayResults(data);
            } catch (error) {
                console.error('Error:', error);
                alert(error.message || 'Failed to compare sequences');
            }
        });
    }

    // Helper functions
    function showLoading() {
        if (loadingOverlay) {
            loadingOverlay.classList.remove('d-none');
        }
    }

    function hideLoading() {
        if (loadingOverlay) {
            loadingOverlay.classList.add('d-none');
        }
    }

    function handleResponse(response) {
        if (!response.ok) {
            return response.json().then(data => {
                throw new Error(data.error || 'Error processing request');
            });
        }
        return response.json();
    }

    function handleError(error) {
        hideLoading();
        alert('Error: ' + error.message);
        console.error('Error:', error);
    }

    function displayResults(data) {
        // Hide input methods and show results
        inputMethods.classList.add('d-none');
        resultsTable.classList.remove('d-none');
        
        // Display the data
        populateOverviewTable(data.sequences);
        createGCContentChart(data.gc_content);
        displayStatistics(data.sequences);
        
        // Scroll to results
        resultsTable.scrollIntoView({ behavior: 'smooth' });
    }
});

function populateOverviewTable(sequences) {
    const tbody = document.getElementById('overviewTableBody');
    if (!tbody) return;

    tbody.innerHTML = sequences.map(seq => `
        <tr>
            <td>${seq.seq_id}</td>
            <td>${seq.length.toLocaleString()} bp</td>
            <td>
                <div class="d-flex align-items-center">
                    <div class="flex-grow-1 me-2">
                        <div class="progress" style="height: 8px;">
                            <div class="progress-bar" role="progressbar" 
                                style="width: ${seq.gc_content}%; background: linear-gradient(90deg, #3a86ff, #60a5fa);">
                            </div>
                        </div>
                    </div>
                    <span class="text-nowrap">${seq.gc_content.toFixed(2)}%</span>
                </div>
            </td>
        </tr>
    `).join('');
}

function createGCContentChart(data) {
    const ctx = document.getElementById('gcContentChart');
    if (!ctx) return;

    // Destroy existing chart if it exists
    if (window.gcChart) {
        window.gcChart.destroy();
    }

    const gradient = ctx.getContext('2d').createLinearGradient(0, 0, 0, 400);
    gradient.addColorStop(0, '#3a86ff');
    gradient.addColorStop(1, '#60a5fa');

    window.gcChart = new Chart(ctx, {
        type: 'bar',
        data: {
            labels: data.labels,
            datasets: [{
                label: 'GC Content (%)',
                data: data.values,
                backgroundColor: gradient,
                borderColor: '#3a86ff',
                borderWidth: 0,
                borderRadius: 8,
                maxBarThickness: 60
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: {
                legend: {
                    display: false
                },
                tooltip: {
                    callbacks: {
                        label: function(context) {
                            return `GC Content: ${context.raw.toFixed(2)}%`;
                        }
                    }
                }
            },
            scales: {
                y: {
                    beginAtZero: true,
                    max: 100,
                    grid: {
                        color: 'rgba(0, 0, 0, 0.05)'
                    },
                    ticks: {
                        callback: function(value) {
                            return value + '%';
                        },
                        color: '#666',
                        font: {
                            family: "'Roboto', sans-serif",
                            size: 12
                        }
                    }
                },
                x: {
                    grid: {
                        display: false
                    },
                    ticks: {
                        color: '#666',
                        font: {
                            family: "'Roboto', sans-serif",
                            size: 12
                        }
                    }
                }
            }
        }
    });
}

function displayStatistics(sequences) {
    const statsContainer = document.getElementById('gcStats');
    if (!statsContainer) return;

    // Calculate statistics
    const gcValues = sequences.map(seq => seq.gc_content);
    const stats = {
        average: calculateAverage(gcValues),
        median: calculateMedian(gcValues),
        min: Math.min(...gcValues),
        max: Math.max(...gcValues),
        range: Math.max(...gcValues) - Math.min(...gcValues)
    };

    // Display statistics
    statsContainer.innerHTML = `
        <div class="stat-item mb-3 p-3 bg-white rounded-3 shadow-sm">
            <div class="text-muted mb-1">Average GC Content</div>
            <div class="h4 mb-0">${stats.average.toFixed(2)}%</div>
        </div>
        <div class="stat-item mb-3 p-3 bg-white rounded-3 shadow-sm">
            <div class="text-muted mb-1">Median GC Content</div>
            <div class="h4 mb-0">${stats.median.toFixed(2)}%</div>
        </div>
        <div class="stat-item mb-3 p-3 bg-white rounded-3 shadow-sm">
            <div class="text-muted mb-1">Range</div>
            <div class="h4 mb-0">${stats.range.toFixed(2)}%</div>
        </div>
        <div class="stat-item mb-3 p-3 bg-white rounded-3 shadow-sm">
            <div class="text-muted mb-1">Minimum</div>
            <div class="h4 mb-0">${stats.min.toFixed(2)}%</div>
        </div>
        <div class="stat-item mb-3 p-3 bg-white rounded-3 shadow-sm">
            <div class="text-muted mb-1">Maximum</div>
            <div class="h4 mb-0">${stats.max.toFixed(2)}%</div>
        </div>
    `;
}

function calculateAverage(values) {
    return values.reduce((a, b) => a + b, 0) / values.length;
}

function calculateMedian(values) {
    const sorted = [...values].sort((a, b) => a - b);
    const middle = Math.floor(sorted.length / 2);
    
    if (sorted.length % 2 === 0) {
        return (sorted[middle - 1] + sorted[middle]) / 2;
    }
    
    return sorted[middle];
}

function showError(message) {
    const main = document.querySelector('main');
    if (!main) return;

    main.innerHTML = `
        <div class="alert alert-danger" role="alert">
            <i class="bi bi-exclamation-triangle me-2"></i>
            ${message}
        </div>
        <div class="text-center mt-4">
            <a href="/analyze" class="btn btn-primary">
                <i class="bi bi-arrow-left me-2"></i>Back to Analysis
            </a>
        </div>
    `;
} 