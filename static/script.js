document.addEventListener('DOMContentLoaded', () => {
    // --- MODAL LOGIC ---
    const modal = document.getElementById('infoModal');
    const infoBtn = document.getElementById('infoBtn');
    const closeBtn = document.querySelector('.close-btn');

    if(infoBtn) infoBtn.onclick = () => modal.style.display = "block";
    if(closeBtn) closeBtn.onclick = () => modal.style.display = "none";
    window.onclick = (e) => { if(e.target == modal) modal.style.display = "none"; }

    let globalData = []; // Store data for CSV export
    let gaugeChartInstance = null;
    let radarChartInstance = null;
    let trendChartInstance = null;

    // --- UPLOAD LOGIC ---
    const dropZone = document.getElementById('drop-zone');
    const fileInput = document.getElementById('file-input');
    const fileList = document.getElementById('file-list');
    const calculateBtn = document.getElementById('calculate-btn');
    
    const uploadPanel = document.getElementById('upload-panel');
    const progressPanel = document.getElementById('progress-panel');
    const resultsPanel = document.getElementById('results-panel');
    
    let selectedFiles = [];

    // Drag & Drop Handlers
    dropZone.addEventListener('dragover', (e) => { e.preventDefault(); dropZone.classList.add('dragover'); });
    dropZone.addEventListener('dragleave', () => dropZone.classList.remove('dragover'));
    dropZone.addEventListener('drop', (e) => {
        e.preventDefault(); dropZone.classList.remove('dragover');
        handleFiles(e.dataTransfer.files);
    });

    fileInput.addEventListener('change', (e) => handleFiles(e.target.files));

    function handleFiles(files) {
        for (let file of files) {
            if (file.name.endsWith('.fasta') || file.name.endsWith('.fa')) {
                selectedFiles.push(file);
                const badge = document.createElement('div');
                badge.className = 'file-badge';
                badge.textContent = file.name;
                fileList.appendChild(badge);
            }
        }
        if (selectedFiles.length > 0) {
            calculateBtn.style.display = 'block';
            const weightPanel = document.getElementById('weight-config-panel');
            if (weightPanel) weightPanel.style.display = 'block';
        }
    }

    calculateBtn.addEventListener('click', async () => {
        if (selectedFiles.length === 0) return;
        
        const formData = new FormData();
        selectedFiles.forEach(file => formData.append('files[]', file));
        
        // Harvest custom weights
        const customWeights = {};
        const weightIds = ['w_mu', 'w_re', 'w_pi', 'w_mb', 'w_dnds', 'w_gd', 'w_ri', 'w_cai'];
        weightIds.forEach(id => {
            const input = document.getElementById(id);
            if (input && input.value !== "") {
                customWeights[id] = parseFloat(input.value);
            }
        });
        formData.append('weights', JSON.stringify(customWeights));
        
        uploadPanel.style.display = 'none';
        progressPanel.style.display = 'block';
        
        try {
            const res = await fetch('/upload', { method: 'POST', body: formData });
            const data = await res.json();
            if (data.job_id) pollStatus(data.job_id);
            else alert(data.error);
        } catch (err) {
            alert('Server error while uploading.');
        }
    });

    async function pollStatus(jobId) {
        const pBar = document.getElementById('progress-bar');
        const pPercent = document.getElementById('progress-percent');
        const cFile = document.getElementById('current-file-lbl');
        const cStep = document.getElementById('current-step-lbl');
        
        const interval = setInterval(async () => {
            try {
                const res = await fetch(`/status/${jobId}`);
                if (!res.ok) return; // Skip if server temporarily unavail
                const status = await res.json();
                
                pBar.style.width = `${status.progress}%`;
                pPercent.textContent = `${status.progress}%`;
                cFile.textContent = `Processing: ${status.current_file}`;
                
                if (cStep && status.current_step) {
                    cStep.textContent = status.current_step;
                }
                
                if (status.status === 'completed') {
                    clearInterval(interval);
                    progressPanel.style.display = 'none';
                    resultsPanel.style.display = 'block';
                    
                    globalData = status.results; // Save for CSV export
                    renderResults(status.results);
                } else if (status.status === 'error') {
                    clearInterval(interval);
                    alert('Error: ' + status.message);
                }
            } catch (err) {
                console.warn('Polling interrupted (e.g. network change). Retrying in 1s...', err);
            }
        }, 1000);
    }

    function renderResults(data) {
        // Render GVI 2.0 Score Cards (Use the most recent/latest row in dataset)
        if (data.length > 0) {
            const latestRow = data[data.length - 1];
            if (latestRow.GVI_Score !== undefined) {
                const gviDisplay = document.getElementById('gvi-score-display');
                if (gviDisplay) gviDisplay.textContent = latestRow.GVI_Score;
            }
        }

        // Render Table
        const tbody = document.getElementById('results-body');
        const years = [], dnds = [], pi = [], mut = [], evo = [];
        
        data.forEach(row => {
            const tr = document.createElement('tr');
            tr.innerHTML = `
                <td>${row.Year}</td>
                <td>${row.Sequence_Count}</td>
                <td>${row.Detected_Model}</td>
                <td>${row.Nucleotide_Diversity_Pi}</td>
                <td>${row.Genetic_Distance_GD}</td>
                <td>${row.dN_dS_Ratio}</td>
                <td>${row.Mutation_Burden_MB}</td>
                <td>${row.Codon_Adaptation_Index}</td>
                <td>${row.Effective_Reproduction_Number_Re !== undefined ? row.Effective_Reproduction_Number_Re : '-'}</td>
                <td>${row.GC_Content}</td>
                <td>${row.Evolutionary_Rate}</td>
                <td>${row.Recombination_Rate}</td>
                <td style="font-weight: 700; color: #0f4c81;">${row.GVI_Score !== undefined ? row.GVI_Score : '-'}</td>
            `;
            tbody.appendChild(tr);
            
            years.push(row.Year);
            dnds.push(row.dN_dS_Ratio);
            pi.push(row.Nucleotide_Diversity_Pi);
            mut.push(row.Mutation_Burden_MB);
            evo.push(parseFloat(row.Evolutionary_Rate));
        });

        // Render Tree Gallery (Interactive Phylocanvas)
        const gallery = document.getElementById('tree-gallery');
        gallery.innerHTML = '';
        
        data.forEach(row => {
            if (row.Newick_Tree && row.Newick_Tree.length > 5) {
                const card = document.createElement('div');
                card.className = 'tree-card';
                
                const canvasId = `tree-canvas-${row.Year}`;
                card.innerHTML = `
                    <div id="${canvasId}" class="tree-canvas"></div>
                    <h4>Year ${row.Year}</h4>
                `;
                gallery.appendChild(card);
                
                // Initialize Phylocanvas inside the div
                setTimeout(() => {
                    console.log("Tree string:", row.Newick_Tree);
                    try {
                        const createFn = Phylocanvas.createTree || (Phylocanvas.default && Phylocanvas.default.createTree);
                        if (!createFn) throw new Error("Could not find createTree function in Phylocanvas");
                        const tree = createFn(canvasId);
                        tree.load(row.Newick_Tree);
                        tree.setTreeType('radial');
                        tree.setNodeSize(3);
                        tree.setTextSize(12);
                        tree.draw();
                    } catch (e) {
                        console.error("Phylocanvas failed to load:", e);
                        document.getElementById(canvasId).innerHTML = "<p style='color:red; margin-top:20px;'>Error rendering interactive tree: " + (e.message || e.toString()) + "</p>";
                    }
                }, 200);
            }
        });

        if (gaugeChartInstance) gaugeChartInstance.destroy();
        if (radarChartInstance) radarChartInstance.destroy();
        if (trendChartInstance) trendChartInstance.destroy();

        let latestRow = data[data.length - 1];
        let score = latestRow.GVI_Score;
        if (score === "Not Calculated" || score === undefined) score = 0;
        
        // 1. Gauge Chart
        const ctxGauge = document.getElementById('gaugeChart').getContext('2d');
        gaugeChartInstance = new Chart(ctxGauge, {
            type: 'doughnut',
            data: {
                labels: ['Score', 'Remaining'],
                datasets: [{
                    data: [score, 100 - score],
                    backgroundColor: ['#3498db', '#222'],
                    borderWidth: 0
                }]
            },
            options: {
                circumference: 180,
                rotation: 270,
                cutout: '80%',
                responsive: true,
                maintainAspectRatio: false,
                plugins: { legend: { display: false }, tooltip: { enabled: false } }
            }
        });

        // 2. Radar Chart
        const ctxRadar = document.getElementById('radarChart').getContext('2d');
        const parseVal = (v) => isNaN(parseFloat(v)) ? 0 : parseFloat(v);
        radarChartInstance = new Chart(ctxRadar, {
            type: 'radar',
            data: {
                labels: ['Pi (x10)', 'GD (x10)', 'dN/dS', 'MB (/100)', 'CAI', 'Re', 'GC', 'Recomb'],
                datasets: [{
                    label: 'Latest Dataset Metrics',
                    data: [
                        parseVal(latestRow.Nucleotide_Diversity_Pi) * 10,
                        parseVal(latestRow.Genetic_Distance_GD) * 10,
                        parseVal(latestRow.dN_dS_Ratio),
                        parseVal(latestRow.Mutation_Burden_MB) / 100,
                        parseVal(latestRow.Codon_Adaptation_Index),
                        parseVal(latestRow.Effective_Reproduction_Number_Re),
                        parseVal(latestRow.GC_Content),
                        parseVal(latestRow.Recombination_Rate)
                    ],
                    backgroundColor: 'rgba(52, 152, 219, 0.4)',
                    borderColor: 'rgba(52, 152, 219, 1)',
                    pointBackgroundColor: '#fff',
                }]
            },
            options: {
                responsive: true, maintainAspectRatio: false,
                scales: { r: { ticks: { display: false }, grid: { color: '#444' }, angleLines: { color: '#444' }, pointLabels: { color: '#bbb' } } },
                plugins: { legend: { display: false } }
            }
        });

        // 3. Trend Chart
        const ctxTrend = document.getElementById('trendChart').getContext('2d');
        const validData = data.filter(r => r.GVI_Score !== "Not Calculated" && r.GVI_Score !== undefined);
        const labels = validData.map(r => r.Year);
        const scores = validData.map(r => parseVal(r.GVI_Score));

        trendChartInstance = new Chart(ctxTrend, {
            type: 'line',
            data: {
                labels: labels,
                datasets: [{
                    label: 'GVI Score',
                    data: scores,
                    borderColor: '#e74c3c',
                    backgroundColor: 'rgba(231, 76, 60, 0.2)',
                    fill: true,
                    tension: 0.4
                }]
            },
            options: {
                responsive: true, maintainAspectRatio: false,
                scales: { 
                    y: { min: 0, max: 100, grid: { color: '#333' }, ticks: { color: '#bbb' } },
                    x: { grid: { color: '#333' }, ticks: { color: '#bbb' } }
                },
                plugins: { legend: { display: false } }
            }
        });
    }

    // CSV Download Logic
    const downloadCsvBtn = document.getElementById('downloadCsvBtn');
    if (downloadCsvBtn) {
        downloadCsvBtn.addEventListener('click', () => {
            if (globalData.length === 0) return;
            
            const headers = Object.keys(globalData[0]).filter(k => k !== 'Newick_Tree');
            const csvRows = [];
            csvRows.push(headers.join(','));
            
            for (const row of globalData) {
                const values = headers.map(header => {
                    const val = row[header];
                    return `"${val}"`;
                });
                csvRows.push(values.join(','));
            }
            
            const csvData = csvRows.join('\n');
            const blob = new Blob([csvData], { type: 'text/csv' });
            const url = window.URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.setAttribute('hidden', '');
            a.setAttribute('href', url);
            a.setAttribute('download', 'GVI_Genomic_Indices_Report.csv');
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
        });
    }
});
