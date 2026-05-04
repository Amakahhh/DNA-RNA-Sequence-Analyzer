/**
 * TideLab — Frontend Logic
 * Ocean-themed variable and function names throughout.
 * All API communication, result rendering, and UI state management lives here.
 */

'use strict';

// ── DOM Anchors ────────────────────────────────────────────────────────────────

const sequenceHarbor   = document.getElementById('sequenceInput');
const fileCurrentInput = document.getElementById('fileInput');
const charTideCounter  = document.getElementById('charCounter');
const strandCompass    = document.getElementById('strandCompass');
const launchButton     = document.getElementById('launchBtn');
const clearButton      = document.getElementById('clearBtn');
const dropHarbor       = document.getElementById('dropHarbor');
const errorBanner      = document.getElementById('errorBanner');
const resultsBay       = document.getElementById('resultsBay');
const analysisLoader   = document.getElementById('analysisLoader');

// ── Sample sequences (one DNA with stop codon, one RNA fragment) ────────────────

const TIDE_SAMPLES = {
  dna_coding: 'ATGAAAGCAATTTTCGTACTGAAAGGTTTTGTTGGTTTTCTTAAAGCATGGCCCAGCCCTGCGTTTGAGCAAATGATCCTTTAG',
  dna_template: 'CTAAAGGATCATTTGCTCAAACGCAGGGCTGGGCCATGCTTTAAGAAAACCAACAAAACCTTTCAGTACGAAAATTGCTTTCAT',
  rna: 'AUGAAAGCAAUUUUCGUACUGAAAGGUUUUGUUGGUUUUCUUAAAGCAUGGCCCAGCCCUGCGUUUGAGCAAAUGAUCCUUUAG',
};

// ── Character counter ─────────────────────────────────────────────────────────

function refreshTideCounter() {
  const charCount = sequenceHarbor.value.replace(/\s/g, '').length;
  charTideCounter.textContent = `${charCount.toLocaleString()} bases`;
}

sequenceHarbor.addEventListener('input', refreshTideCounter);

// ── Sample sequence injection ──────────────────────────────────────────────────

document.querySelectorAll('.btn-sample').forEach(sampleBtn => {
  sampleBtn.addEventListener('click', () => {
    const sampleKey = sampleBtn.dataset.sample;
    sequenceHarbor.value = TIDE_SAMPLES[sampleKey] || '';
    refreshTideCounter();

    // Pre-select the matching strand type if relevant
    if (sampleKey === 'dna_template') {
      document.getElementById('strandTemplate').checked = true;
    } else {
      document.getElementById('strandCoding').checked = true;
    }
  });
});

// ── File upload ───────────────────────────────────────────────────────────────

document.getElementById('uploadBtn').addEventListener('click', () => {
  fileCurrentInput.click();
});

fileCurrentInput.addEventListener('change', (evt) => {
  const beachFile = evt.target.files[0];
  if (beachFile) castNetFile(beachFile);
  fileCurrentInput.value = ''; // reset so the same file can be re-loaded
});

function castNetFile(beachFile) {
  const allowedTypes = ['.txt', '.fasta', '.fa', '.seq'];
  const fileName     = beachFile.name.toLowerCase();
  const isAllowed    = allowedTypes.some(ext => fileName.endsWith(ext));

  // Also allow files with no extension (common for FASTA)
  if (!isAllowed && fileName.includes('.')) {
    showErrorBanner(`Unsupported file type. Please upload a .txt or .fasta file.`);
    return;
  }

  const tideCurrent = new FileReader();
  tideCurrent.onload = (e) => {
    sequenceHarbor.value = e.target.result;
    refreshTideCounter();
    hideErrorBanner();
  };
  tideCurrent.onerror = () => {
    showErrorBanner('Could not read the file. Please try again.');
  };
  tideCurrent.readAsText(beachFile);
}

// ── Drag and drop ─────────────────────────────────────────────────────────────

dropHarbor.addEventListener('dragenter', (e) => { e.preventDefault(); dropHarbor.classList.add('drag-over'); });
dropHarbor.addEventListener('dragover',  (e) => { e.preventDefault(); });
dropHarbor.addEventListener('dragleave', (e) => {
  // Only remove class if leaving the drop zone itself (not a child element)
  if (!dropHarbor.contains(e.relatedTarget)) {
    dropHarbor.classList.remove('drag-over');
  }
});

dropHarbor.addEventListener('drop', (e) => {
  e.preventDefault();
  dropHarbor.classList.remove('drag-over');
  const tideSurge = e.dataTransfer.files[0];
  if (tideSurge) {
    castNetFile(tideSurge);
  } else {
    // Handle dragged plain text
    const rawText = e.dataTransfer.getData('text');
    if (rawText) {
      sequenceHarbor.value = rawText;
      refreshTideCounter();
    }
  }
});

// ── Clear ─────────────────────────────────────────────────────────────────────

clearButton.addEventListener('click', () => {
  sequenceHarbor.value = '';
  refreshTideCounter();
  hideErrorBanner();
  resultsBay.classList.remove('visible');
  strandCompass.classList.remove('visible');
  analysisLoader.classList.remove('visible');
  clearBlastPollTimer();
});

// ── Launch analysis ───────────────────────────────────────────────────────────

launchButton.addEventListener('click', harborAnalyse);
sequenceHarbor.addEventListener('keydown', (e) => {
  if (e.ctrlKey && e.key === 'Enter') harborAnalyse();
});

async function harborAnalyse() {
  const rawTide = sequenceHarbor.value.trim();
  if (!rawTide) {
    showErrorBanner('Please enter a nucleotide sequence before launching the analysis.');
    return;
  }

  const selectedStrand = document.querySelector('input[name="strandPolarity"]:checked');
  const strandPolarity = selectedStrand ? selectedStrand.value : 'coding';

  hideErrorBanner();
  launchButton.disabled = true;
  resultsBay.classList.remove('visible');
  analysisLoader.classList.add('visible');
  clearBlastPollTimer();

  try {
    const tidewaterResponse = await fetch('/api/analyse', {
      method:  'POST',
      headers: { 'Content-Type': 'application/json' },
      body:    JSON.stringify({ sequence: rawTide, strand_type: strandPolarity }),
    });

    const reefData = await tidewaterResponse.json();

    if (!tidewaterResponse.ok) {
      throw new Error(reefData.error || 'Analysis failed — please check your sequence.');
    }

    analysisLoader.classList.remove('visible');
    renderAllTidalSteps(reefData);

    // Show the DNA strand selector based on what was detected
    if (reefData.detection?.ocean_layer === 'DNA') {
      strandCompass.classList.add('visible');
    }

  } catch (reefFault) {
    analysisLoader.classList.remove('visible');
    showErrorBanner(reefFault.message || 'An unexpected error occurred. Please try again.');
  } finally {
    launchButton.disabled = false;
  }
}

// ── Render all result cards ────────────────────────────────────────────────────

function renderAllTidalSteps(reefData) {
  resultsBay.classList.add('visible');

  paintDetectionBay(reefData.detection);
  paintTranscriptionCurrent(reefData.transcription, reefData.detection);
  paintCodonWaveLog(reefData.translation);
  paintPearlNecklace(reefData.translation);
  paintReefProtein(reefData.protein, reefData.blast_job);

  // Stagger the card float-in animations
  const tidalCards = document.querySelectorAll('.tidal-card');
  tidalCards.forEach((card, waveIndex) => {
    card.classList.remove('floated');
    setTimeout(() => card.classList.add('floated'), waveIndex * 120);
  });

  // Scroll to results
  setTimeout(() => {
    resultsBay.scrollIntoView({ behavior: 'smooth', block: 'start' });
  }, 200);
}

// ── Step 1: Detection ─────────────────────────────────────────────────────────

function paintDetectionBay(detection) {
  const oceanLayer  = detection.ocean_layer;
  const badgeClass  = oceanLayer === 'DNA' ? 'dna-badge' : 'rna-badge';
  const badgeIcon   = oceanLayer === 'DNA' ? '🧬' : '🌀';

  // Base tally chips
  const baseOrder = oceanLayer === 'DNA' ? ['A', 'T', 'C', 'G', 'N'] : ['A', 'U', 'C', 'G', 'N'];
  const baseChips = baseOrder
    .filter(base => detection.base_tally[base] !== undefined)
    .map(base => {
      const count = detection.base_tally[base];
      const pct   = ((count / detection.length) * 100).toFixed(1);
      return `<span class="base-chip base-${base}">${base}: ${count.toLocaleString()} (${pct}%)</span>`;
    })
    .join('');

  document.getElementById('detectionContent').innerHTML = `
    <div class="detection-badge ${badgeClass}">
      ${badgeIcon} <span>${oceanLayer}</span>
    </div>
    <div class="sequence-stats-grid">
      <div class="stat-chip">
        <span class="stat-value">${detection.length.toLocaleString()}</span>
        <span class="stat-label">Total Bases</span>
      </div>
      <div class="stat-chip">
        <span class="stat-value">${detection.gc_content}%</span>
        <span class="stat-label">GC Content</span>
      </div>
      <div class="stat-chip">
        <span class="stat-value">${Object.keys(detection.base_tally).length}</span>
        <span class="stat-label">Unique Bases</span>
      </div>
    </div>
    <div class="base-tally">${baseChips}</div>
    <div class="explanation-panel">${detection.explanation}</div>
  `;
}

// ── Step 2: Transcription ─────────────────────────────────────────────────────

function paintTranscriptionCurrent(transcription, detection) {
  const inputSeq  = detection.sequence;
  const mrnaSeq   = transcription.surface_ripple;
  const isRNA     = detection.ocean_layer === 'RNA';

  const inputLabel  = isRNA ? 'RNA Input' : 'DNA Input';
  const inputClass  = isRNA ? 'mrna-strand' : 'dna-strand';
  const arrowHtml   = isRNA ? '' : `<div class="transcription-arrow">↓ Transcription</div>`;
  const mrnaSection = isRNA ? '' : `
    ${arrowHtml}
    <div class="seq-row">
      <span class="seq-label">mRNA</span>
      <div class="seq-display mrna-strand">${escapeHtml(mrnaSeq)}</div>
    </div>
  `;

  document.getElementById('transcriptionContent').innerHTML = `
    <div class="sequence-comparison">
      <div class="seq-row">
        <span class="seq-label">${inputLabel}</span>
        <div class="seq-display ${inputClass}">${escapeHtml(inputSeq)}</div>
      </div>
      ${mrnaSection}
    </div>
    <div class="explanation-panel">${transcription.explanation}</div>
  `;
}

// ── Step 3: Translation — Codon table ─────────────────────────────────────────

function paintCodonWaveLog(translation) {
  const { codon_log, explanation } = translation;

  const rowsHtml = codon_log.map(entry => {
    const statusClass  = `row-${entry.status}`;
    const badgeClass   = `badge-${entry.status}`;
    const badgeText    = statusBadgeLabel(entry.status);

    const singleHtml = entry.status === 'coding' || entry.status === 'start'
      ? `<span class="single-letter-aa">${entry.single_letter}</span>`
      : `<span style="color:var(--text-muted)">${entry.single_letter}</span>`;

    return `
      <tr class="${statusClass}">
        <td style="color:var(--text-muted);font-variant-numeric:tabular-nums">#${entry.position}</td>
        <td><span class="codon-letters">${entry.codon}</span></td>
        <td>${singleHtml}</td>
        <td>${entry.three_letter}</td>
        <td>${entry.amino_acid}</td>
        <td><span class="codon-badge ${badgeClass}">${badgeText}</span></td>
      </tr>
    `;
  }).join('');

  document.getElementById('translationContent').innerHTML = `
    <div class="codon-table-wrap">
      <table class="codon-table">
        <thead>
          <tr>
            <th>#</th>
            <th>Codon</th>
            <th>1-Letter</th>
            <th>3-Letter</th>
            <th>Amino Acid</th>
            <th>Type</th>
          </tr>
        </thead>
        <tbody>${rowsHtml}</tbody>
      </table>
    </div>
    <div class="explanation-panel">${explanation}</div>
  `;
}

function statusBadgeLabel(status) {
  const waveLabels = {
    start:      'Start',
    stop:       'Stop',
    coding:     'Coding',
    utr:        "5' UTR",
    incomplete: 'Incomplete',
    unknown:    'Unknown',
  };
  return waveLabels[status] || status;
}

// ── Step 4: Amino Acids — Pearl necklace ──────────────────────────────────────

function paintPearlNecklace(translation) {
  const { pearl_necklace, explanation } = translation;

  if (!pearl_necklace || pearl_necklace.length === 0) {
    document.getElementById('aminoAcidsContent').innerHTML = `
      <p style="color:var(--text-muted);font-size:0.88rem">
        No amino acids were produced from this sequence.
      </p>
      <div class="explanation-panel">${explanation}</div>
    `;
    return;
  }

  // Bead chain
  const beadsHtml = pearl_necklace.map((bead, waveIndex) => `
    <div class="pearl-bead" title="${bead.name} (${bead.three} / ${bead.one})">
      <div class="bead-circle">${bead.one}</div>
      <span class="bead-index">${waveIndex + 1}</span>
    </div>
  `).join('');

  // Unique amino acids table (sorted by frequency)
  const tallied = {};
  pearl_necklace.forEach(bead => {
    if (!tallied[bead.one]) tallied[bead.one] = { ...bead, count: 0 };
    tallied[bead.one].count++;
  });
  const sortedEntries = Object.values(tallied).sort((a, b) => b.count - a.count);

  const tableRowsHtml = sortedEntries.map((entry, idx) => `
    <tr>
      <td style="color:var(--text-muted)">${idx + 1}</td>
      <td><span class="aa-one-letter">${entry.one}</span></td>
      <td style="font-weight:600">${entry.name}</td>
      <td style="font-family:monospace;color:var(--text-secondary)">${entry.three}</td>
      <td style="color:var(--wave-blue)">${entry.count}</td>
    </tr>
  `).join('');

  document.getElementById('aminoAcidsContent').innerHTML = `
    <div class="pearl-chain">${beadsHtml}</div>
    <div class="amino-acid-table-wrap">
      <table class="amino-acid-table">
        <thead>
          <tr>
            <th>#</th>
            <th>Code</th>
            <th>Full Name</th>
            <th>3-Letter</th>
            <th>Count</th>
          </tr>
        </thead>
        <tbody>${tableRowsHtml}</tbody>
      </table>
    </div>
    <div class="explanation-panel">${explanation}</div>
  `;
}

// ── Step 5: Protein characterisation ─────────────────────────────────────────

function paintReefProtein(protein, blastJob) {
  if (!protein || !protein.valid) {
    document.getElementById('proteinContent').innerHTML = `
      <div class="blast-error">
        ${protein?.reason || 'No protein could be characterised from this sequence.'}
      </div>
    `;
    return;
  }

  // Composition bar chart (top 8)
  const topComposition = protein.composition_table.slice(0, 8);
  const maxCompCount   = topComposition[0]?.count || 1;
  const compBarsHtml   = topComposition.map(entry => `
    <div class="comp-row">
      <span class="comp-aa-name">${entry.name} <span style="color:var(--text-muted)">(${entry.three})</span></span>
      <div class="comp-bar-track">
        <div class="comp-bar-fill" style="width:${(entry.count / maxCompCount * 100).toFixed(1)}%"></div>
      </div>
      <span class="comp-pct">${entry.pct}%</span>
    </div>
  `).join('');

  document.getElementById('proteinContent').innerHTML = `
    <div class="protein-sequence-box">${escapeHtml(protein.clean_sequence)}</div>

    <div class="property-grid">
      <div class="property-chip">
        <span class="prop-value">${protein.total_residues}</span>
        <span class="prop-label">Amino Acids</span>
      </div>
      <div class="property-chip">
        <span class="prop-value">${(protein.estimated_mass / 1000).toFixed(2)} kDa</span>
        <span class="prop-label">Est. Mass</span>
      </div>
      <div class="property-chip">
        <span class="prop-value">${protein.hydrophobic_pct}%</span>
        <span class="prop-label">Hydrophobic</span>
      </div>
      <div class="property-chip">
        <span class="prop-value" style="font-size:0.85rem">${protein.charge_nature.split(' ')[0]}</span>
        <span class="prop-label">Charge</span>
      </div>
    </div>

    <div class="composition-chart">
      <h4>Amino Acid Composition (top residues)</h4>
      ${compBarsHtml}
    </div>

    <div class="explanation-panel">${protein.explanation}</div>

    <div class="blast-section" id="blastSection">
      <h3>🔍 Database Search — UniProtKB/Swiss-Prot</h3>
      <div id="blastResultsZone"></div>
    </div>
  `;

  // Kick off BLAST result display
  anchorBlastDisplay(blastJob);
}

// ── BLAST polling ─────────────────────────────────────────────────────────────

let blastPollTimer = null;

function clearBlastPollTimer() {
  if (blastPollTimer !== null) {
    clearTimeout(blastPollTimer);
    blastPollTimer = null;
  }
}

function anchorBlastDisplay(blastJob) {
  const blastZone = document.getElementById('blastResultsZone');
  if (!blastZone) return;

  if (!blastJob || blastJob.error) {
    blastZone.innerHTML = `
      <div class="blast-error">${blastJob?.error || 'Database search not available for this sequence.'}</div>
    `;
    return;
  }

  if (!blastJob.job_id) {
    blastZone.innerHTML = `<div class="blast-error">BLAST job ID was not returned by the server.</div>`;
    return;
  }

  // Show spinner while polling
  blastZone.innerHTML = `
    <div class="blast-status-bar">
      <div class="blast-spinner"></div>
      <span>Submitting sequence to UniProtKB/Swiss-Prot via BLAST… this may take 30–90 seconds.</span>
    </div>
  `;

  // Start polling with increasing interval (3s → 5s → 8s → 10s → …)
  const pollSchedule = [3000, 5000, 8000, 10000, 10000, 10000, 10000, 10000, 10000, 10000];
  let pollDepth = 0;

  function scheduleNextTide() {
    if (pollDepth >= pollSchedule.length) {
      blastZone.innerHTML = `
        <div class="blast-error">
          BLAST is taking longer than expected. The EBI server may be busy.
          Try re-analysing your sequence in a few minutes.
        </div>
      `;
      return;
    }
    blastPollTimer = setTimeout(() => pollBlastTide(blastJob.job_id, scheduleNextTide), pollSchedule[pollDepth++]);
  }

  scheduleNextTide();
}

async function pollBlastTide(jobId, onPending) {
  const blastZone = document.getElementById('blastResultsZone');
  if (!blastZone) return;

  try {
    const waveResponse = await fetch(`/api/blast_result/${encodeURIComponent(jobId)}`);
    const reefResult   = await waveResponse.json();

    if (reefResult.status === 'pending') {
      blastZone.innerHTML = `
        <div class="blast-status-bar">
          <div class="blast-spinner"></div>
          <span>${reefResult.message}</span>
        </div>
      `;
      onPending();
      return;
    }

    if (reefResult.status === 'error') {
      blastZone.innerHTML = `<div class="blast-error">${reefResult.message}</div>`;
      return;
    }

    if (reefResult.status === 'complete') {
      renderBlastMatches(reefResult, blastZone);
    }

  } catch (networkFault) {
    blastZone.innerHTML = `
      <div class="blast-error">Network error while checking BLAST status: ${networkFault.message}</div>
    `;
  }
}

function renderBlastMatches(reefResult, blastZone) {
  if (!reefResult.matches || reefResult.matches.length === 0) {
    blastZone.innerHTML = `
      <div class="blast-status-bar" style="border-color:rgba(255,169,77,0.3)">
        ⚠️ No significant matches found in Swiss-Prot.
        The sequence may be novel, very short, or highly mutated relative to known proteins.
      </div>
    `;
    return;
  }

  const hitsHtml = reefResult.matches.map(hit => {
    const identStr  = typeof hit.identity === 'number' ? `${hit.identity.toFixed(1)}%` : String(hit.identity);
    const evalStr   = formatEValue(hit.evalue);

    return `
      <div class="blast-hit-card">
        <div class="hit-header">
          <span class="hit-protein-name">${escapeHtml(hit.protein_name)}</span>
          <span class="hit-accession">${escapeHtml(hit.accession)}</span>
        </div>
        <div class="hit-organism">📍 ${escapeHtml(hit.organism)}</div>
        <div class="hit-metrics">
          <div class="hit-metric">
            <span class="metric-val">${hit.score}</span>
            <span class="metric-key">Score</span>
          </div>
          <div class="hit-metric">
            <span class="metric-val">${identStr}</span>
            <span class="metric-key">Identity</span>
          </div>
          <div class="hit-metric">
            <span class="metric-val">${evalStr}</span>
            <span class="metric-key">E-value</span>
          </div>
        </div>
      </div>
    `;
  }).join('');

  blastZone.innerHTML = `
    ${hitsHtml}
    <div class="explanation-panel mt-12">${reefResult.explanation}</div>
  `;
}

// ── Utilities ─────────────────────────────────────────────────────────────────

function escapeHtml(str) {
  if (typeof str !== 'string') return String(str);
  return str
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;');
}

function formatEValue(evalue) {
  if (typeof evalue === 'number') {
    if (evalue < 1e-99) return '< 1e-99';
    if (evalue < 0.001) return evalue.toExponential(2);
    return evalue.toFixed(4);
  }
  return String(evalue);
}

function showErrorBanner(message) {
  errorBanner.innerHTML = `⚠️ ${message}`;
  errorBanner.classList.add('visible');
}

function hideErrorBanner() {
  errorBanner.classList.remove('visible');
}
