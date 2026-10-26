/**
 * PhyloGap Studio — Precision Scientific Web Controller
 * Zero emojis, clean genomics tracks, interactive SVG risk modeling,
 * single-base resolution provenance map, and client-side FASTA/audit export.
 */

let cachedPresets = [];
let lastReconstructionResult = null;

// Tab Navigation
function switchTab(tabId) {
  document.querySelectorAll('.nav-tab').forEach(btn => btn.classList.remove('active'));
  document.querySelectorAll('.tab-pane').forEach(pane => pane.classList.remove('active'));

  const targetPane = document.getElementById(`tab-${tabId}`);
  if (targetPane) {
    targetPane.classList.add('active');
  }

  const activeBtn = document.getElementById(`tab-btn-${tabId}`);
  if (activeBtn) {
    activeBtn.classList.add('active');
  }

  if (tabId === 'calculator') {
    const currentMya = document.getElementById('mya-slider').value || 11.5;
    updateRiskCalc(currentMya);
  }
}

// Preset Loader & Segmented Pill Handlers
async function fetchPresets() {
  try {
    const res = await fetch('/api/presets');
    if (!res.ok) throw new Error(`HTTP error ${res.status}`);
    const data = await res.json();
    cachedPresets = data.presets || [];
    
    // Auto load default preset (Tier 2: Zebra Finch)
    selectPresetByTier(2);
  } catch (err) {
    console.error('Failed to load presets:', err);
  }
}

function selectPresetByTier(tierId) {
  // Update UI pills
  document.querySelectorAll('.preset-pill').forEach(pill => pill.classList.remove('active'));
  const activePill = document.getElementById(`pill-tier-${tierId}`);
  if (activePill) {
    activePill.classList.add('active');
  }

  const selector = document.getElementById('preset-selector');
  if (selector) {
    selector.value = tierId;
    loadSelectedPreset();
  }
}

function loadSelectedPreset() {
  const selector = document.getElementById('preset-selector');
  if (!selector) return;
  const tierId = parseInt(selector.value, 10);
  if (isNaN(tierId)) return;

  const preset = cachedPresets.find(p => p.tier_id === tierId);
  if (!preset) return;

  document.getElementById('flank-left').value = preset.flank_left;
  document.getElementById('gap-len').value = preset.expected_gap_len;
  document.getElementById('div-mya').value = preset.divergence_mya;
  document.getElementById('flank-right').value = preset.flank_right;
  document.getElementById('ref-locus').value = preset.homology_reference;
  document.getElementById('sparse-reads').value = (preset.sparse_reads || []).join('\n');
  document.getElementById('ground-truth').value = preset.gap_truth || '';

  // Update Schematic Display
  const speciesDisplay = document.getElementById('schematic-species-display');
  if (speciesDisplay) {
    speciesDisplay.textContent = `Target: ${preset.target_species} | Homology Ref: ${preset.reference_species} (${preset.divergence_mya} Mya)`;
  }
  const gapLabel = document.getElementById('schematic-gap-label');
  if (gapLabel) {
    gapLabel.textContent = `TARGET GAP [ ${preset.expected_gap_len} bp ]`;
  }

  // Update Field Character Counters
  updateAllFieldCounters();

  // Sync risk calculator slider
  const slider = document.getElementById('mya-slider');
  if (slider) {
    slider.value = preset.divergence_mya;
    updateRiskCalc(preset.divergence_mya);
  }
}

function updateAllFieldCounters() {
  const count = (id) => (document.getElementById(id)?.value || '').trim().length;
  
  const cFlankL = document.getElementById('counter-flank-left');
  if (cFlankL) cFlankL.textContent = `${count('flank-left')} bp`;

  const cFlankR = document.getElementById('counter-flank-right');
  if (cFlankR) cFlankR.textContent = `${count('flank-right')} bp`;

  const cRef = document.getElementById('counter-ref-locus');
  if (cRef) cRef.textContent = `${count('ref-locus')} bp`;

  const cReads = document.getElementById('counter-sparse-reads');
  if (cReads) {
    const lines = (document.getElementById('sparse-reads')?.value || '').split('\n').filter(s => s.trim().length > 0);
    cReads.textContent = `${lines.length} reads`;
  }

  const cGt = document.getElementById('counter-ground-truth');
  if (cGt) cGt.textContent = `${count('ground-truth')} bp`;
}

// Bind live counter listeners
['flank-left', 'flank-right', 'ref-locus', 'sparse-reads', 'ground-truth'].forEach(id => {
  const el = document.getElementById(id);
  if (el) {
    el.addEventListener('input', updateAllFieldCounters);
  }
});

// Reconstruct Form Submission
async function handleReconstruct(event) {
  event.preventDefault();
  const btnRun = document.getElementById('btn-run');
  const origBtnContent = btnRun.innerHTML;
  btnRun.disabled = true;
  btnRun.innerHTML = `
    <svg class="icon-btn animate-spin" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
      <line x1="12" y1="2" x2="12" y2="6"></line>
      <line x1="12" y1="18" x2="12" y2="22"></line>
      <line x1="4.93" y1="4.93" x2="7.76" y2="7.76"></line>
      <line x1="16.24" y1="16.24" x2="19.07" y2="19.07"></line>
      <line x1="2" y1="12" x2="6" y2="12"></line>
      <line x1="18" y1="12" x2="22" y2="12"></line>
      <line x1="4.93" y1="19.07" x2="7.76" y2="16.24"></line>
      <line x1="16.24" y1="7.76" x2="19.07" y2="4.93"></line>
    </svg>
    <span>Synthesizing Multimodal Evidence...</span>
  `;

  const flankL = document.getElementById('flank-left').value.trim();
  const flankR = document.getElementById('flank-right').value.trim();
  const ref = document.getElementById('ref-locus').value.trim();
  const gapLen = parseInt(document.getElementById('gap-len').value, 10);
  const divMya = parseFloat(document.getElementById('div-mya').value);
  const readsText = document.getElementById('sparse-reads').value.trim();
  const sparseReads = readsText ? readsText.split('\n').map(s => s.trim()).filter(s => s.length > 0) : [];
  const gt = document.getElementById('ground-truth').value.trim() || null;

  try {
    const response = await fetch('/api/reconstruct', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        flank_left: flankL,
        flank_right: flankR,
        homology_reference: ref,
        divergence_mya: divMya,
        expected_gap_len: gapLen,
        sparse_reads: sparseReads,
        ground_truth_gap: gt
      })
    });

    if (!response.ok) {
      const errData = await response.json().catch(() => ({ detail: 'Unknown server error' }));
      throw new Error(errData.detail || `Server error: ${response.status}`);
    }

    const data = await response.json();
    lastReconstructionResult = data;
    renderReconstructionResults(data);

  } catch (err) {
    alert(`Reconstruction Error: ${err.message}`);
    console.error(err);
  } finally {
    btnRun.disabled = false;
    btnRun.innerHTML = origBtnContent;
  }
}

function renderReconstructionResults(res) {
  // 1. Verdict Badge
  const verdictContainer = document.getElementById('verdict-container');
  if (res.is_accepted) {
    verdictContainer.innerHTML = '<span class="status-badge badge-accepted">ACCEPTED [ EVIDENCE VERIFIED ]</span>';
  } else {
    verdictContainer.innerHTML = '<span class="status-badge badge-rejected">REJECTED [ GATING ACTIVE ]</span>';
  }

  // 2. Metrics Grid
  document.getElementById('metric-conf').textContent = res.confidence_score !== null ? res.confidence_score.toFixed(3) : '--';
  document.getElementById('metric-whom').textContent = res.homology_weight !== null ? res.homology_weight.toFixed(3) : '--';
  document.getElementById('metric-flank-id').textContent = res.flank_identity_mean !== null ? (res.flank_identity_mean * 100).toFixed(1) + '%' : '--';
  document.getElementById('metric-read-corr').textContent = (res.read_corrections_applied || 0) + ' bp';

  // 3. Coordinate Ruler & Provenance Track
  const ribbon = document.getElementById('provenance-ribbon');
  const ruler = document.getElementById('sequence-ruler');
  ribbon.innerHTML = '';
  ruler.innerHTML = '';

  const provenance = res.provenance || [];
  const totalBases = provenance.length;

  if (totalBases > 0) {
    // Render coordinate ticks (e.g. 1, 11, 21, 31, ...)
    let tickHtml = '';
    for (let i = 0; i < totalBases; i += 10) {
      tickHtml += `<span class="ruler-tick" style="width: ${10 * 22 + 6}px; text-align: left;">${i + 1}</span>`;
    }
    ruler.innerHTML = tickHtml;

    // Render sequence in grouped chunks of 10
    let currentGroup = null;
    provenance.forEach((item, idx) => {
      if (idx % 10 === 0) {
        currentGroup = document.createElement('div');
        currentGroup.className = 'base-block-group';
        ribbon.appendChild(currentGroup);
      }

      const span = document.createElement('span');
      span.className = `base-tag source-${item.source}`;
      span.textContent = item.base;
      span.title = `Pos ${idx + 1}: ${item.base} | Source: ${item.source} (Original: ${item.original})`;
      currentGroup.appendChild(span);
    });
  } else {
    ribbon.innerHTML = `<span class="track-placeholder">${res.final_sequence}</span>`;
  }

  // 4. Selective Gating Callout
  const rationaleBox = document.getElementById('gating-rationale-box');
  const rationaleText = document.getElementById('gating-rationale-text');
  if (!res.is_accepted) {
    rationaleBox.style.display = 'block';
    rationaleText.textContent = res.rejection_reason || 'Composite confidence score fell below dynamic threshold C_min(t). Unresolved gap preserved as masked Ns to prevent assembly false-closures.';
  } else {
    rationaleBox.style.display = 'none';
  }

  // 5. Ground Truth Evaluation Callout
  const gtBox = document.getElementById('gt-eval-box');
  const gtMetrics = document.getElementById('gt-eval-metrics');
  const gtEval = res.ground_truth_eval || res.evaluation;
  
  if (gtEval) {
    gtBox.style.display = 'block';
    const accPct = (gtEval.nucleotide_accuracy * 100).toFixed(1);
    const emStatus = gtEval.exact_match ? '<span class="text-success font-semibold">YES (100% IDENTITY)</span>' : '<span class="text-warning font-semibold">SUBSTITUTION DRIFT</span>';
    const fcStatus = gtEval.is_false_closure ? '<span class="text-danger font-semibold">FALSE CLOSURE DETECTED</span>' : '<span class="text-success font-semibold">SAFE RECONSTRUCTION</span>';

    gtMetrics.innerHTML = `
      <div style="display: flex; gap: 1.5rem; flex-wrap: wrap; margin-bottom: 0.5rem;">
        <div>Accuracy: <strong class="text-primary">${accPct}%</strong></div>
        <div>Exact Match: ${emStatus}</div>
        <div>Safety Status: ${fcStatus}</div>
      </div>
      <div style="font-family: var(--font-mono); font-size: 0.72rem; color: var(--text-dim);">
        Levenshtein Distance: ${gtEval.edit_distance} bp | Substitutions: ${gtEval.substitutions} bp | Indels: ${gtEval.indels} bp
      </div>
    `;
  } else {
    gtBox.style.display = 'none';
  }
}

// Sequence Clipboard Copy
function copyReconstructedSequence() {
  if (!lastReconstructionResult || !lastReconstructionResult.final_sequence) {
    alert('Please execute reconstruction first.');
    return;
  }
  navigator.clipboard.writeText(lastReconstructionResult.final_sequence).then(() => {
    const copyBtnText = document.getElementById('btn-copy-seq-text');
    if (copyBtnText) {
      const orig = copyBtnText.textContent;
      copyBtnText.textContent = 'Copied!';
      setTimeout(() => { copyBtnText.textContent = orig; }, 2000);
    }
  });
}

// Client-side FASTA File Export
function exportFasta() {
  if (!lastReconstructionResult || !lastReconstructionResult.final_sequence) {
    alert('Please execute reconstruction first.');
    return;
  }
  const res = lastReconstructionResult;
  const fastaHeader = `>PhyloGap_Reconstructed_Gap|len=${res.final_sequence.length}|status=${res.status}|conf=${res.confidence_score !== null ? res.confidence_score.toFixed(3) : 'NA'}`;
  const fastaContent = `${fastaHeader}\n${res.final_sequence}\n`;

  const blob = new Blob([fastaContent], { type: 'text/plain;charset=utf-8' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = `phylogap_reconstruction_${Date.now()}.fa`;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

// Tab 2: Evolutionary Risk Calculator
function setSliderMya(val) {
  const slider = document.getElementById('mya-slider');
  if (slider) {
    slider.value = val;
    updateRiskCalc(val);
  }
}

async function updateRiskCalc(mya) {
  const t = Math.max(0, parseFloat(mya));
  const sliderDisplay = document.getElementById('slider-val-display');
  if (sliderDisplay) sliderDisplay.textContent = `${t.toFixed(1)} Mya`;

  try {
    const res = await fetch(`/api/risk-calc?mya=${t}`);
    if (!res.ok) throw new Error(`HTTP error ${res.status}`);
    const data = await res.json();

    // Update risk UI
    const riskBadge = document.getElementById('risk-badge');
    const riskExpectedId = document.getElementById('risk-expected-id');
    const riskGuidance = document.getElementById('risk-guidance-text');

    if (riskBadge) {
      riskBadge.textContent = data.risk_label;
      riskBadge.style.color = data.risk_color;
    }
    if (riskExpectedId) {
      riskExpectedId.textContent = `${data.expected_nucleotide_identity.toFixed(1)}%`;
      riskExpectedId.style.color = data.risk_color;
    }
    if (riskGuidance) {
      riskGuidance.innerHTML = `
        <div style="margin-bottom: 0.4rem;">${data.guidance}</div>
        <div style="font-size: 0.75rem; color: var(--text-dim); font-family: var(--font-mono);">
          Dynamic Weight: <code style="color: var(--primary);">w_hom(${t}) = ${data.homology_weight_w_hom}</code>
          ${data.is_exceeding_t_crit ? ' — <strong style="color: #F43F5E;">CRITICAL THRESHOLD EXCEEDED (w_hom = 0)</strong>' : ''}
        </div>
      `;
    }

    renderSvgDecayChart(t, data.expected_nucleotide_identity, data.critical_threshold_t_crit);

  } catch (err) {
    console.error('Failed to update risk calculator:', err);
  }
}

function renderSvgDecayChart(currentMya, currentIdentity, tCrit) {
  const svg = document.getElementById('decay-chart');
  if (!svg) return;

  const w = 760;
  const h = 240;
  const padLeft = 60;
  const padRight = 30;
  const padTop = 30;
  const padBottom = 40;

  const plotW = w - padLeft - padRight;
  const plotH = h - padTop - padBottom;

  const minMya = 0;
  const maxMya = 90;
  const minAcc = 50;
  const maxAcc = 100;

  const scaleX = (m) => padLeft + ((m - minMya) / (maxMya - minMya)) * plotW;
  const scaleY = (acc) => padTop + plotH - ((acc - minAcc) / (maxAcc - minAcc)) * plotH;

  let svgContent = '';

  // Defs for gradients
  svgContent += `
    <defs>
      <linearGradient id="curveGradient" x1="0" y1="0" x2="0" y2="1">
        <stop offset="0%" stop-color="#38BDF8" stop-opacity="0.25" />
        <stop offset="100%" stop-color="#38BDF8" stop-opacity="0.0" />
      </linearGradient>
      <pattern id="dangerHatch" width="8" height="8" patternTransform="rotate(45 0 0)" patternUnits="userSpaceOnUse">
        <line x1="0" y1="0" x2="0" y2="8" stroke="rgba(244, 63, 94, 0.25)" stroke-width="1.5" />
      </pattern>
    </defs>
  `;

  // Grid Lines & Labels
  for (let acc = 60; acc <= 100; acc += 10) {
    const y = scaleY(acc);
    svgContent += `<line x1="${padLeft}" y1="${y}" x2="${w - padRight}" y2="${y}" stroke="rgba(255,255,255,0.06)" stroke-dasharray="3,3" />`;
    svgContent += `<text x="${padLeft - 10}" y="${y + 4}" fill="#64748B" font-size="10" text-anchor="end" font-family="JetBrains Mono">${acc}%</text>`;
  }

  for (let m = 0; m <= 90; m += 15) {
    const x = scaleX(m);
    svgContent += `<line x1="${x}" y1="${padTop}" x2="${x}" y2="${h - padBottom}" stroke="rgba(255,255,255,0.06)" stroke-dasharray="3,3" />`;
    svgContent += `<text x="${x}" y="${h - padBottom + 18}" fill="#64748B" font-size="10" text-anchor="middle" font-family="JetBrains Mono">${m} Mya</text>`;
  }

  // Red Danger Zone Shading (t >= tCrit)
  const xCrit = scaleX(tCrit);
  const dangerWidth = (w - padRight) - xCrit;
  if (dangerWidth > 0) {
    svgContent += `
      <rect x="${xCrit}" y="${padTop}" width="${dangerWidth}" height="${plotH}" fill="url(#dangerHatch)" />
      <text x="${xCrit + 12}" y="${padTop + 24}" fill="#F43F5E" font-size="10" font-weight="600" font-family="Inter" letter-spacing="0.05em">DANGER ZONE (FCR &gt; 90%)</text>
    `;
  }

  // Vertical Line at tCrit
  svgContent += `
    <line x1="${xCrit}" y1="${padTop}" x2="${xCrit}" y2="${h - padBottom}" stroke="#F43F5E" stroke-width="1.5" stroke-dasharray="4,4" />
    <text x="${xCrit}" y="${padTop - 10}" fill="#F43F5E" font-size="10" font-weight="700" text-anchor="middle" font-family="JetBrains Mono">t_crit = ${tCrit} Mya</text>
  `;

  // Curve and Area under Curve
  const lambda = 0.006395;
  let pathD = '';
  let areaD = `M ${scaleX(0)} ${scaleY(minAcc)}`;
  
  for (let m = 0; m <= 90; m += 1) {
    const acc = Math.exp(-lambda * m) * 100.0;
    const x = scaleX(m);
    const y = scaleY(acc);
    if (m === 0) {
      pathD += `M ${x} ${y}`;
      areaD += ` L ${x} ${y}`;
    } else {
      pathD += ` L ${x} ${y}`;
      areaD += ` L ${x} ${y}`;
    }
  }
  areaD += ` L ${scaleX(90)} ${scaleY(minAcc)} Z`;

  svgContent += `<path d="${areaD}" fill="url(#curveGradient)" />`;
  svgContent += `<path d="${pathD}" fill="none" stroke="#38BDF8" stroke-width="2.5" />`;

  // Current Position Cursor
  const curX = scaleX(currentMya);
  const curY = scaleY(currentIdentity);
  svgContent += `
    <line x1="${curX}" y1="${padTop}" x2="${curX}" y2="${h - padBottom}" stroke="#F59E0B" stroke-width="1.2" stroke-dasharray="2,2" />
    <circle cx="${curX}" cy="${curY}" r="5" fill="#F59E0B" stroke="#080C14" stroke-width="2" />
    <circle cx="${curX}" cy="${curY}" r="10" fill="none" stroke="rgba(245, 158, 11, 0.4)" stroke-width="1.5">
      <animate attributeName="r" values="6;12;6" dur="2s" repeatCount="indefinite"/>
      <animate attributeName="opacity" values="0.8;0.2;0.8" dur="2s" repeatCount="indefinite"/>
    </circle>
    <text x="${curX + (currentMya > 70 ? -10 : 10)}" y="${curY - 10}" fill="#F59E0B" font-size="11" font-weight="700" font-family="JetBrains Mono" text-anchor="${currentMya > 70 ? 'end' : 'start'}">
      ${currentIdentity.toFixed(1)}% (${currentMya.toFixed(1)} Mya)
    </text>
  `;

  svg.innerHTML = svgContent;
}

// Initialize on DOM ready
document.addEventListener('DOMContentLoaded', () => {
  fetchPresets();
  updateRiskCalc(11.5);
});
