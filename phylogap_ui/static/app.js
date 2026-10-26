/**
 * PhyloGap Studio — Interactive Frontend Controller
 * Connects UI elements to FastAPI backend endpoints:
 * /api/presets, /api/reconstruct, /api/risk-calc
 */

let cachedPresets = [];

// Tab Navigation
function switchTab(tabId) {
  document.querySelectorAll('.nav-btn').forEach(btn => btn.classList.remove('active'));
  document.querySelectorAll('.tab-pane').forEach(pane => pane.classList.remove('active'));

  const targetPane = document.getElementById(`tab-${tabId}`);
  if (targetPane) {
    targetPane.classList.add('active');
  }

  // Find active button
  const buttons = document.querySelectorAll('.nav-btn');
  if (tabId === 'studio' && buttons[0]) buttons[0].classList.add('active');
  if (tabId === 'calculator') {
    if (buttons[1]) buttons[1].classList.add('active');
    // Ensure SVG chart is rendered cleanly
    const currentMya = document.getElementById('mya-slider').value || 11.5;
    updateRiskCalc(currentMya);
  }
  if (tabId === 'benchmarks' && buttons[2]) buttons[2].classList.add('active');
}

// Preset Loader
async function fetchPresets() {
  try {
    const res = await fetch('/api/presets');
    if (!res.ok) throw new Error(`HTTP error ${res.status}`);
    const data = await res.json();
    cachedPresets = data.presets || [];
    
    // Auto load default preset (Tier 2: Zebra finch)
    const selector = document.getElementById('preset-selector');
    if (selector && selector.value) {
      loadSelectedPreset();
    }
  } catch (err) {
    console.error('Failed to load presets:', err);
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

  // Also update risk calculator slider to match divergence time
  const slider = document.getElementById('mya-slider');
  if (slider) {
    slider.value = preset.divergence_mya;
    updateRiskCalc(preset.divergence_mya);
  }
}

// Reconstruct Form Submission
async function handleReconstruct(event) {
  event.preventDefault();
  const btnRun = document.getElementById('btn-run');
  const origBtnContent = btnRun.innerHTML;
  btnRun.disabled = true;
  btnRun.innerHTML = '<span>⏳ Synthesizing Multimodal Evidence...</span>';

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
    verdictContainer.innerHTML = '<span class="badge-verdict badge-accepted">✅ ACCEPTED (EVIDENCE VERIFIED)</span>';
  } else {
    verdictContainer.innerHTML = '<span class="badge-verdict badge-rejected">🛡️ REJECTED (PROTECTED BY GATING)</span>';
  }

  // 2. Metrics Row
  document.getElementById('metric-conf').textContent = res.confidence_score !== null ? res.confidence_score.toFixed(3) : '--';
  document.getElementById('metric-whom').textContent = res.homology_weight !== null ? res.homology_weight.toFixed(3) : '--';
  document.getElementById('metric-flank-id').textContent = res.flank_identity_mean !== null ? (res.flank_identity_mean * 100).toFixed(1) + '%' : '--';
  document.getElementById('metric-read-corr').textContent = (res.read_corrections_applied || 0) + ' bp';

  // 3. Provenance Ribbon
  const ribbon = document.getElementById('provenance-ribbon');
  ribbon.innerHTML = '';
  if (res.provenance && res.provenance.length > 0) {
    res.provenance.forEach((item, idx) => {
      const span = document.createElement('span');
      span.className = `base-tag source-${item.source}`;
      span.textContent = item.base;
      span.title = `Pos ${idx + 1}: ${item.base} | Source: ${item.source} (Orig: ${item.original})`;
      ribbon.appendChild(span);
    });
  } else {
    ribbon.innerHTML = `<span style="color: var(--text-muted); font-family: var(--font-mono);">${res.final_sequence}</span>`;
  }

  // 4. Gating Rationale Box
  const rationaleBox = document.getElementById('gating-rationale-box');
  const rationaleText = document.getElementById('gating-rationale-text');
  if (!res.is_accepted) {
    rationaleBox.style.display = 'block';
    rationaleText.textContent = res.rejection_reason || 'Confidence score fell below dynamic threshold C_min(t). Gap protected to prevent False Closure.';
  } else {
    rationaleBox.style.display = 'none';
  }

  // 5. Ground Truth Audit Box
  const gtBox = document.getElementById('gt-eval-box');
  const gtMetrics = document.getElementById('gt-eval-metrics');
  if (res.ground_truth_eval) {
    const gt = res.ground_truth_eval;
    gtBox.style.display = 'block';
    const accPct = (gt.nucleotide_accuracy * 100).toFixed(1);
    const emBadge = gt.exact_match ? '<span style="color: #34D399; font-weight: 700;">✅ EXACT MATCH (100%)</span>' : '<span style="color: #FBBF24;">SUBSTITUTION DRIFT</span>';
    const fcBadge = gt.is_false_closure ? '<span style="color: #EF4444; font-weight: 700;">⚠️ FALSE CLOSURE</span>' : '<span style="color: #34D399; font-weight: 700;">🛡️ SAFE RECONSTRUCTION</span>';

    gtMetrics.innerHTML = `
      <div style="display: flex; gap: 1rem; flex-wrap: wrap; margin-bottom: 0.4rem;">
        <div>Acc_nuc: <strong style="color: #38BDF8;">${accPct}%</strong></div>
        <div>Status: ${emBadge}</div>
        <div>Safety: ${fcBadge}</div>
      </div>
      <div style="font-size: 0.8rem; color: var(--text-muted); font-family: var(--font-mono);">
        Levenshtein Distance: ${gt.edit_distance} bp | Substitutions: ${gt.substitutions} bp | Indels: ${gt.indels} bp
      </div>
    `;
  } else {
    gtBox.style.display = 'none';
  }
}

// Tab 2: Evolutionary Risk Calculator
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
        <div style="margin-bottom: 0.35rem;">${data.guidance}</div>
        <div style="font-size: 0.8rem; color: var(--text-muted);">
          Dynamic Homology Weight: <code style="color: var(--primary);">w_hom(${t}) = ${data.homology_weight_w_hom}</code>
          ${data.is_exceeding_t_crit ? ' — <strong style="color: #EF4444;">VƯỢT NGƯỠNG TỚI HẠN (w_hom = 0)</strong>' : ''}
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

  const w = 700;
  const h = 220;
  const padLeft = 60;
  const padRight = 30;
  const padTop = 25;
  const padBottom = 35;

  const plotW = w - padLeft - padRight;
  const plotH = h - padTop - padBottom;

  const minMya = 0;
  const maxMya = 90;
  const minAcc = 50;
  const maxAcc = 100;

  const scaleX = (m) => padLeft + ((m - minMya) / (maxMya - minMya)) * plotW;
  const scaleY = (acc) => padTop + plotH - ((acc - minAcc) / (maxAcc - minAcc)) * plotH;

  // Build SVG Content
  let svgContent = '';

  // Background Grid Lines
  for (let acc = 60; acc <= 100; acc += 10) {
    const y = scaleY(acc);
    svgContent += `<line x1="${padLeft}" y1="${y}" x2="${w - padRight}" y2="${y}" stroke="rgba(255,255,255,0.08)" stroke-dasharray="3,3" />`;
    svgContent += `<text x="${padLeft - 10}" y="${y + 4}" fill="#64748B" font-size="11" text-anchor="end" font-family="JetBrains Mono">${acc}%</text>`;
  }

  for (let m = 0; m <= 90; m += 15) {
    const x = scaleX(m);
    svgContent += `<line x1="${x}" y1="${padTop}" x2="${x}" y2="${h - padBottom}" stroke="rgba(255,255,255,0.06)" stroke-dasharray="3,3" />`;
    svgContent += `<text x="${x}" y="${h - padBottom + 18}" fill="#64748B" font-size="11" text-anchor="middle" font-family="JetBrains Mono">${m} Mya</text>`;
  }

  // Red Danger Zone Shading (t >= tCrit)
  const xCrit = scaleX(tCrit);
  const dangerWidth = (w - padRight) - xCrit;
  if (dangerWidth > 0) {
    svgContent += `
      <rect x="${xCrit}" y="${padTop}" width="${dangerWidth}" height="${plotH}" fill="rgba(239, 68, 68, 0.08)" />
      <text x="${xCrit + 10}" y="${padTop + 20}" fill="#F87171" font-size="11" font-weight="600" font-family="Outfit">⚠️ VÙNG NGUY HIỂM (FCR > 90%)</text>
    `;
  }

  // Vertical Line at tCrit
  svgContent += `
    <line x1="${xCrit}" y1="${padTop}" x2="${xCrit}" y2="${h - padBottom}" stroke="#EF4444" stroke-width="1.8" stroke-dasharray="4,4" />
    <text x="${xCrit}" y="${padTop - 8}" fill="#EF4444" font-size="11" font-weight="700" text-anchor="middle" font-family="JetBrains Mono">t_crit = ${tCrit} Mya</text>
  `;

  // Theoretical Exponential Decay Curve
  const lambda = 0.006395;
  let pathD = '';
  for (let m = 0; m <= 90; m += 1) {
    const acc = Math.exp(-lambda * m) * 100.0;
    const x = scaleX(m);
    const y = scaleY(acc);
    if (m === 0) pathD += `M ${x} ${y}`;
    else pathD += ` L ${x} ${y}`;
  }
  svgContent += `<path d="${pathD}" fill="none" stroke="#38BDF8" stroke-width="3" />`;

  // Current Position Highlight
  const curX = scaleX(currentMya);
  const curY = scaleY(currentIdentity);
  svgContent += `
    <line x1="${curX}" y1="${padTop}" x2="${curX}" y2="${h - padBottom}" stroke="#F59E0B" stroke-width="1.5" />
    <circle cx="${curX}" cy="${curY}" r="7" fill="#F59E0B" stroke="#0F172A" stroke-width="2" />
    <circle cx="${curX}" cy="${curY}" r="12" fill="none" stroke="rgba(245, 158, 11, 0.5)" stroke-width="1.5">
      <animate attributeName="r" values="7;14;7" dur="2s" repeatCount="indefinite"/>
      <animate attributeName="opacity" values="0.8;0.2;0.8" dur="2s" repeatCount="indefinite"/>
    </circle>
    <text x="${curX + (currentMya > 70 ? -12 : 12)}" y="${curY - 12}" fill="#F59E0B" font-size="12" font-weight="700" font-family="JetBrains Mono" text-anchor="${currentMya > 70 ? 'end' : 'start'}">
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
