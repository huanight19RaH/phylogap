/**
 * PhyloMorph 3D Studio - Evolutionary Genomic Synthesis & Particle Morphing Engine.
 * Features:
 * 1. Two-Way 3D Particle Morphing: Phenotype (Avian Mesh) <-> Genotype (DNA Double Helix).
 * 2. Parametric DNA Helix with nucleotide base-pair color mapping (A: Crimson, T: Gold, G: Emerald, C: Cyan).
 * 3. Five Procedural Avian Phenotypes: Hummingbird, Finch, Junglefowl, Dodo, Archaeopteryx.
 * 4. Cross-Species In-Silico Hybridization & Ancestral Sequence Prediction.
 * 5. Interactive 3D Gap Healing Simulator with dynamic donor particle injection.
 *
 * Technology: Three.js WebGL Particle Shaders with GPU-accelerated vortex morphing.
 * Author: PhyloGap Research & Product Team
 * Date: 2026-10-28
 */

class PhyloMorphEngine {
  constructor(containerId) {
    this.container = document.getElementById(containerId);
    if (!this.container) return;

    this.particleCount = 6400;
    this.currentSpecies = 'c_anna';
    this.morphProgress = 0.0; // 0.0 = Phenotype, 1.0 = Double Helix
    this.targetMorphProgress = 0.0;
    this.isGapActive = false;
    this.isHealing = false;

    this.speciesProfiles = {
      c_anna: { name: "Calypte anna", common: "Anna's Hummingbird", morph: "hummingbird", color: [0.0, 0.9, 0.9] },
      t_guttata: { name: "Taeniopygia guttata", common: "Zebra Finch", morph: "finch", color: [1.0, 0.6, 0.2] },
      g_gallus: { name: "Gallus gallus", common: "Red Junglefowl", morph: "fowl", color: [1.0, 0.2, 0.2] },
      r_cucullatus: { name: "Raphus cucullatus", common: "Dodo (Extinct)", morph: "dodo", color: [0.3, 0.8, 0.4] },
      archaeopteryx: { name: "Archaeopteryx lithographica", common: "Urvogel", morph: "archaeopteryx", color: [0.8, 0.3, 0.9] }
    };

    this.initScene();
    this.generateParticleAttributes();
    this.createParticleSystem();
    this.bindEvents();
    this.animate();
  }

  initScene() {
    const width = this.container.clientWidth || 800;
    const height = this.container.clientHeight || 550;

    this.scene = new THREE.Scene();
    this.scene.fog = new THREE.FogExp2(0x060913, 0.015);

    this.camera = new THREE.PerspectiveCamera(45, width / height, 0.1, 1000);
    this.camera.position.set(0, 5, 28);

    this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    this.renderer.setSize(width, height);
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    this.container.appendChild(this.renderer.domElement);

    if (window.THREE && THREE.OrbitControls) {
      this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
      this.controls.enableDamping = true;
      this.controls.dampingFactor = 0.05;
      this.controls.autoRotate = true;
      this.controls.autoRotateSpeed = 0.8;
      this.controls.maxDistance = 60;
      this.controls.minDistance = 5;
    }

    this.clock = new THREE.Clock();
  }

  createSoftGlowTexture() {
    const canvas = document.createElement('canvas');
    canvas.width = 64;
    canvas.height = 64;
    const ctx = canvas.getContext('2d');
    const grad = ctx.createRadialGradient(32, 32, 0, 32, 32, 32);
    grad.addColorStop(0, 'rgba(255, 255, 255, 1)');
    grad.addColorStop(0.3, 'rgba(255, 255, 255, 0.8)');
    grad.addColorStop(0.7, 'rgba(255, 255, 255, 0.15)');
    grad.addColorStop(1, 'rgba(255, 255, 255, 0)');
    ctx.fillStyle = grad;
    ctx.fillRect(0, 0, 64, 64);
    return new THREE.CanvasTexture(canvas);
  }

  generateDnaPositions(count) {
    const pos = new Float32Array(count * 3);
    const colors = new Float32Array(count * 3);
    const heightSpan = 22.0;
    const radius = 3.5;
    const turns = 4.5;

    // Nucleotide Palette: A: Crimson, T: Gold, G: Emerald, C: Cyan, N: Pulsing Red
    const basePalette = [
      [1.0, 0.2, 0.4], // A (Crimson)
      [1.0, 0.8, 0.1], // T (Gold)
      [0.0, 0.95, 0.55], // G (Emerald)
      [0.0, 0.8, 1.0]  // C (Cyan)
    ];

    const backboneCount = Math.floor(count * 0.45);
    const rungsCount = count - backboneCount;

    // 1. Double Helix Backbones
    for (let i = 0; i < backboneCount; i++) {
      const t = (i / backboneCount) * Math.PI * 2 * turns;
      const strand = (i % 2 === 0) ? 1.0 : -1.0;
      const y = ((i / backboneCount) - 0.5) * heightSpan;
      const x = Math.cos(t) * radius * strand;
      const z = Math.sin(t) * radius * strand;

      const idx = i * 3;
      pos[idx] = x;
      pos[idx + 1] = y;
      pos[idx + 2] = z;

      // Backbone electric cyan/magenta
      if (strand > 0) {
        colors[idx] = 0.0; colors[idx + 1] = 0.85; colors[idx + 2] = 1.0;
      } else {
        colors[idx] = 0.95; colors[idx + 1] = 0.2; colors[idx + 2] = 0.85;
      }
    }

    // 2. Base-Pair Rungs
    const offset = backboneCount * 3;
    const basePairSteps = 36;
    for (let i = 0; i < rungsCount; i++) {
      const stepIdx = Math.floor((i / rungsCount) * basePairSteps);
      const frac = (stepIdx / basePairSteps);
      const theta = frac * Math.PI * 2 * turns;
      const y = (frac - 0.5) * heightSpan;

      const p1 = new THREE.Vector3(Math.cos(theta) * radius, y, Math.sin(theta) * radius);
      const p2 = new THREE.Vector3(Math.cos(theta + Math.PI) * radius, y, Math.sin(theta + Math.PI) * radius);

      // Lerp position along rung
      const rungAlpha = (i % (rungsCount / basePairSteps)) / (rungsCount / basePairSteps);
      const pt = new THREE.Vector3().lerpVectors(p1, p2, rungAlpha);

      const idx = offset + i * 3;
      pos[idx] = pt.x;
      pos[idx + 1] = pt.y;
      pos[idx + 2] = pt.z;

      // Base color assignment
      const baseCol = basePalette[(stepIdx * 3 + (rungAlpha > 0.5 ? 1 : 0)) % 4];
      colors[idx] = baseCol[0];
      colors[idx + 1] = baseCol[1];
      colors[idx + 2] = baseCol[2];
    }

    return { positions: pos, colors: colors };
  }

  generateAvianMorphPositions(morphType, count) {
    const pos = new Float32Array(count * 3);
    const rng = Math.random;

    for (let i = 0; i < count; i++) {
      let x = 0, y = 0, z = 0;
      const p = rng();

      if (morphType === 'hummingbird') {
        // Slender torpedo body + needle bill + swept back wings
        if (p < 0.35) { // Torpedo body
          const u = rng() * Math.PI * 2;
          const v = (rng() - 0.5) * Math.PI;
          const r = 2.0;
          x = r * 0.7 * Math.cos(v) * Math.cos(u);
          y = (r * 1.5 * Math.sin(v)) - 0.5;
          z = r * 0.7 * Math.cos(v) * Math.sin(u);
        } else if (p < 0.50) { // Long needle bill
          const len = rng() * 4.5;
          x = (rng() - 0.5) * 0.15;
          y = 1.8 + len * 0.6;
          z = 0.8 + len * 0.8;
        } else if (p < 0.85) { // Swept-back wings
          const side = rng() > 0.5 ? 1 : -1;
          const wingSpan = rng() * 8.0;
          x = side * (1.2 + wingSpan);
          y = 0.5 - (wingSpan * 0.3) + (rng() - 0.5) * 0.4;
          z = -1.0 - (wingSpan * 0.7) + (rng() - 0.5) * 0.4;
        } else { // Tail feathers
          const tailL = rng() * 4.0;
          x = (rng() - 0.5) * (1.0 + tailL * 0.3);
          y = -2.0 - tailL;
          z = -1.5 - tailL * 0.4;
        }
      } else if (morphType === 'finch') {
        // Plump rounded body + short conical bill + perching wings
        if (p < 0.55) { // Plump round body
          const u = rng() * Math.PI * 2;
          const v = (rng() - 0.5) * Math.PI;
          const r = 2.8;
          x = r * 1.0 * Math.cos(v) * Math.cos(u);
          y = r * 1.1 * Math.sin(v);
          z = r * 0.9 * Math.cos(v) * Math.sin(u);
        } else if (p < 0.65) { // Short conical seed beak
          const len = rng() * 1.6;
          x = (rng() - 0.5) * (1.2 - len * 0.6);
          y = 1.5 + len * 0.5;
          z = 2.0 + len * 0.9;
        } else if (p < 0.88) { // Folded compact wings
          const side = rng() > 0.5 ? 1 : -1;
          const w = rng() * 3.5;
          x = side * (2.2 + w * 0.4);
          y = 0.5 - w * 0.8;
          z = -0.5 - w * 0.6;
        } else { // Short tail
          const tL = rng() * 2.5;
          x = (rng() - 0.5) * 0.8;
          y = -2.5 - tL * 0.8;
          z = -1.8 - tL * 0.6;
        }
      } else if (morphType === 'fowl') {
        // Proud upright rooster body + arched tail + head comb
        if (p < 0.45) { // Heavy chest and body
          const u = rng() * Math.PI * 2;
          const v = (rng() - 0.5) * Math.PI;
          const r = 3.2;
          x = r * 0.85 * Math.cos(v) * Math.cos(u);
          y = r * 1.3 * Math.sin(v) - 0.5;
          z = r * 0.85 * Math.cos(v) * Math.sin(u);
        } else if (p < 0.58) { // Head comb and wattles
          const cH = rng() * 2.2;
          x = (rng() - 0.5) * 0.35;
          y = 3.2 + cH;
          z = 0.8 + Math.sin(cH * 3.0) * 0.8;
        } else if (p < 0.85) { // High arching sickle tail
          const tArc = rng() * Math.PI * 0.7;
          x = (rng() - 0.5) * 1.5;
          y = 0.5 + Math.sin(tArc) * 5.0;
          z = -2.0 - Math.cos(tArc) * 4.5;
        } else { // Legs and spurs
          const side = rng() > 0.5 ? 1 : -1;
          const lH = rng() * 3.5;
          x = side * 1.2;
          y = -3.5 - lH;
          z = (rng() - 0.5) * 0.5;
        }
      } else if (morphType === 'dodo') {
        // Giant rotund body + heavy hooked beak + tiny vestigial wings
        if (p < 0.65) { // Rotund heavy torso
          const u = rng() * Math.PI * 2;
          const v = (rng() - 0.5) * Math.PI;
          const r = 4.0;
          x = r * 1.1 * Math.cos(v) * Math.cos(u);
          y = r * 1.0 * Math.sin(v) - 0.5;
          z = r * 1.1 * Math.cos(v) * Math.sin(u);
        } else if (p < 0.80) { // Large bulbous hooked beak
          const bL = rng() * 3.5;
          x = (rng() - 0.5) * 1.0;
          y = 2.0 - Math.pow(bL / 3.5, 2) * 1.8;
          z = 2.5 + bL * 1.2;
        } else { // Stubby wings
          const side = rng() > 0.5 ? 1 : -1;
          const wL = rng() * 1.8;
          x = side * (3.5 + wL * 0.3);
          y = -wL * 0.5;
          z = -wL * 0.4;
        }
      } else { // Archaeopteryx
        // Long reptilian feathered tail + wing claws + toothed jaw
        if (p < 0.35) { // Sleek dinosaurian body
          const u = rng() * Math.PI * 2;
          const v = (rng() - 0.5) * Math.PI;
          const r = 2.2;
          x = r * 0.7 * Math.cos(v) * Math.cos(u);
          y = r * 1.1 * Math.sin(v);
          z = r * 0.7 * Math.cos(v) * Math.sin(u);
        } else if (p < 0.65) { // Long multi-vertebrae feathered tail
          const tailLen = rng() * 8.5;
          x = (rng() - 0.5) * (0.4 + Math.sin(tailLen * 0.8) * 1.2);
          y = -1.5 - tailLen * 0.7;
          z = -1.8 - tailLen * 1.1;
        } else { // Broad feathered wings with claws
          const side = rng() > 0.5 ? 1 : -1;
          const w = rng() * 6.5;
          x = side * (1.5 + w);
          y = 0.8 - w * 0.2;
          z = -w * 0.5;
        }
      }

      const idx = i * 3;
      pos[idx] = x;
      pos[idx + 1] = y;
      pos[idx + 2] = z;
    }

    return pos;
  }

  generateParticleAttributes() {
    this.dnaData = this.generateDnaPositions(this.particleCount);
    this.organismPos = this.generateAvianMorphPositions(this.speciesProfiles[this.currentSpecies].morph, this.particleCount);

    this.geometry = new THREE.BufferGeometry();
    this.geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(this.organismPos), 3));
    this.geometry.setAttribute('aDnaPosition', new THREE.BufferAttribute(this.dnaData.positions, 3));
    this.geometry.setAttribute('aDnaColor', new THREE.BufferAttribute(this.dnaData.colors, 3));

    // Base organism colors (iridescent cyan/amber)
    const orgColors = new Float32Array(this.particleCount * 3);
    const spCol = this.speciesProfiles[this.currentSpecies].color;
    for (let i = 0; i < this.particleCount; i++) {
      const idx = i * 3;
      orgColors[idx] = spCol[0] + (Math.random() - 0.5) * 0.15;
      orgColors[idx + 1] = spCol[1] + (Math.random() - 0.5) * 0.15;
      orgColors[idx + 2] = spCol[2] + (Math.random() - 0.5) * 0.15;
    }
    this.geometry.setAttribute('aOrgColor', new THREE.BufferAttribute(orgColors, 3));

    // Random phases for organic particle floating
    const phases = new Float32Array(this.particleCount);
    for (let i = 0; i < this.particleCount; i++) phases[i] = Math.random() * Math.PI * 2;
    this.geometry.setAttribute('aPhase', new THREE.BufferAttribute(phases, 1));
  }

  createParticleSystem() {
    const vertexShader = `
      uniform float uProgress;
      uniform float uTime;
      uniform float uNoise;
      uniform float uGapAlarm;
      uniform float uHealing;
      
      attribute vec3 aDnaPosition;
      attribute vec3 aDnaColor;
      attribute vec3 aOrgColor;
      attribute float aPhase;
      
      varying vec3 vColor;
      varying float vAlpha;
      
      void main() {
        // Organic vortex interpolation
        float p = smoothstep(0.0, 1.0, uProgress);
        
        // Swirling particle vortex during morphing
        float vortex = sin(p * 3.14159) * 2.8;
        float angle = uTime * 1.5 + aDnaPosition.y * 0.4;
        vec3 vortexOffset = vec3(sin(angle) * vortex, 0.0, cos(angle) * vortex);
        
        vec3 mixedPos = mix(position, aDnaPosition, p) + vortexOffset;
        
        // Idle breathing wave
        mixedPos.y += sin(uTime * 2.0 + aPhase) * 0.08;
        mixedPos.x += cos(uTime * 1.5 + aPhase) * 0.05;
        
        // Color transition
        vec3 col = mix(aOrgColor, aDnaColor, p);
        
        // Gap Alarm (Pulsing Red Break at middle locus)
        if (p > 0.4 && uGapAlarm > 0.1) {
          if (abs(aDnaPosition.y) < 2.5) {
            float pulse = 0.5 + 0.5 * sin(uTime * 10.0);
            col = mix(col, vec3(1.0, 0.05, 0.1), pulse * uGapAlarm);
            mixedPos.x += (sin(uTime * 20.0 + aPhase) * 0.4) * uGapAlarm;
          }
        }
        
        // Healing green stream
        if (uHealing > 0.1) {
          if (abs(aDnaPosition.y) < 3.0) {
            col = mix(col, vec3(0.0, 1.0, 0.6), uHealing);
          }
        }
        
        vColor = col;
        vAlpha = 0.85;
        
        vec4 mvPosition = modelViewMatrix * vec4(mixedPos, 1.0);
        gl_Position = projectionMatrix * mvPosition;
        
        // Perspective particle sizing
        gl_PointSize = (18.0 / -mvPosition.z) * (1.0 + sin(aPhase + uTime) * 0.2);
      }
    `;

    const fragmentShader = `
      uniform sampler2D uTexture;
      varying vec3 vColor;
      varying float vAlpha;
      
      void main() {
        vec4 tex = texture2D(uTexture, gl_PointCoord);
        if (tex.a < 0.05) discard;
        gl_FragColor = vec4(vColor, vAlpha * tex.a);
      }
    `;

    this.material = new THREE.ShaderMaterial({
      vertexShader: vertexShader,
      fragmentShader: fragmentShader,
      uniforms: {
        uProgress: { value: 0.0 },
        uTime: { value: 0.0 },
        uNoise: { value: 1.0 },
        uGapAlarm: { value: 0.0 },
        uHealing: { value: 0.0 },
        uTexture: { value: this.createSoftGlowTexture() }
      },
      transparent: true,
      depthWrite: false,
      blending: THREE.AdditiveBlending
    });

    this.particleMesh = new THREE.Points(this.geometry, this.material);
    this.scene.add(this.particleMesh);
  }

  setMorphProgress(progress) {
    this.targetMorphProgress = Math.max(0.0, Math.min(1.0, progress));
  }

  setSpecies(speciesId) {
    if (!this.speciesProfiles[speciesId]) return;
    this.currentSpecies = speciesId;
    const morphType = this.speciesProfiles[speciesId].morph;

    // Recalculate organism positions
    const newPositions = this.generateAvianMorphPositions(morphType, this.particleCount);
    const posAttr = this.geometry.attributes.position;
    posAttr.copyArray(newPositions);
    posAttr.needsUpdate = true;

    // Update base organism color
    const spCol = this.speciesProfiles[speciesId].color;
    const orgColorAttr = this.geometry.attributes.aOrgColor;
    for (let i = 0; i < this.particleCount; i++) {
      const idx = i * 3;
      orgColorAttr.array[idx] = spCol[0] + (Math.random() - 0.5) * 0.15;
      orgColorAttr.array[idx + 1] = spCol[1] + (Math.random() - 0.5) * 0.15;
      orgColorAttr.array[idx + 2] = spCol[2] + (Math.random() - 0.5) * 0.15;
    }
    orgColorAttr.needsUpdate = true;
  }

  setHybridMorph(speciesA, speciesB, ratioB) {
    const morphA = this.generateAvianMorphPositions(this.speciesProfiles[speciesA].morph, this.particleCount);
    const morphB = this.generateAvianMorphPositions(this.speciesProfiles[speciesB].morph, this.particleCount);

    const posAttr = this.geometry.attributes.position;
    for (let i = 0; i < this.particleCount * 3; i++) {
      posAttr.array[i] = THREE.MathUtils.lerp(morphA[i], morphB[i], ratioB);
    }
    posAttr.needsUpdate = true;

    // Blend Colors
    const colA = this.speciesProfiles[speciesA].color;
    const colB = this.speciesProfiles[speciesB].color;
    const orgColorAttr = this.geometry.attributes.aOrgColor;
    for (let i = 0; i < this.particleCount; i++) {
      const idx = i * 3;
      orgColorAttr.array[idx] = THREE.MathUtils.lerp(colA[0], colB[0], ratioB);
      orgColorAttr.array[idx + 1] = THREE.MathUtils.lerp(colA[1], colB[1], ratioB);
      orgColorAttr.array[idx + 2] = THREE.MathUtils.lerp(colA[2], colB[2], ratioB);
    }
    orgColorAttr.needsUpdate = true;
  }

  triggerGapAlarm(enable = true) {
    this.isGapActive = enable;
  }

  triggerHealingPulse() {
    this.isHealing = true;
    setTimeout(() => {
      this.isHealing = false;
      this.isGapActive = false;
    }, 2200);
  }

  bindEvents() {
    window.addEventListener('resize', () => {
      if (!this.container) return;
      const width = this.container.clientWidth;
      const height = this.container.clientHeight;
      if (width === 0 || height === 0) return;
      this.camera.aspect = width / height;
      this.camera.updateProjectionMatrix();
      this.renderer.setSize(width, height);
    });
  }

  animate() {
    requestAnimationFrame(() => this.animate());

    const delta = this.clock.getDelta();
    const elapsedTime = this.clock.getElapsedTime();

    // Smooth Lerp for morph progress
    this.morphProgress = THREE.MathUtils.lerp(this.morphProgress, this.targetMorphProgress, delta * 3.5);
    this.material.uniforms.uProgress.value = this.morphProgress;
    this.material.uniforms.uTime.value = elapsedTime;

    // Smooth Gap and Healing Transitions
    const targetGap = this.isGapActive ? 1.0 : 0.0;
    this.material.uniforms.uGapAlarm.value = THREE.MathUtils.lerp(this.material.uniforms.uGapAlarm.value, targetGap, delta * 6.0);

    const targetHealing = this.isHealing ? 1.0 : 0.0;
    this.material.uniforms.uHealing.value = THREE.MathUtils.lerp(this.material.uniforms.uHealing.value, targetHealing, delta * 5.0);

    if (this.controls) this.controls.update();
    this.renderer.render(this.scene, this.camera);
  }
}

// Global hook
window.PhyloMorphEngine = PhyloMorphEngine;
