/**
 * PhyloMorph 3D Studio - High-Luminescence Particle Morphing Engine.
 * Inspired by ThreeUI and BlendSwap Particle Simulations.
 *
 * Core Features:
 * 1. Self-contained mathematical GLSL circular particle shader:
 *    - Anti-aliased circular aperture (zero texture dependencies, zero WebGL binding failures).
 *    - Multi-tier radiant energy: soft atmospheric halo + dense inner glow + white-hot nucleus.
 * 2. 12,000 High-Density GPU particles with helical vortex morphing.
 * 3. Sculpted 3D Models:
 *    - Avian (Calypte anna - Hummingbird with needle bill, gorget, arched swept wings with wing-beat motion, and fanned tail).
 *    - Mammalia (Quadruped vertebrate with articulated spine, 4 limbs with paws, skull, muzzle, ears, and curved tail).
 *    - Viridiplantae (Arabidopsis thaliana - 10 radiating rosette leaves with venation, cauline stem, and apical inflorescence).
 *    - Finch (Taeniopygia guttata - Plump songbird with conical seed beak).
 *    - DNA Double Helix (Dual antiparallel backbones, 60 base-pair rungs in A-T, G-C, and ambient ion dust).
 * 4. Interactive HeroDnaBackground:
 *    - Autonomous 3D floating DNA helix with mouse drag interaction, glowing base pairs, and ResizeObserver.
 * 5. Camera Presets (Isometric 3/4, Front, Side, Top) and auto-rotation toggle.
 *
 * Author: PhyloGap Research & Product Team
 * Date: 2026-10-28
 */

class PhyloMorphEngine {
  constructor(containerId) {
    this.container = document.getElementById(containerId);
    if (!this.container) {
      console.warn(`PhyloMorphEngine: Container #${containerId} not found.`);
      return;
    }

    this.particleCount = 12000;
    this.currentSpecies = 'c_anna';
    this.morphProgress = 0.0; // 0.0 = Organism Phenotype, 1.0 = DNA Helix
    this.targetMorphProgress = 0.0;
    this.isGapActive = false;
    this.isHealing = false;
    this.autoRotate = true;
    this.particleBaseSize = 85.0;

    this.speciesProfiles = {
      c_anna: {
        name: "Calypte anna",
        common: "Anna's Hummingbird (Chim Ruồi)",
        morph: "hummingbird",
        color: [0.06, 0.95, 0.65], // Bio Emerald
        secondaryColor: [0.08, 0.85, 1.0] // Cyan
      },
      mammal: {
        name: "Mammalia (Thú)",
        common: "Vertebrate Model (Động Vật Có Vú)",
        morph: "mammal",
        color: [0.2, 0.92, 0.55], // Electric Mint
        secondaryColor: [0.05, 0.75, 0.95]
      },
      plant: {
        name: "Arabidopsis thaliana",
        common: "Thale Cress (Thực Vật)",
        morph: "plant",
        color: [0.3, 1.0, 0.45], // Botanical Neon Green
        secondaryColor: [0.1, 0.9, 0.75]
      },
      t_guttata: {
        name: "Taeniopygia guttata",
        common: "Zebra Finch (Chim Sẻ)",
        morph: "finch",
        color: [0.08, 0.95, 0.8], // Cyan-Mint
        secondaryColor: [0.25, 1.0, 0.55]
      }
    };

    this.initScene();
    this.generateParticleAttributes();
    this.createParticleSystem();
    this.bindEvents();
    this.animate();
  }

  initScene() {
    const width = this.container.clientWidth || 800;
    const height = this.container.clientHeight || 560;

    this.scene = new THREE.Scene();
    this.scene.fog = new THREE.FogExp2(0x030303, 0.010);

    this.camera = new THREE.PerspectiveCamera(45, width / height, 0.1, 1000);
    // Dynamic 3/4 isometric perspective so models show full 3D silhouette
    this.camera.position.set(18, 9, 24);

    this.renderer = new THREE.WebGLRenderer({
      antialias: true,
      alpha: true,
      powerPreference: "high-performance"
    });
    this.renderer.setSize(width, height);
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    this.renderer.domElement.style.width = '100%';
    this.renderer.domElement.style.height = '100%';
    this.renderer.domElement.style.display = 'block';

    // Clear any previous canvas
    while (this.container.firstChild && this.container.firstChild.tagName === 'CANVAS') {
      this.container.removeChild(this.container.firstChild);
    }
    this.container.appendChild(this.renderer.domElement);

    if (window.THREE && THREE.OrbitControls) {
      this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
      this.controls.enableDamping = true;
      this.controls.dampingFactor = 0.05;
      this.controls.autoRotate = this.autoRotate;
      this.controls.autoRotateSpeed = 0.8;
      this.controls.maxDistance = 75;
      this.controls.minDistance = 6;
      this.controls.target.set(0, 0, 0);
    }

    this.clock = new THREE.Clock();

    // ResizeObserver ensures flawless sizing on all displays
    if (window.ResizeObserver) {
      this.resizeObserver = new ResizeObserver(() => this.handleResize());
      this.resizeObserver.observe(this.container);
    }
  }

  handleResize() {
    if (!this.container || !this.renderer || !this.camera) return;
    const w = this.container.clientWidth;
    const h = this.container.clientHeight;
    if (w === 0 || h === 0) return;
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
    this.renderer.setSize(w, h);
  }

  generateDnaPositions(count) {
    const pos = new Float32Array(count * 3);
    const colors = new Float32Array(count * 3);
    const heightSpan = 24.0;
    const radius = 4.4;
    const turns = 5.0;

    // Vibrant Base Colors:
    // A: Bio-Emerald, T: Cyber Cyan, G: Electric Mint, C: Deep Aqua
    const basePalette = [
      [0.06, 0.95, 0.65], // A (Emerald)
      [0.05, 0.85, 1.0],  // T (Cyan)
      [0.25, 1.0, 0.50],  // G (Mint)
      [0.10, 0.70, 0.95]   // C (Aqua)
    ];

    const backboneCount = Math.floor(count * 0.45);
    const rungsCount = Math.floor(count * 0.45);
    const ionCount = count - backboneCount - rungsCount;

    // 1. Double Helix Backbones (Antiparallel double strands)
    for (let i = 0; i < backboneCount; i++) {
      const frac = i / backboneCount;
      const theta = frac * Math.PI * 2 * turns;
      const strand = (i % 2 === 0) ? 1.0 : -1.0;
      const y = (frac - 0.5) * heightSpan;

      const rJitter = radius + (Math.sin(i * 12.0) * 0.18);
      const x = Math.cos(theta) * rJitter * strand;
      const z = Math.sin(theta) * rJitter * strand;

      const idx = i * 3;
      pos[idx] = x;
      pos[idx + 1] = y;
      pos[idx + 2] = z;

      if (strand > 0) {
        colors[idx] = 0.06; colors[idx + 1] = 0.95; colors[idx + 2] = 0.65; // Emerald
      } else {
        colors[idx] = 0.05; colors[idx + 1] = 0.85; colors[idx + 2] = 1.0; // Cyan
      }
    }

    // 2. Base-Pair Cross Rungs
    const offset = backboneCount * 3;
    const baseSteps = 60; // 60 paired rungs
    const particlesPerRung = Math.floor(rungsCount / baseSteps);

    for (let s = 0; s < baseSteps; s++) {
      const frac = s / baseSteps;
      const theta = frac * Math.PI * 2 * turns;
      const y = (frac - 0.5) * heightSpan;

      const p1 = new THREE.Vector3(Math.cos(theta) * radius, y, Math.sin(theta) * radius);
      const p2 = new THREE.Vector3(Math.cos(theta + Math.PI) * radius, y, Math.sin(theta + Math.PI) * radius);

      const colLeft = basePalette[(s * 2) % 4];
      const colRight = basePalette[(s * 2 + 1) % 4];

      for (let j = 0; j < particlesPerRung; j++) {
        const rungIdx = s * particlesPerRung + j;
        if (rungIdx >= rungsCount) break;

        const alpha = j / particlesPerRung;
        // Hydrogen bond gap at center
        let ptAlpha = alpha;
        if (Math.abs(alpha - 0.5) < 0.05) {
          ptAlpha = alpha < 0.5 ? 0.44 : 0.56;
        }
        const pt = new THREE.Vector3().lerpVectors(p1, p2, ptAlpha);

        const idx = offset + rungIdx * 3;
        pos[idx] = pt.x + (Math.random() - 0.5) * 0.12;
        pos[idx + 1] = pt.y + (Math.random() - 0.5) * 0.12;
        pos[idx + 2] = pt.z + (Math.random() - 0.5) * 0.12;

        const c = alpha < 0.5 ? colLeft : colRight;
        colors[idx] = c[0];
        colors[idx + 1] = c[1];
        colors[idx + 2] = c[2];
      }
    }

    // 3. Ambient Molecular Ion Cloud
    const ionOffset = (backboneCount + rungsCount) * 3;
    for (let k = 0; k < ionCount; k++) {
      const u = Math.random() * Math.PI * 2;
      const r = radius * (1.1 + Math.random() * 1.5);
      const y = (Math.random() - 0.5) * (heightSpan + 4.0);

      const idx = ionOffset + k * 3;
      pos[idx] = Math.cos(u) * r;
      pos[idx + 1] = y;
      pos[idx + 2] = Math.sin(u) * r;

      colors[idx] = 0.2;
      colors[idx + 1] = 0.95;
      colors[idx + 2] = 0.85;
    }

    return { positions: pos, colors: colors };
  }

  generateOrganismPositions(morphType, count) {
    const pos = new Float32Array(count * 3);
    const rng = Math.random;

    for (let i = 0; i < count; i++) {
      let x = 0, y = 0, z = 0;
      const p = rng();

      if (morphType === 'hummingbird') {
        // Anatomical Hummingbird (Calypte anna)
        // Oriented diagonally at a 35° angle for dramatic 3D silhouette
        if (p < 0.28) {
          // Torpedo Fuselage & Chest
          const u = rng() * Math.PI * 2;
          const v = (rng() - 0.5) * Math.PI;
          const r = 2.4;
          const lx = r * 0.75 * Math.cos(v) * Math.cos(u);
          const ly = (r * 1.6 * Math.sin(v)) - 0.4;
          const lz = r * 0.85 * Math.cos(v) * Math.sin(u);
          x = lx; y = ly; z = lz;
        } else if (p < 0.42) {
          // Head, Gorget & Needle Beak
          if (rng() < 0.38) {
            // Head cranium
            const u = rng() * Math.PI * 2;
            const v = (rng() - 0.5) * Math.PI;
            x = 0.95 * Math.cos(v) * Math.cos(u);
            y = 2.4 + 0.95 * Math.sin(v);
            z = 1.2 + 0.95 * Math.cos(v) * Math.sin(u);
          } else {
            // Long needle-thin beak pointing forward and slightly downward
            const bL = rng() * 5.2;
            const taper = 1.0 - (bL / 5.2) * 0.7;
            x = (rng() - 0.5) * 0.12 * taper;
            y = 2.3 + bL * 0.38;
            z = 2.0 + bL * 1.1;
          }
        } else if (p < 0.82) {
          // Swept-back articulated wings with primary flight feathers
          const side = rng() > 0.5 ? 1.0 : -1.0;
          const span = rng(); // 0 (wing root) to 1 (wing tip)
          const wingLen = span * 9.2;
          const chord = (1.0 - span * 0.7) * 2.4 * (rng() - 0.5);
          const camber = Math.sin(span * Math.PI) * 1.4;

          x = side * (1.1 + wingLen);
          y = 0.8 + camber - (span * 1.8) + chord * 0.25;
          z = -0.4 - (span * 4.2) + chord;
        } else {
          // Fanned rectrix tail feathers
          const side = (rng() - 0.5) * 2.0;
          const tL = rng() * 4.8;
          x = side * (tL * 0.45);
          y = -2.8 - tL * 0.85;
          z = -1.8 - tL * 0.75;
        }
      } else if (morphType === 'mammal') {
        // Anatomical Mammalia (Vertebrate Quadruped)
        if (p < 0.30) {
          // Arched Spine & Thoracic Ribcage
          const sT = (rng() - 0.5) * 8.6; // -4.3 to +4.3 along z
          const arch = Math.cos((sT / 4.3) * (Math.PI / 2.3)) * 1.6;
          const ribR = 2.1 * Math.sin(rng() * Math.PI);
          const ribA = rng() * Math.PI * 2;
          x = Math.cos(ribA) * ribR * 0.82;
          y = arch + Math.sin(ribA) * ribR * 0.88;
          z = sT;
        } else if (p < 0.46) {
          // Head, Muzzle & Upright Ears
          const hP = rng();
          if (hP < 0.45) {
            // Skull cranium
            const u = rng() * Math.PI * 2;
            const v = (rng() - 0.5) * Math.PI;
            x = 1.15 * Math.cos(v) * Math.cos(u);
            y = 2.0 + 1.15 * Math.sin(v);
            z = 4.4 + 1.15 * Math.cos(v) * Math.sin(u);
          } else if (hP < 0.75) {
            // Elongated Snout/Muzzle
            const mFrac = rng();
            x = (rng() - 0.5) * 0.85 * (1.0 - mFrac * 0.4);
            y = 1.5 + (rng() - 0.5) * 0.6;
            z = 5.4 + mFrac * 1.6;
          } else {
            // Upright Triangular Ears
            const side = rng() > 0.5 ? 1 : -1;
            const eH = rng() * 1.4;
            x = side * (0.85 + eH * 0.35);
            y = 2.8 + eH;
            z = 4.2 - eH * 0.25;
          }
        } else if (p < 0.84) {
          // 4 Articulated Limbs with Paws planted on ground plane
          const isFront = rng() > 0.5;
          const side = rng() > 0.5 ? 1 : -1;
          const legFrac = rng();
          const legHeight = legFrac * 5.4;

          // Muscular haunches/shoulders
          const jointBulge = Math.sin(legFrac * Math.PI) * 0.4;
          x = side * (1.65 + jointBulge);
          y = 0.6 - legHeight; // Ground level is at y = -4.8
          z = isFront ? (3.0 - legFrac * 0.6) : (-3.2 + legFrac * 0.7);
        } else {
          // Long Curving Tail
          const tL = rng() * 5.8;
          x = Math.sin(tL * 0.7) * 0.9;
          y = -0.2 + Math.sin(tL * 0.5) * 1.9;
          z = -4.3 - tL * 0.95;
        }
      } else if (morphType === 'plant') {
        // Anatomical Viridiplantae (Arabidopsis thaliana)
        if (p < 0.22) {
          // Central Erect Flowering Stem (Cauline axis)
          const sH = (rng() - 0.5) * 12.0; // -6.0 to +6.0
          const r = 0.35 + (rng() - 0.5) * 0.12;
          const a = rng() * Math.PI * 2;
          x = Math.cos(a) * r;
          y = sH;
          z = Math.sin(a) * r;
        } else if (p < 0.76) {
          // 10 Rosette Leaves spreading in spiral phyllotaxy
          const leafIdx = Math.floor(rng() * 10);
          const goldenAngle = 2.39996; // 137.5 degrees
          const leafAngle = leafIdx * goldenAngle + (rng() - 0.5) * 0.15;
          const leafLen = rng();
          const dist = 0.8 + leafLen * 7.2;
          const leafWidth = Math.sin(leafLen * Math.PI) * 2.1 * (rng() - 0.5);

          // Horizontal leaf blade with gentle downward curve
          const lx = Math.cos(leafAngle) * dist - Math.sin(leafAngle) * leafWidth;
          const lz = Math.sin(leafAngle) * dist + Math.cos(leafAngle) * leafWidth;
          const ly = -3.2 - Math.pow(leafLen, 1.8) * 1.3 + (rng() - 0.5) * 0.2;

          x = lx; y = ly; z = lz;
        } else {
          // Apical Flower Inflorescence (Raceme with buds & petals)
          const fAngle = rng() * Math.PI * 2;
          const fDist = rng() * 2.4;
          x = Math.cos(fAngle) * fDist;
          y = 5.4 + (rng() - 0.5) * 2.0;
          z = Math.sin(fAngle) * fDist;
        }
      } else {
        // Finch (Taeniopygia guttata) - Compact Songbird
        if (p < 0.52) {
          // Chubby passerine fuselage
          const u = rng() * Math.PI * 2;
          const v = (rng() - 0.5) * Math.PI;
          const r = 3.0;
          x = r * 0.95 * Math.cos(v) * Math.cos(u);
          y = r * 1.15 * Math.sin(v) - 0.2;
          z = r * 0.92 * Math.cos(v) * Math.sin(u);
        } else if (p < 0.68) {
          // Thick conical seed bill
          const bL = rng() * 2.0;
          x = (rng() - 0.5) * (1.1 - bL * 0.45);
          y = 1.9 + bL * 0.4;
          z = 2.4 + bL * 0.85;
        } else {
          // Folded Wings & Tail
          const side = rng() > 0.5 ? 1 : -1;
          const wL = rng() * 4.6;
          x = side * (2.3 + wL * 0.35);
          y = 0.3 - wL * 0.85;
          z = -0.9 - wL * 0.65;
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
    this.organismPos = this.generateOrganismPositions(this.speciesProfiles[this.currentSpecies].morph, this.particleCount);

    this.geometry = new THREE.BufferGeometry();
    this.geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(this.organismPos), 3));
    this.geometry.setAttribute('aDnaPosition', new THREE.BufferAttribute(this.dnaData.positions, 3));
    this.geometry.setAttribute('aDnaColor', new THREE.BufferAttribute(this.dnaData.colors, 3));

    // Bio-Emerald & Cyan Organism Coloring
    const orgColors = new Float32Array(this.particleCount * 3);
    const sp = this.speciesProfiles[this.currentSpecies];
    for (let i = 0; i < this.particleCount; i++) {
      const idx = i * 3;
      const alpha = Math.random();
      orgColors[idx] = THREE.MathUtils.lerp(sp.color[0], sp.secondaryColor[0], alpha);
      orgColors[idx + 1] = THREE.MathUtils.lerp(sp.color[1], sp.secondaryColor[1], alpha);
      orgColors[idx + 2] = THREE.MathUtils.lerp(sp.color[2], sp.secondaryColor[2], alpha);
    }
    this.geometry.setAttribute('aOrgColor', new THREE.BufferAttribute(orgColors, 3));

    // Random phases for organic breathing wave
    const phases = new Float32Array(this.particleCount);
    for (let i = 0; i < this.particleCount; i++) phases[i] = Math.random() * Math.PI * 2;
    this.geometry.setAttribute('aPhase', new THREE.BufferAttribute(phases, 1));
  }

  createParticleSystem() {
    const vertexShader = `
      uniform float uProgress;
      uniform float uTime;
      uniform float uGapAlarm;
      uniform float uHealing;
      uniform float uParticleBaseSize;
      
      attribute vec3 aDnaPosition;
      attribute vec3 aDnaColor;
      attribute vec3 aOrgColor;
      attribute float aPhase;
      
      varying vec3 vColor;
      varying float vAlpha;
      
      void main() {
        // High-order smoothstep for organic morphing
        float p = smoothstep(0.0, 1.0, uProgress);
        
        // Helical particle vortex wave during transition
        float vortex = sin(p * 3.14159) * 3.8;
        float angle = uTime * 2.2 + aDnaPosition.y * 0.45;
        vec3 vortexOffset = vec3(sin(angle) * vortex, sin(p * 3.14159) * 0.9, cos(angle) * vortex);
        
        vec3 mixedPos = mix(position, aDnaPosition, p) + vortexOffset;
        
        // Gentle bio-luminescent breathing wave
        mixedPos.y += sin(uTime * 2.2 + aPhase) * 0.10;
        mixedPos.x += cos(uTime * 1.8 + aPhase) * 0.07;
        
        // Color blending
        vec3 col = mix(aOrgColor, aDnaColor, p);
        
        // Broken Gap Alarm: Pulsing Neon Red at middle locus
        if (p > 0.35 && uGapAlarm > 0.05) {
          if (abs(aDnaPosition.y) < 2.8) {
            float pulse = 0.6 + 0.4 * sin(uTime * 14.0);
            col = mix(col, vec3(1.0, 0.05, 0.20), pulse * uGapAlarm);
            mixedPos.x += (sin(uTime * 28.0 + aPhase) * 0.6) * uGapAlarm;
          }
        }
        
        // Healing Stream Pulse
        if (uHealing > 0.05 && abs(aDnaPosition.y) < 3.2) {
          float healWave = sin(uTime * 16.0 - aDnaPosition.y * 3.2);
          col = mix(col, vec3(0.06, 1.0, 0.85), clamp(healWave * uHealing, 0.0, 1.0));
        }

        vColor = col;
        vAlpha = 0.88 + 0.12 * sin(uTime * 3.0 + aPhase);

        vec4 mvPosition = modelViewMatrix * vec4(mixedPos, 1.0);
        // Robust point size calculation with distance attenuation and clamping
        float pSize = (uParticleBaseSize / -mvPosition.z) * (1.0 + 0.25 * sin(uTime * 3.0 + aPhase));
        gl_PointSize = clamp(pSize, 2.5, 24.0);
        gl_Position = projectionMatrix * mvPosition;
      }
    `;

    // Self-contained mathematical circular shader with white-hot nucleus and halo
    const fragmentShader = `
      varying vec3 vColor;
      varying float vAlpha;
      
      void main() {
        vec2 coord = gl_PointCoord - vec2(0.5);
        float distSq = dot(coord, coord);
        if (distSq > 0.25) discard;
        
        float dist = sqrt(distSq) * 2.0; // 0.0 center -> 1.0 edge
        
        // Three-tier bioluminescent glow physics:
        float halo = pow(1.0 - dist, 1.8);
        float radiance = pow(1.0 - dist, 4.0);
        float whiteCore = pow(1.0 - dist, 10.0) * 1.6;
        
        vec3 finalColor = vColor * (halo * 0.75 + radiance * 1.15) + vec3(1.0) * whiteCore;
        float alpha = vAlpha * halo;
        
        gl_FragColor = vec4(finalColor, alpha);
      }
    `;

    this.material = new THREE.ShaderMaterial({
      vertexShader: vertexShader,
      fragmentShader: fragmentShader,
      uniforms: {
        uProgress: { value: 0.0 },
        uTime: { value: 0.0 },
        uGapAlarm: { value: 0.0 },
        uHealing: { value: 0.0 },
        uParticleBaseSize: { value: this.particleBaseSize }
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
    const sp = this.speciesProfiles[speciesId];

    // Recalculate organism positions
    const newPositions = this.generateOrganismPositions(sp.morph, this.particleCount);
    const posAttr = this.geometry.attributes.position;
    posAttr.array.set(newPositions);
    posAttr.needsUpdate = true;

    // Recalculate vibrant colors
    const orgColorAttr = this.geometry.attributes.aOrgColor;
    for (let i = 0; i < this.particleCount; i++) {
      const idx = i * 3;
      const alpha = Math.random();
      orgColorAttr.array[idx] = THREE.MathUtils.lerp(sp.color[0], sp.secondaryColor[0], alpha);
      orgColorAttr.array[idx + 1] = THREE.MathUtils.lerp(sp.color[1], sp.secondaryColor[1], alpha);
      orgColorAttr.array[idx + 2] = THREE.MathUtils.lerp(sp.color[2], sp.secondaryColor[2], alpha);
    }
    orgColorAttr.needsUpdate = true;
  }

  setHybridMorph(speciesA, speciesB, ratioB) {
    const spA = this.speciesProfiles[speciesA] || this.speciesProfiles.c_anna;
    const spB = this.speciesProfiles[speciesB] || this.speciesProfiles.mammal;

    const morphA = this.generateOrganismPositions(spA.morph, this.particleCount);
    const morphB = this.generateOrganismPositions(spB.morph, this.particleCount);

    const posAttr = this.geometry.attributes.position;
    for (let i = 0; i < this.particleCount * 3; i++) {
      posAttr.array[i] = THREE.MathUtils.lerp(morphA[i], morphB[i], ratioB);
    }
    posAttr.needsUpdate = true;

    // Blend Colors
    const orgColorAttr = this.geometry.attributes.aOrgColor;
    for (let i = 0; i < this.particleCount; i++) {
      const idx = i * 3;
      orgColorAttr.array[idx] = THREE.MathUtils.lerp(spA.color[0], spB.color[0], ratioB);
      orgColorAttr.array[idx + 1] = THREE.MathUtils.lerp(spA.color[1], spB.color[1], ratioB);
      orgColorAttr.array[idx + 2] = THREE.MathUtils.lerp(spA.color[2], spB.color[2], ratioB);
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

  setParticleBaseSize(size) {
    this.particleBaseSize = size;
    if (this.material && this.material.uniforms.uParticleBaseSize) {
      this.material.uniforms.uParticleBaseSize.value = size;
    }
  }

  setCameraPreset(mode) {
    if (!this.camera || !this.controls) return;
    this.controls.autoRotate = false;
    if (mode === 'isometric') {
      this.camera.position.set(18, 9, 24);
    } else if (mode === 'front') {
      this.camera.position.set(0, 3, 30);
    } else if (mode === 'side') {
      this.camera.position.set(30, 2, 0);
    } else if (mode === 'top') {
      this.camera.position.set(0, 32, 0.1);
    }
    this.controls.target.set(0, 0, 0);
    this.controls.update();
  }

  toggleAutoRotate() {
    if (this.controls) {
      this.autoRotate = !this.controls.autoRotate;
      this.controls.autoRotate = this.autoRotate;
    }
    return this.autoRotate;
  }

  bindEvents() {
    window.addEventListener('resize', () => this.handleResize());
  }

  animate() {
    requestAnimationFrame(() => this.animate());

    const delta = this.clock.getDelta();
    const elapsedTime = this.clock.getElapsedTime();

    // Smooth lerp to target morph
    this.morphProgress = THREE.MathUtils.lerp(this.morphProgress, this.targetMorphProgress, delta * 4.2);
    this.material.uniforms.uProgress.value = this.morphProgress;
    this.material.uniforms.uTime.value = elapsedTime;

    // Smooth gap alarm & healing transitions
    const targetGap = this.isGapActive ? 1.0 : 0.0;
    this.material.uniforms.uGapAlarm.value = THREE.MathUtils.lerp(this.material.uniforms.uGapAlarm.value, targetGap, delta * 6.0);

    const targetHealing = this.isHealing ? 1.0 : 0.0;
    this.material.uniforms.uHealing.value = THREE.MathUtils.lerp(this.material.uniforms.uHealing.value, targetHealing, delta * 5.0);

    if (this.controls) this.controls.update();
    this.renderer.render(this.scene, this.camera);
  }
}

/**
 * Lightweight Autonomous 3D Floating DNA Helix for Landing Page Hero
 * Enhanced with interactive mouse drag rotation, self-contained GLSL shader, and ResizeObserver.
 */
class HeroDnaBackground {
  constructor(canvasContainerId) {
    this.container = document.getElementById(canvasContainerId);
    if (!this.container) {
      console.warn(`HeroDnaBackground: Container #${canvasContainerId} not found.`);
      return;
    }

    const width = this.container.clientWidth || 480;
    const height = this.container.clientHeight || 320;

    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(40, width / height, 0.1, 100);
    this.camera.position.set(0, 0, 22);

    this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    this.renderer.setSize(width, height);
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    this.renderer.domElement.style.width = '100%';
    this.renderer.domElement.style.height = '100%';
    this.renderer.domElement.style.display = 'block';

    // Clear previous canvas
    while (this.container.firstChild && this.container.firstChild.tagName === 'CANVAS') {
      this.container.removeChild(this.container.firstChild);
    }
    this.container.appendChild(this.renderer.domElement);

    this.buildHelix();
    this.clock = new THREE.Clock();
    this.bindInteractions();
    this.animate();

    if (window.ResizeObserver) {
      this.resizeObserver = new ResizeObserver(() => this.handleResize());
      this.resizeObserver.observe(this.container);
    }
  }

  handleResize() {
    if (!this.container || !this.renderer || !this.camera) return;
    const w = this.container.clientWidth;
    const h = this.container.clientHeight;
    if (w === 0 || h === 0) return;
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
    this.renderer.setSize(w, h);
  }

  buildHelix() {
    const particleCount = 3600;
    const geometry = new THREE.BufferGeometry();
    const positions = new Float32Array(particleCount * 3);
    const colors = new Float32Array(particleCount * 3);
    const phases = new Float32Array(particleCount);

    const heightSpan = 18.0;
    const radius = 3.6;
    const turns = 4.2;

    for (let i = 0; i < particleCount; i++) {
      const frac = i / particleCount;
      const theta = frac * Math.PI * 2 * turns;
      const strand = (i % 2 === 0) ? 1.0 : -1.0;
      const y = (frac - 0.5) * heightSpan;

      const idx = i * 3;
      phases[i] = Math.random() * Math.PI * 2;

      if (i % 3 === 0) {
        // Cross rung particle
        const alpha = (i % 32) / 32.0;
        const p1 = new THREE.Vector3(Math.cos(theta) * radius, y, Math.sin(theta) * radius);
        const p2 = new THREE.Vector3(Math.cos(theta + Math.PI) * radius, y, Math.sin(theta + Math.PI) * radius);
        const pt = new THREE.Vector3().lerpVectors(p1, p2, alpha);
        positions[idx] = pt.x;
        positions[idx + 1] = pt.y;
        positions[idx + 2] = pt.z;

        // Base pairs
        if (alpha < 0.5) {
          colors[idx] = 0.05; colors[idx + 1] = 0.85; colors[idx + 2] = 1.0; // Cyan
        } else {
          colors[idx] = 0.06; colors[idx + 1] = 0.95; colors[idx + 2] = 0.65; // Emerald
        }
      } else {
        // Strand particle
        positions[idx] = Math.cos(theta) * radius * strand;
        positions[idx + 1] = y;
        positions[idx + 2] = Math.sin(theta) * radius * strand;

        if (strand > 0) {
          colors[idx] = 0.06; colors[idx + 1] = 0.95; colors[idx + 2] = 0.65; // Emerald strand
        } else {
          colors[idx] = 0.05; colors[idx + 1] = 0.85; colors[idx + 2] = 1.0; // Cyan strand
        }
      }
    }

    geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
    geometry.setAttribute('aPhase', new THREE.BufferAttribute(phases, 1));

    const vertexShader = `
      attribute vec3 color;
      attribute float aPhase;
      varying vec3 vColor;
      uniform float uTime;
      
      void main() {
        vColor = color;
        vec3 pos = position;
        pos.y += sin(uTime * 2.0 + aPhase) * 0.08;
        vec4 mvPosition = modelViewMatrix * vec4(pos, 1.0);
        gl_PointSize = clamp(60.0 / -mvPosition.z, 2.5, 14.0);
        gl_Position = projectionMatrix * mvPosition;
      }
    `;

    const fragmentShader = `
      varying vec3 vColor;
      
      void main() {
        vec2 coord = gl_PointCoord - vec2(0.5);
        float distSq = dot(coord, coord);
        if (distSq > 0.25) discard;
        
        float dist = sqrt(distSq) * 2.0;
        float halo = pow(1.0 - dist, 1.6);
        float core = pow(1.0 - dist, 8.0) * 1.5;
        
        vec3 col = vColor * halo + vec3(1.0) * core;
        gl_FragColor = vec4(col, halo * 0.92);
      }
    `;

    this.helixMaterial = new THREE.ShaderMaterial({
      vertexShader: vertexShader,
      fragmentShader: fragmentShader,
      uniforms: {
        uTime: { value: 0.0 }
      },
      transparent: true,
      depthWrite: false,
      blending: THREE.AdditiveBlending
    });

    this.helixMesh = new THREE.Points(geometry, this.helixMaterial);
    this.scene.add(this.helixMesh);
  }

  bindInteractions() {
    this.isDragging = false;
    this.previousMousePosition = { x: 0, y: 0 };
    this.dragRotationVelocity = { x: 0, y: 0 };

    const el = this.renderer.domElement;

    el.addEventListener('mousedown', (e) => {
      this.isDragging = true;
      this.previousMousePosition = { x: e.clientX, y: e.clientY };
    });

    window.addEventListener('mouseup', () => {
      this.isDragging = false;
    });

    window.addEventListener('mousemove', (e) => {
      if (!this.isDragging || !this.helixMesh) return;
      const deltaX = e.clientX - this.previousMousePosition.x;
      const deltaY = e.clientY - this.previousMousePosition.y;

      this.dragRotationVelocity.y = deltaX * 0.006;
      this.dragRotationVelocity.x = deltaY * 0.006;

      this.helixMesh.rotation.y += this.dragRotationVelocity.y;
      this.helixMesh.rotation.x += this.dragRotationVelocity.x;

      this.previousMousePosition = { x: e.clientX, y: e.clientY };
    });

    // Touch support for mobile
    el.addEventListener('touchstart', (e) => {
      if (e.touches.length === 1) {
        this.isDragging = true;
        this.previousMousePosition = { x: e.touches[0].clientX, y: e.touches[0].clientY };
      }
    }, { passive: true });

    window.addEventListener('touchend', () => {
      this.isDragging = false;
    });

    window.addEventListener('touchmove', (e) => {
      if (!this.isDragging || !this.helixMesh || e.touches.length !== 1) return;
      const deltaX = e.touches[0].clientX - this.previousMousePosition.x;
      const deltaY = e.touches[0].clientY - this.previousMousePosition.y;

      this.helixMesh.rotation.y += deltaX * 0.008;
      this.helixMesh.rotation.x += deltaY * 0.008;

      this.previousMousePosition = { x: e.touches[0].clientX, y: e.touches[0].clientY };
    }, { passive: true });
  }

  animate() {
    requestAnimationFrame(() => this.animate());
    const elapsedTime = this.clock.getElapsedTime();

    if (this.helixMaterial && this.helixMaterial.uniforms.uTime) {
      this.helixMaterial.uniforms.uTime.value = elapsedTime;
    }

    if (this.helixMesh) {
      if (!this.isDragging) {
        this.helixMesh.rotation.y += 0.006;
        this.helixMesh.rotation.z = Math.sin(elapsedTime * 0.3) * 0.12;
      }
      this.helixMesh.position.y = Math.sin(elapsedTime * 0.8) * 0.4;
    }

    this.renderer.render(this.scene, this.camera);
  }
}

// Global hooks
window.PhyloMorphEngine = PhyloMorphEngine;
window.HeroDnaBackground = HeroDnaBackground;
