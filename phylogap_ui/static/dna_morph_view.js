/**
 * PhyloMorph 3D Studio - High-Luminescence Particle Morphing Engine.
 * Features:
 * 1. Ultra-Vivid Bio-Luminescent Particle Shader (White-hot core + Emerald/Cyan aura).
 * 2. 10,000 Particles with GPU-Accelerated Vortex Morphing.
 * 3. Sculpted 3D Models:
 *    - Avian (Calypte anna - Hummingbird with aerodynamic body, needle bill, swept wings)
 *    - Mammalia (Quadruped Vertebrate with spine, 4 limbs, head, paws, tail)
 *    - Viridiplantae (Arabidopsis thaliana - Rosette leaves with venation and flower stalk)
 *    - Finch (Taeniopygia guttata - Compact songbird)
 *    - DNA Double Helix (Parametric dual strands + 48 nucleotide rungs in A, T, G, C)
 * 4. HeroDnaBackground: Lightweight autonomous floating 3D DNA helix for the landing page hero.
 * 5. Interactive Gap Injection & Healing simulation with streaming repair particles.
 *
 * Author: PhyloGap Research & Product Team
 * Date: 2026-10-28
 */

class PhyloMorphEngine {
  constructor(containerId) {
    this.container = document.getElementById(containerId);
    if (!this.container) return;

    this.particleCount = 9600;
    this.currentSpecies = 'c_anna';
    this.morphProgress = 0.0; // 0.0 = Organism Phenotype, 1.0 = DNA Helix
    this.targetMorphProgress = 0.0;
    this.isGapActive = false;
    this.isHealing = false;

    this.speciesProfiles = {
      c_anna: {
        name: "Calypte anna",
        common: "Anna's Hummingbird (Chim Ruồi)",
        morph: "hummingbird",
        color: [0.0, 1.0, 0.65], // Bio Emerald
        secondaryColor: [0.1, 0.85, 1.0]
      },
      mammal: {
        name: "Mammalia (Thú)",
        common: "Vertebrate Model (Động Vật Có Vú)",
        morph: "mammal",
        color: [0.2, 0.9, 0.5], // Electric Mint
        secondaryColor: [0.0, 0.7, 0.9]
      },
      plant: {
        name: "Arabidopsis thaliana",
        common: "Thale Cress (Thực Vật)",
        morph: "plant",
        color: [0.3, 1.0, 0.4], // Neon Botanical Green
        secondaryColor: [0.1, 0.9, 0.7]
      },
      t_guttata: {
        name: "Taeniopygia guttata",
        common: "Zebra Finch (Chim Sẻ)",
        morph: "finch",
        color: [0.05, 0.95, 0.8], // Cyan-Emerald
        secondaryColor: [0.2, 1.0, 0.5]
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
    this.scene.fog = new THREE.FogExp2(0x030303, 0.012);

    this.camera = new THREE.PerspectiveCamera(45, width / height, 0.1, 1000);
    this.camera.position.set(0, 4, 30);

    this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true, powerPreference: "high-performance" });
    this.renderer.setSize(width, height);
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    this.container.appendChild(this.renderer.domElement);

    if (window.THREE && THREE.OrbitControls) {
      this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
      this.controls.enableDamping = true;
      this.controls.dampingFactor = 0.05;
      this.controls.autoRotate = true;
      this.controls.autoRotateSpeed = 0.9;
      this.controls.maxDistance = 65;
      this.controls.minDistance = 6;
    }

    this.clock = new THREE.Clock();
  }

  createHighIntensityGlowTexture() {
    const canvas = document.createElement('canvas');
    canvas.width = 128;
    canvas.height = 128;
    const ctx = canvas.getContext('2d');

    const grad = ctx.createRadialGradient(64, 64, 0, 64, 64, 64);
    grad.addColorStop(0, 'rgba(255, 255, 255, 1.0)');     // Brilliant white-hot nucleus
    grad.addColorStop(0.18, 'rgba(255, 255, 255, 0.95)');
    grad.addColorStop(0.45, 'rgba(16, 185, 129, 0.55)');   // Bio-emerald halo
    grad.addColorStop(0.75, 'rgba(6, 182, 212, 0.20)');    // Cyan outer aura
    grad.addColorStop(1, 'rgba(0, 0, 0, 0)');
    ctx.fillStyle = grad;
    ctx.fillRect(0, 0, 128, 128);

    const texture = new THREE.CanvasTexture(canvas);
    texture.needsUpdate = true;
    return texture;
  }

  generateDnaPositions(count) {
    const pos = new Float32Array(count * 3);
    const colors = new Float32Array(count * 3);
    const heightSpan = 24.0;
    const radius = 4.2;
    const turns = 5.0;

    // Vibrant Base Colors:
    // A: Crimson/Rose, T: Gold/Amber, G: Emerald Green, C: Neon Cyan
    const basePalette = [
      [1.0, 0.2, 0.45],  // A (Rose-Crimson)
      [1.0, 0.85, 0.15], // T (Gold-Amber)
      [0.0, 1.0, 0.55],  // G (Emerald Green)
      [0.0, 0.85, 1.0]   // C (Cyber Cyan)
    ];

    const backboneCount = Math.floor(count * 0.48);
    const rungsCount = count - backboneCount;

    // 1. Double Helix Backbones (Thick, luminous double strands)
    for (let i = 0; i < backboneCount; i++) {
      const frac = (i / backboneCount);
      const theta = frac * Math.PI * 2 * turns;
      const strand = (i % 2 === 0) ? 1.0 : -1.0;
      const y = (frac - 0.5) * heightSpan;

      // Slight thickness jitter
      const rJitter = radius + (Math.sin(i * 12.0) * 0.15);
      const x = Math.cos(theta) * rJitter * strand;
      const z = Math.sin(theta) * rJitter * strand;

      const idx = i * 3;
      pos[idx] = x;
      pos[idx + 1] = y;
      pos[idx + 2] = z;

      // Bio-emerald & Cyber-cyan backbones
      if (strand > 0) {
        colors[idx] = 0.0; colors[idx + 1] = 1.0; colors[idx + 2] = 0.7; // Emerald
      } else {
        colors[idx] = 0.0; colors[idx + 1] = 0.85; colors[idx + 2] = 1.0; // Cyan
      }
    }

    // 2. Base-Pair Cross Rungs
    const offset = backboneCount * 3;
    const baseSteps = 48; // 48 paired rungs along the helix
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
        const pt = new THREE.Vector3().lerpVectors(p1, p2, alpha);

        const idx = offset + rungIdx * 3;
        pos[idx] = pt.x + (Math.random() - 0.5) * 0.1;
        pos[idx + 1] = pt.y + (Math.random() - 0.5) * 0.1;
        pos[idx + 2] = pt.z + (Math.random() - 0.5) * 0.1;

        const c = alpha < 0.5 ? colLeft : colRight;
        colors[idx] = c[0];
        colors[idx + 1] = c[1];
        colors[idx + 2] = c[2];
      }
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
        // Anatomical Hummingbird (Calypte anna):
        // Slender aerodynamic fuselage, spherical head, needle-thin straight beak, swept wings, fanned tail
        if (p < 0.30) {
          // Torpedo Torso & Breast
          const u = rng() * Math.PI * 2;
          const v = (rng() - 0.5) * Math.PI;
          const r = 2.4;
          x = r * 0.7 * Math.cos(v) * Math.cos(u);
          y = (r * 1.6 * Math.sin(v)) - 0.6;
          z = r * 0.75 * Math.cos(v) * Math.sin(u);
        } else if (p < 0.44) {
          // Cranial Head Dome + Needle Beak
          if (rng() < 0.4) {
            // Head sphere
            const u = rng() * Math.PI * 2;
            const v = (rng() - 0.5) * Math.PI;
            x = 0.9 * Math.cos(v) * Math.cos(u);
            y = 2.2 + 0.9 * Math.sin(v);
            z = 1.0 + 0.9 * Math.cos(v) * Math.sin(u);
          } else {
            // Needle-thin long bill
            const bL = rng() * 4.8;
            x = (rng() - 0.5) * 0.12;
            y = 2.1 + bL * 0.45;
            z = 1.8 + bL * 0.9;
          }
        } else if (p < 0.82) {
          // Swept-back articulated wings with primary feathers
          const side = rng() > 0.5 ? 1 : -1;
          const span = rng(); // 0 to 1
          const wingL = span * 8.5;
          const chord = (1.0 - span * 0.65) * 2.2 * (rng() - 0.5);

          x = side * (1.2 + wingL);
          y = 0.6 - (span * 2.2) + chord * 0.3;
          z = -0.5 - (span * 4.5) + chord;
        } else {
          // Split fork tail
          const side = rng() > 0.5 ? 1 : -1;
          const tL = rng() * 4.5;
          x = side * (tL * 0.35) + (rng() - 0.5) * 0.3;
          y = -2.6 - tL * 0.9;
          z = -1.6 - tL * 0.6;
        }
      } else if (morphType === 'mammal') {
        // Anatomical Mammalia (Vertebrate Quadruped):
        // Arched spine, thoracic ribcage, 4 articulated legs, skull with ears, tail
        if (p < 0.32) {
          // Arched Spine & Ribcage
          const sT = (rng() - 0.5) * 8.0; // -4 to +4 along z
          const arch = Math.cos((sT / 4.0) * (Math.PI / 2.5)) * 1.5;
          const ribR = 2.0 * Math.sin(rng() * Math.PI);
          const ribA = rng() * Math.PI * 2;
          x = Math.cos(ribA) * ribR * 0.8;
          y = arch + Math.sin(ribA) * ribR * 0.85;
          z = sT;
        } else if (p < 0.48) {
          // Cranial Head, Muzzle, and Ears
          const hP = rng();
          if (hP < 0.5) {
            // Skull
            const u = rng() * Math.PI * 2;
            const v = (rng() - 0.5) * Math.PI;
            x = 1.1 * Math.cos(v) * Math.cos(u);
            y = 1.8 + 1.1 * Math.sin(v);
            z = 4.2 + 1.1 * Math.cos(v) * Math.sin(u);
          } else if (hP < 0.8) {
            // Muzzle
            x = (rng() - 0.5) * 0.8;
            y = 1.4 + (rng() - 0.5) * 0.6;
            z = 5.2 + rng() * 1.2;
          } else {
            // Triangular Ears
            const side = rng() > 0.5 ? 1 : -1;
            const eH = rng() * 1.2;
            x = side * (0.8 + eH * 0.3);
            y = 2.6 + eH;
            z = 4.0 - eH * 0.2;
          }
        } else if (p < 0.84) {
          // 4 Articulated Limbs
          const isFront = rng() > 0.5;
          const side = rng() > 0.5 ? 1 : -1;
          const legFrac = rng();
          const legH = legFrac * 5.2;

          x = side * (1.6 + Math.sin(legFrac * Math.PI) * 0.3);
          y = 0.5 - legH;
          z = isFront ? (2.8 - legFrac * 0.6) : (-3.0 + legFrac * 0.8);
        } else {
          // Graceful Curved Tail
          const tL = rng() * 5.5;
          x = Math.sin(tL * 0.8) * 0.8;
          y = -0.2 + Math.sin(tL * 0.6) * 1.8;
          z = -4.0 - tL * 0.9;
        }
      } else if (morphType === 'plant') {
        // Anatomical Viridiplantae (Arabidopsis thaliana):
        // Central botanical stem, 8 radiating elliptical rosette leaves with central veins, apical inflorescence
        if (p < 0.20) {
          // Central Vertical Stem
          const sH = (rng() - 0.5) * 12.0;
          const r = 0.35 + (rng() - 0.5) * 0.1;
          const a = rng() * Math.PI * 2;
          x = Math.cos(a) * r;
          y = sH;
          z = Math.sin(a) * r;
        } else if (p < 0.75) {
          // 8 Rosette Leaves radiating symmetrically
          const leafIdx = Math.floor(rng() * 8);
          const leafAngle = (leafIdx / 8.0) * Math.PI * 2 + (rng() - 0.5) * 0.2;
          const leafLen = rng(); // 0 to 1 along leaf axis
          const dist = 1.0 + leafLen * 6.5;
          const leafWidth = Math.sin(leafLen * Math.PI) * 1.8 * (rng() - 0.5);

          // Horizontal leaf blade with gentle downward arch
          const lx = Math.cos(leafAngle) * dist - Math.sin(leafAngle) * leafWidth;
          const lz = Math.sin(leafAngle) * dist + Math.cos(leafAngle) * leafWidth;
          const ly = -2.5 - Math.pow(leafLen, 1.8) * 1.2 + (rng() - 0.5) * 0.2;

          x = lx;
          y = ly;
          z = lz;
        } else {
          // Apical Flower Inflorescence (Top cluster of petals and buds)
          const fAngle = rng() * Math.PI * 2;
          const fDist = rng() * 2.2;
          x = Math.cos(fAngle) * fDist;
          y = 5.2 + (rng() - 0.5) * 1.8;
          z = Math.sin(fAngle) * fDist;
        }
      } else {
        // Finch (Taeniopygia guttata):
        // Plump, round compact songbird
        if (p < 0.55) {
          const u = rng() * Math.PI * 2;
          const v = (rng() - 0.5) * Math.PI;
          const r = 2.9;
          x = r * 0.95 * Math.cos(v) * Math.cos(u);
          y = r * 1.1 * Math.sin(v) - 0.2;
          z = r * 0.9 * Math.cos(v) * Math.sin(u);
        } else if (p < 0.70) {
          // Short conical seed beak
          const bL = rng() * 1.8;
          x = (rng() - 0.5) * (1.1 - bL * 0.5);
          y = 1.8 + bL * 0.4;
          z = 2.2 + bL * 0.8;
        } else {
          // Folded Perching Wings & Short Tail
          const side = rng() > 0.5 ? 1 : -1;
          const wL = rng() * 4.2;
          x = side * (2.2 + wL * 0.35);
          y = 0.2 - wL * 0.8;
          z = -0.8 - wL * 0.6;
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

    // Random phases for organic particle floating
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
        float vortex = sin(p * 3.14159) * 3.4;
        float angle = uTime * 2.0 + aDnaPosition.y * 0.45;
        vec3 vortexOffset = vec3(sin(angle) * vortex, sin(p * 3.14159) * 0.8, cos(angle) * vortex);
        
        vec3 mixedPos = mix(position, aDnaPosition, p) + vortexOffset;
        
        // Gentle bio-luminescent breathing wave
        mixedPos.y += sin(uTime * 2.2 + aPhase) * 0.09;
        mixedPos.x += cos(uTime * 1.8 + aPhase) * 0.06;
        
        // Color blending
        vec3 col = mix(aOrgColor, aDnaColor, p);
        
        // Broken Gap Alarm: Pulsing Neon Red at middle locus
        if (p > 0.35 && uGapAlarm > 0.05) {
          if (abs(aDnaPosition.y) < 2.8) {
            float pulse = 0.6 + 0.4 * sin(uTime * 12.0);
            col = mix(col, vec3(1.0, 0.05, 0.15), pulse * uGapAlarm);
            mixedPos.x += (sin(uTime * 25.0 + aPhase) * 0.5) * uGapAlarm;
          }
        }
        
        // Healing Stream: Emerald/Cyan restoration wave
        if (uHealing > 0.05) {
          if (abs(aDnaPosition.y) < 3.2) {
            col = mix(col, vec3(0.1, 1.0, 0.65), uHealing);
          }
        }
        
        vColor = col;
        vAlpha = 0.95;
        
        vec4 mvPosition = modelViewMatrix * vec4(mixedPos, 1.0);
        gl_Position = projectionMatrix * mvPosition;
        
        // High-impact point sizing with perspective scaling
        gl_PointSize = (32.0 / -mvPosition.z) * (1.0 + sin(aPhase + uTime * 2.0) * 0.25);
      }
    `;

    const fragmentShader = `
      uniform sampler2D uTexture;
      varying vec3 vColor;
      varying float vAlpha;
      
      void main() {
        vec4 tex = texture2D(uTexture, gl_PointCoord);
        if (tex.a < 0.02) discard;
        
        // White-hot core + intense bio-luminescent aura
        float dist = length(gl_PointCoord - vec2(0.5));
        float glow = pow(max(0.0, 1.0 - dist * 2.0), 1.6);
        float core = pow(max(0.0, 1.0 - dist * 2.0), 4.5) * 1.6;
        
        vec3 finalColor = vColor * glow + vec3(1.0) * core;
        gl_FragColor = vec4(finalColor, vAlpha * glow);
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
        uTexture: { value: this.createHighIntensityGlowTexture() }
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
    posAttr.copyArray(newPositions);
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

/**
 * Lightweight Autonomous 3D Floating DNA Helix for Landing Page Hero
 */
class HeroDnaBackground {
  constructor(canvasContainerId) {
    this.container = document.getElementById(canvasContainerId);
    if (!this.container) return;

    const width = this.container.clientWidth || 480;
    const height = this.container.clientHeight || 480;

    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(40, width / height, 0.1, 100);
    this.camera.position.set(0, 0, 22);

    this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    this.renderer.setSize(width, height);
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    this.container.appendChild(this.renderer.domElement);

    this.buildHelix();
    this.clock = new THREE.Clock();
    this.animate();

    window.addEventListener('resize', () => {
      if (!this.container) return;
      const w = this.container.clientWidth;
      const h = this.container.clientHeight;
      if (w === 0 || h === 0) return;
      this.camera.aspect = w / h;
      this.camera.updateProjectionMatrix();
      this.renderer.setSize(w, h);
    });
  }

  buildHelix() {
    const particleCount = 2800;
    const geometry = new THREE.BufferGeometry();
    const positions = new Float32Array(particleCount * 3);
    const colors = new Float32Array(particleCount * 3);

    const heightSpan = 18.0;
    const radius = 3.6;
    const turns = 4.0;

    for (let i = 0; i < particleCount; i++) {
      const frac = i / particleCount;
      const theta = frac * Math.PI * 2 * turns;
      const strand = (i % 2 === 0) ? 1.0 : -1.0;
      const y = (frac - 0.5) * heightSpan;

      const idx = i * 3;
      if (i % 3 === 0) {
        // Cross rung particle
        const alpha = (i % 30) / 30.0;
        const p1 = new THREE.Vector3(Math.cos(theta) * radius, y, Math.sin(theta) * radius);
        const p2 = new THREE.Vector3(Math.cos(theta + Math.PI) * radius, y, Math.sin(theta + Math.PI) * radius);
        const pt = new THREE.Vector3().lerpVectors(p1, p2, alpha);
        positions[idx] = pt.x;
        positions[idx + 1] = pt.y;
        positions[idx + 2] = pt.z;

        // Base pairs
        colors[idx] = 0.06; colors[idx + 1] = 0.85; colors[idx + 2] = 0.95; // Cyan
      } else {
        // Strand particle
        positions[idx] = Math.cos(theta) * radius * strand;
        positions[idx + 1] = y;
        positions[idx + 2] = Math.sin(theta) * radius * strand;

        // Emerald strand
        colors[idx] = 0.05; colors[idx + 1] = 0.95; colors[idx + 2] = 0.55;
      }
    }

    geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));

    const material = new THREE.PointsMaterial({
      size: 0.16,
      vertexColors: true,
      transparent: true,
      opacity: 0.9,
      blending: THREE.AdditiveBlending
    });

    this.helixMesh = new THREE.Points(geometry, material);
    this.scene.add(this.helixMesh);
  }

  animate() {
    requestAnimationFrame(() => this.animate());
    const elapsedTime = this.clock.getElapsedTime();

    if (this.helixMesh) {
      this.helixMesh.rotation.y = elapsedTime * 0.45;
      this.helixMesh.rotation.z = Math.sin(elapsedTime * 0.3) * 0.12;
      this.helixMesh.position.y = Math.sin(elapsedTime * 0.8) * 0.4;
    }

    this.renderer.render(this.scene, this.camera);
  }
}

// Global hooks
window.PhyloMorphEngine = PhyloMorphEngine;
window.HeroDnaBackground = HeroDnaBackground;
