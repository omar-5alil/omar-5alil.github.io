// public/sim-worker.js
/* eslint-disable no-restricted-globals */
/**
 * Atmospheric entry worker (indefinite until impact or safety cap).
 * Returns typed array buffers for: x distance, lat, lon, altitude, speed, vertical velocity, times.
 * Only appends an impact sample (alt=0) when trajectory actually crosses ground.
 * Adds meta diagnostics: impacted, steps, impactTime, lastAlt, lastTime.
 *
 * Message in: { cmd:'run', params:{ lat0, lon0, bearing, diameter, density, v0, angleDeg, alt0?, dt? } }
 * Message out: { type:'simResult', xs,lats,lons,alts,vs,wzs,times, crater_km,mw,tsunami_m, meta }
 */

/**
 * @typedef {{ lat0:number, lon0:number, bearing:number, diameter:number, density:number, v0:number, angleDeg:number, alt0?:number, dt?:number }} SimParams
 */

self.addEventListener('message', (ev) => {
  const data = ev.data;
  if (!data || data.cmd !== 'run') return;
  const p = /** @type {SimParams} */ (data.params);
  const out = runSimulation(p);
  self.postMessage({
    type: 'simResult',
    xs: out.xs,
    lats: out.lats,
    lons: out.lons,
    alts: out.alts,
    vs: out.vs,
    wzs: out.wzs,
    times: out.times,
    ms: out.ms,
    crater_km: out.crater_km,
    mw: out.mw,
    tsunami_m: out.tsunami_m,
    meta: out.meta
  }, [out.xs, out.lats, out.lons, out.alts, out.vs, out.wzs, out.times, out.ms]);
});

function toRad(d) { return d * Math.PI / 180; }
function toDeg(r) { return r * 180 / Math.PI; }

function destinationLatLon(lat1, lon1, bearingDeg, distance_m) {
  const R = 6371000;
  const φ1 = toRad(lat1);
  const λ1 = toRad(lon1);
  const θ = toRad(bearingDeg);
  const δ = distance_m / R;
  const sinφ2 = Math.sin(φ1) * Math.cos(δ) + Math.cos(φ1) * Math.sin(δ) * Math.cos(θ);
  const φ2 = Math.asin(Math.max(-1, Math.min(1, sinφ2)));
  const y = Math.sin(θ) * Math.sin(δ) * Math.cos(φ1);
  const x = Math.cos(δ) - Math.sin(φ1) * Math.sin(φ2);
  const λ2 = λ1 + Math.atan2(y, x);
  return [toDeg(φ2), ((toDeg(λ2) + 540) % 360) - 180];
}

function rhoAtmosphere(z) {
  if (z > 86000) { // very thin
    return 1.225 * Math.exp(-z / 40000) * 1e-6;
  } else if (z > 50000) { // mesosphere
    const rho50 = 1.225 * Math.exp(-50000 / 8500);
    return rho50 * Math.exp(-(z - 50000) / 7000);
  } else { // troposphere/stratosphere simplified exponential
    const H = 8500;
    return 1.225 * Math.exp(-z / H);
  }
}

function derivsVec(state, p) {
  const vx = state.vx, wz = state.wz, m = state.m;
  const v = Math.sqrt(vx*vx + wz*wz) || 1e-6;
  const A = Math.PI * (p.radius * p.radius);
  const rho = rhoAtmosphere(Math.max(0, state.z));

  // Drag (unchanged)
  let Cd = p.Cd;
  if (v > 15000) Cd = 0.3; else if (v > 8000) Cd = 0.5; else Cd = 1.0;
  const dragFactor = 0.5 * Cd * rho * A * v * v / Math.max(1e-9, m);
  const ax = -dragFactor * (vx / v);
  const az = -dragFactor * (wz / v) - p.g;

  // --- Moderated ablation model ---
  // Original: dmdt_raw = -(p.Ch * rho * A * v^2) / (2 * p.Q)
  // We apply multipliers so near-ground mass loss is less abrupt.
  const HEAT = {
    vCut: 4000,      // below this speed ablation efficiency ~0
    vRef: 12000,     // above this speed velocity factor ~1
    altShield: 10000, // m, altitude scale controlling near-ground taper
    dynP0: 2e5,      // Pa, soft reference for dynamic pressure ramp
    maxRho: 0.5      // kg/m^3 cap for ablation heating (dense boundary layer mitigation)
  };

  const clamp01 = (x) => x < 0 ? 0 : (x > 1 ? 1 : x);

  // Effective density (prevents explosive growth at very low altitude)
  const rhoEff = Math.min(rho, HEAT.maxRho);

  // Velocity efficiency (smooth 0→1)
  let effV = (v - HEAT.vCut) / (HEAT.vRef - HEAT.vCut);
  effV = Math.pow(clamp01(effV), 1.3); // slight bias toward higher speeds

  // Altitude shielding (→0 at ground, ~1 high up)
  const zPos = Math.max(0, state.z);
  const effAlt = zPos / (zPos + HEAT.altShield); // simple rational taper

  // Dynamic pressure ramp
  const dynP = 0.5 * rhoEff * v * v;
  const effDynP = dynP / (dynP + HEAT.dynP0);

  const heatMult = effV * effAlt * effDynP;

  const dmdt_base = - (p.Ch * rhoEff * A * v * v) / (2 * p.Q);
  const dmdt = dmdt_base * heatMult;

  return { dxdt: vx, dzdt: wz, ax, az, dmdt, speed: v };
}

function rk4Vec(state, dt, p) {
  const s1 = derivsVec(state, p);
  const st2 = { x: state.x + s1.dxdt*dt/2, z: state.z + s1.dzdt*dt/2, vx: state.vx + s1.ax*dt/2, wz: state.wz + s1.az*dt/2, m: state.m + s1.dmdt*dt/2 };
  const s2 = derivsVec(st2, p);
  const st3 = { x: state.x + s2.dxdt*dt/2, z: state.z + s2.dzdt*dt/2, vx: state.vx + s2.ax*dt/2, wz: state.wz + s2.az*dt/2, m: state.m + s2.dmdt*dt/2 };
  const s3 = derivsVec(st3, p);
  const st4 = { x: state.x + s3.dxdt*dt, z: state.z + s3.dzdt*dt, vx: state.vx + s3.ax*dt, wz: state.wz + s3.az*dt, m: state.m + s3.dmdt*dt };
  const s4 = derivsVec(st4, p);
  const vxNew = state.vx + dt * (s1.ax + 2*s2.ax + 2*s3.ax + s4.ax) / 6;
  const wzNew = state.wz + dt * (s1.az + 2*s2.az + 2*s3.az + s4.az) / 6;
  return {
    x: state.x + dt * (s1.dxdt + 2*s2.dxdt + 2*s3.dxdt + s4.dxdt) / 6,
    z: state.z + dt * (s1.dzdt + 2*s2.dzdt + 2*s3.dzdt + s4.dzdt) / 6,
    vx: vxNew,
    wz: wzNew,
    m: state.m + dt * (s1.dmdt + 2*s2.dmdt + 2*s3.dmdt + s4.dmdt) / 6,
    speed: Math.sqrt(vxNew*vxNew + wzNew*wzNew)
  };
}

function adaptiveTimeStep(state, p, dt_base) {
  const v = Math.sqrt(state.vx*state.vx + state.wz*state.wz);
  let dt_adaptive = dt_base;
  if (v > 25000) dt_adaptive = Math.min(dt_adaptive, 0.0005);
  else if (v > 20000) dt_adaptive = Math.min(dt_adaptive, 0.001);
  else if (v > 15000) dt_adaptive = Math.min(dt_adaptive, 0.002);
  else if (v > 10000) dt_adaptive = Math.min(dt_adaptive, 0.005);
  if (state.z > 80000) dt_adaptive = Math.min(dt_adaptive, 0.01);
  else if (state.z > 50000) dt_adaptive = Math.min(dt_adaptive, 0.005);
  else if (state.z > 20000) dt_adaptive = Math.min(dt_adaptive, 0.002);
  if (Math.abs(state.wz) > 0) {
    const dt_alt = 50.0 / Math.abs(state.wz); // limit altitude change per step
    dt_adaptive = Math.min(dt_adaptive, dt_alt);
  }
  return Math.max(1e-6, dt_adaptive);
}

function lerp(a,b,f){ return a + (b-a)*f; }

function runSimulation(pIn) {
  const lat0 = pIn.lat0;
  const lon0 = pIn.lon0;
  const bearing = pIn.bearing;
  const diameter = pIn.diameter;
  const density = pIn.density;
  const v0 = pIn.v0;
  const angleDeg = pIn.angleDeg;
  const startAlt = (typeof pIn.alt0 === 'number') ? pIn.alt0 : 100000;
  const dtBase = (pIn.dt && pIn.dt > 0) ? pIn.dt : 0.02;

  // Physical parameters
  const radius = diameter / 2;
  const mass0 = (4/3) * Math.PI * Math.pow(radius, 3) * density;
  // Adjusted ablation parameters: lower heating efficiency (Ch) and higher heat of ablation (Q)
  // to prevent unrealistically rapid mass loss.
  const params = { radius, Cd: 1.0, Ch: 0.005, Q: 5e7, g: 9.81 };
  const strength = 10e6; // fragmentation threshold reference

  // Initial state
  const ang = toRad(Math.abs(angleDeg));
  let vx = v0 * Math.cos(ang);
  let wz = -v0 * Math.sin(ang); // negative downward
  let x = 0;
  let z = startAlt;
  let m = mass0;
  let t = 0;

  // Dynamic storage
  const xs = [], lats = [], lons = [], alts = [], vs = [], wzs = [], times = [], ms = [];
  const HARD_STEP_LIMIT = 2_000_000; // safeguard ~2M samples
  const REPORT_INTERVAL = 5000;
  const z_adapt_threshold = 10000;
  let steps = 0; let impacted = false; let impactTime = null;

  let disintegrated = false;
  let disintegrationIndex = -1;
  let disintegrationTime = null;
  const MIN_FRACTION_FOR_CONTINUATION = 0.0001; // 0.01% of original mass
  while (steps < HARD_STEP_LIMIT) {
    if (z <= 0) { impacted = true; break; }

    const [plat, plon] = destinationLatLon(lat0, lon0, bearing, x);
    xs.push(x); lats.push(plat); lons.push(plon); alts.push(Math.max(0,z));
    const speedCurr = Math.sqrt(vx*vx + wz*wz); vs.push(speedCurr); wzs.push(wz); times.push(t); ms.push(m);

    // After recording this sample, check disintegration (so UI can at least display this frame).
    if (!disintegrated && m <= MIN_FRACTION_FOR_CONTINUATION * mass0) {
      disintegrated = true;
      disintegrationIndex = xs.length - 1;
      disintegrationTime = t;
      // We stop generating further dynamics; break out leaving last sample as final.
      m = 0;
      break;
    }

    if (steps % REPORT_INTERVAL === 0) {
      const massPct = (m / mass0 * 100).toFixed(2);
      console.log(`[worker] step ${steps} t=${t.toFixed(1)}s alt=${(z/1000).toFixed(2)}km speed=${(speedCurr/1000).toFixed(2)}km/s mass=${massPct}%`);
    }

    const state = { x, z, vx, wz, m };
    let dt_use = adaptiveTimeStep(state, params, dtBase);
    if (z < z_adapt_threshold) dt_use = Math.min(dt_use, 0.001);
    let next = rk4Vec(state, dt_use, params);
    let retries = 0; const MAX_RETRIES = 6; const speedSafe = Math.max(speedCurr,1e-6);
    while (retries < MAX_RETRIES) {
      let unstable = false;
      if (next.z > state.z + Math.max(100, Math.abs(state.z)*0.01)) unstable = true;
      const nextSpeed = Math.sqrt(next.vx*next.vx + next.wz*next.wz);
      const ratio = nextSpeed / speedSafe; if (ratio > 1.5 || ratio < 0.5) unstable = true;
      const accelMag = Math.sqrt(Math.pow((next.vx - state.vx)/dt_use,2)+Math.pow((next.wz - state.wz)/dt_use,2));
      if (accelMag > 5000) unstable = true;
      if (!isFinite(next.x)||!isFinite(next.z)||!isFinite(next.vx)||!isFinite(next.wz)||!isFinite(next.m)) unstable = true;
      if (!unstable) break;
      dt_use *= 0.3; if (dt_use < 1e-6) { console.warn('[worker] dt underflow terminating'); break; }
      next = rk4Vec(state, dt_use, params); retries++;
    }
    if (next.z <= 0) {
      const frac = state.z / (state.z - next.z);
      const impX = lerp(state.x, next.x, frac);
      const impVx = lerp(state.vx, next.vx, frac);
      const impWz = lerp(state.wz, next.wz, frac);
      const impSpeed = Math.sqrt(impVx*impVx + impWz*impWz);
      const [ilat, ilon] = destinationLatLon(lat0, lon0, bearing, impX);
      const impTime = t + dt_use * frac;
  xs.push(impX); lats.push(ilat); lons.push(ilon); alts.push(0); vs.push(impSpeed); wzs.push(impWz); times.push(impTime); ms.push(m);
      t = impTime; impacted = true; impactTime = impTime; break;
    }
    x = next.x; z = next.z; vx = next.vx; wz = next.wz; m = Math.max(0, next.m); t += dt_use; steps++;
    // Simple rare fragmentation
    const rhoLocal = rhoAtmosphere(Math.max(0,z));
    const dynP = 0.5 * rhoLocal * (vx*vx + wz*wz);
    if (dynP > strength && Math.random() < 0.02) m *= 0.95;
  }

  const n = xs.length;
  const xsBuf = new Float64Array(xs).buffer;
  const latsBuf = new Float64Array(lats).buffer;
  const lonsBuf = new Float64Array(lons).buffer;
  const altsBuf = new Float32Array(alts).buffer;
  const vsBuf = new Float32Array(vs).buffer;
  const wzsBuf = new Float32Array(wzs).buffer;
  const timesBuf = new Float64Array(times).buffer;
  const msBuf = new Float64Array(ms).buffer;

  // Consequence modeling now uses RESIDUAL kinetic energy at impact (or last sample if no impact).
  // We still compute initial energy for diagnostics and ratio.
  const initialEk = 0.5 * mass0 * v0 * v0;
  const lastSpeed = vs.length ? vs[vs.length - 1] : 0;
  // 'm' here is the final (possibly partially ablated) mass from the integration loop scope.
  const residualEk = 0.5 * m * lastSpeed * lastSpeed;
  const EkForEffects = residualEk > 0 ? residualEk : initialEk; // fallback safety
  const crater_km = 1.435e-5 * Math.pow(EkForEffects, 1/3.4);
  const crater_m = crater_km * 1000;
  const tsunami_m = crater_m * 5;
  const k_eff = 5.7e-5; const Eseis = EkForEffects * k_eff; const mw = (Math.log10(Eseis) - 5.24)/1.44;

  return {
    xs: xsBuf, lats: latsBuf, lons: lonsBuf, alts: altsBuf, vs: vsBuf, wzs: wzsBuf, times: timesBuf, ms: msBuf,
    crater_km: disintegrated ? 0 : crater_km,
    mw: disintegrated ? 0 : mw,
    tsunami_m: disintegrated ? 0 : tsunami_m,
    meta: {
      impacted,
      disintegrated,
      disintegrationIndex,
      disintegrationTime,
      steps,
      impactTime,
      lastAlt: alts[n-1],
      lastTime: times[times.length-1],
      initialEkJ: initialEk,
      residualEkJ: residualEk,
      energyFraction: initialEk > 0 ? (residualEk / initialEk) : 1,
      initialMass: mass0,
      finalMass: m
    }
  };
}
