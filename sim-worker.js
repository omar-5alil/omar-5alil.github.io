// public/sim-worker.js
/* eslint-disable no-restricted-globals */
/**
 * Stable atmospheric entry worker using velocity components (vx, wz).
 * - writes arrays for x, lats, lons, alts, speed, wz
 * - interpolates a single final impact point at z=0 (no bounce)
 * - defensive retries on unrealistic upward steps, dt adapt near ground
 *
 * Use: const worker = new Worker('/sim-worker.js'); worker.postMessage({cmd:'run', params:{...}});
 *
 * params: { lat0, lon0, bearing, diameter, density, v0, angleDeg, alt0?, dt?, maxTime? }
 */

/**
 * @typedef {{ lat0:number, lon0:number, bearing:number, diameter:number, density:number, v0:number, angleDeg:number, alt0?:number, dt?:number, maxTime?:number }} SimParams
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
    crater_km: out.crater_km,
    mw: out.mw,
    tsunami_m: out.tsunami_m
  }, [out.xs, out.lats, out.lons, out.alts, out.vs, out.wzs]);
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
  // More realistic atmospheric model with multiple layers
  if (z > 86000) {
    // Thermosphere: very thin
    return 1.225 * Math.exp(-z / 40000) * 1e-6;
  } else if (z > 50000) {
    // Mesosphere: scale height ~7km
    const rho50 = 1.225 * Math.exp(-50000 / 8500);
    return rho50 * Math.exp(-(z - 50000) / 7000);
  } else {
    // Troposphere/Stratosphere: standard model
    const H = 8500;
    return 1.225 * Math.exp(-z / H);
  }
}

/**
 * derivsVec: for state {x,z,vx,wz,m}
 * returns dxdt, dzdt, ax, az, dmdt, speed
 */
function derivsVec(state, p) {
  const x = state.x, z = state.z, vx = state.vx, wz = state.wz, m = state.m;
  const v = Math.sqrt(vx*vx + wz*wz) || 1e-6;
  const A = Math.PI * (p.radius * p.radius);
  const rho = rhoAtmosphere(Math.max(0, z));
  
  // Velocity-dependent drag coefficient (more realistic)
  let Cd = p.Cd;
  if (v > 15000) {
    Cd = 0.3; // Lower drag at very high speeds
  } else if (v > 8000) {
    Cd = 0.5; // Moderate drag at high speeds
  } else {
    Cd = 1.0; // Standard drag at lower speeds
  }
  
  const Ch = p.Ch;
  const Q = p.Q;
  const g = p.g;

  // drag acceleration components
  const dragFactor = 0.5 * Cd * rho * A * v * v / m;
  const ax = - dragFactor * (vx / v);
  const az = - dragFactor * (wz / v) - g;

  // Realistic ablation model based on actual meteor physics
  // Most meteors lose only 10-90% of their mass during atmospheric entry
  // Only very small meteors completely ablate
  const dmdt = - (Ch * rho * A * v * v) / (2 * Q); // Changed from v³ to v²

  return { dxdt: vx, dzdt: wz, ax, az, dmdt, speed: v };
}

function rk4Vec(state, dt, p) {
  const s1 = derivsVec(state, p);
  const st2 = {
    x: state.x + s1.dxdt * dt/2,
    z: state.z + s1.dzdt * dt/2,
    vx: state.vx + s1.ax * dt/2,
    wz: state.wz + s1.az * dt/2,
    m: state.m + s1.dmdt * dt/2
  };
  const s2 = derivsVec(st2, p);

  const st3 = {
    x: state.x + s2.dxdt * dt/2,
    z: state.z + s2.dzdt * dt/2,
    vx: state.vx + s2.ax * dt/2,
    wz: state.wz + s2.az * dt/2,
    m: state.m + s2.dmdt * dt/2
  };
  const s3 = derivsVec(st3, p);

  const st4 = {
    x: state.x + s3.dxdt * dt,
    z: state.z + s3.dzdt * dt,
    vx: state.vx + s3.ax * dt,
    wz: state.wz + s3.az * dt,
    m: state.m + s3.dmdt * dt
  };
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

function clamp(value, lo, hi) { return Math.max(lo, Math.min(hi, value)); }
function lerp(a, b, f) { return a + (b - a) * f; }

// Adaptive time step calculation based on physics
function adaptiveTimeStep(state, p, dt_base) {
  const v = Math.sqrt(state.vx * state.vx + state.wz * state.wz);
  const rho = rhoAtmosphere(Math.max(0, state.z));
  
  // Much more conservative time stepping
  let dt_adaptive = dt_base;
  
  // Very small time steps for high velocities
  if (v > 25000) {
    dt_adaptive = Math.min(dt_adaptive, 0.0005);
  } else if (v > 20000) {
    dt_adaptive = Math.min(dt_adaptive, 0.001);
  } else if (v > 15000) {
    dt_adaptive = Math.min(dt_adaptive, 0.002);
  } else if (v > 10000) {
    dt_adaptive = Math.min(dt_adaptive, 0.005);
  }
  
  // Altitude-based time step constraints
  if (state.z > 80000) {
    dt_adaptive = Math.min(dt_adaptive, 0.01);
  } else if (state.z > 50000) {
    dt_adaptive = Math.min(dt_adaptive, 0.005);
  } else if (state.z > 20000) {
    dt_adaptive = Math.min(dt_adaptive, 0.002);
  }
  
  // Ensure we don't change altitude by more than 50m per step
  if (Math.abs(state.wz) > 0) {
    const dt_alt = 50.0 / Math.abs(state.wz);
    dt_adaptive = Math.min(dt_adaptive, dt_alt);
  }
  
  return Math.max(1e-6, dt_adaptive);
}

function runSimulation(pIn) {
  // unpack
  const lat0 = pIn.lat0;
  const lon0 = pIn.lon0;
  const bearing = pIn.bearing;
  const diameter = pIn.diameter;
  const density = pIn.density;
  let v0 = pIn.v0;
  const angleDeg = pIn.angleDeg;
  let dt = (pIn.dt && pIn.dt > 0) ? pIn.dt : 0.02;
  const maxTime = pIn.maxTime || 400;

  // params
  const radius = diameter / 2;
  const area = Math.PI * radius * radius;
  const mass0 = (4/3) * Math.PI * Math.pow(radius, 3) * density;
  const Cd = 1.0;
  const Ch = 0.01; // Much lower ablation coefficient (was 0.05)
  const Q = 2e7; // Higher heat of ablation (was 1e7)
  const g = 9.81;
  const strength = 10e6; // Higher material strength (was 5e6)

  const startAlt = (typeof pIn.alt0 === 'number') ? pIn.alt0 : 100000;
  const ang = toRad(Math.abs(angleDeg));
  let vx = v0 * Math.cos(ang);
  let wz = - v0 * Math.sin(ang); // vertical velocity, negative downward
  let x = 0;
  let z = startAlt;
  let m = mass0;

  const maxSteps = Math.ceil(maxTime / Math.max(dt, 1e-6)) + 40;
  const xsArr = new Float64Array(maxSteps);
  const latsArr = new Float64Array(maxSteps);
  const lonsArr = new Float64Array(maxSteps);
  const altsArr = new Float32Array(maxSteps);
  const vsArr = new Float32Array(maxSteps);
  const wzsArr = new Float32Array(maxSteps);
  let step = 0;

  const z_adapt_threshold = 10000; // Increased from 2000m to 10000m
  const dt_near_surface = Math.min(dt, 0.001); // Smaller time step near surface

  for (let iter = 0; iter < maxSteps; iter++) {
    // Only stop if altitude reaches ground, don't stop for low mass
    if (z <= 0) break;
    
    // Prevent mass from going to zero but don't stop simulation
    if (m < 0.001 * mass0) {
      m = 0.001 * mass0; // Keep minimum 0.1% of original mass
    }

    // store
    const [plat, plon] = destinationLatLon(lat0, lon0, bearing, x);
    xsArr[step] = x;
    latsArr[step] = plat;
    lonsArr[step] = plon;
    altsArr[step] = Math.max(0, z);
    const speedCurr = Math.sqrt(vx*vx + wz*wz);
    vsArr[step] = speedCurr;
    wzsArr[step] = wz;
    
    // Debug mass loss every 100 steps
    if (step % 100 === 0) {
      const massPercent = (m / mass0 * 100).toFixed(1);
      console.log(`Step ${step}: Alt=${(z/1000).toFixed(1)}km, Speed=${(speedCurr/1000).toFixed(1)}km/s, Mass=${massPercent}%`);
    }

    // Adaptive time step calculation
    const state = { x, z, vx, wz, m };
    let dt_use = adaptiveTimeStep(state, { radius, Cd, Ch, Q, g }, dt);
    
    // Further reduce near surface
    if (z < z_adapt_threshold) {
      dt_use = Math.min(dt_use, dt_near_surface);
    }

    let next = rk4Vec(state, dt_use, { radius, Cd, Ch, Q, g });
    let retries = 0;
    const maxRetries = 5;

    // Stability checks with multiple retries
    while (retries < maxRetries) {
      let unstable = false;
      
      // Check for unrealistic upward jumps (much more conservative)
      if (next.z > state.z + Math.max(100, Math.abs(state.z)*0.01)) {
        unstable = true;
      }
      
      // Check for unrealistic speed changes (more conservative)
      const nextSpeed = Math.sqrt(next.vx*next.vx + next.wz*next.wz);
      const speedRatio = nextSpeed / Math.max(speedCurr, 1e-6);
      if (speedRatio > 1.5 || speedRatio < 0.5) {
        unstable = true;
      }
      
      // Check for excessive acceleration
      const accelMag = Math.sqrt(
        Math.pow((next.vx - state.vx)/dt_use, 2) + 
        Math.pow((next.wz - state.wz)/dt_use, 2)
      );
      if (accelMag > 5000) { // 5000 m/s² max acceleration
        unstable = true;
      }
      
      // Check for NaN or infinite values
      if (!isFinite(next.x) || !isFinite(next.z) || !isFinite(next.vx) || 
          !isFinite(next.wz) || !isFinite(next.m)) {
        unstable = true;
      }
      
      if (!unstable) break;
      
      // Reduce time step more aggressively and retry
      dt_use *= 0.3;
      if (dt_use < 1e-6) {
        console.warn('Time step too small, ending simulation');
        break;
      }
      next = rk4Vec(state, dt_use, { radius, Cd, Ch, Q, g });
      retries++;
    }

    // if next crosses below ground -> interpolate exact impact and append single final hit
    if (next.z <= 0) {
      const frac = state.z / (state.z - next.z); // between 0 and 1
      const impX = lerp(state.x, next.x, frac);
      const impVx = lerp(state.vx, next.vx, frac);
      const impWz = lerp(state.wz, next.wz, frac);
      const impM = lerp(state.m, next.m, frac);
      const impSpeed = Math.sqrt(impVx*impVx + impWz*impWz);
      const [ilat, ilon] = destinationLatLon(lat0, lon0, bearing, impX);

      // append final
      xsArr[step + 1] = impX;
      latsArr[step + 1] = ilat;
      lonsArr[step + 1] = ilon;
      altsArr[step + 1] = 0;
      vsArr[step + 1] = impSpeed;
      wzsArr[step + 1] = impWz;
      step = step + 2;
      break;
    }

    // accept step
    x = next.x; z = next.z; vx = next.vx; wz = next.wz; m = Math.max(0.01 * mass0, next.m);

    // Much more conservative fragmentation model
    // Only occurs at very high dynamic pressures and with low probability
    const rho = rhoAtmosphere(Math.max(0, z));
    const v_total = Math.sqrt(vx*vx + wz*wz);
    const dynP = 0.5 * rho * v_total * v_total;
    
    // Only fragment if dynamic pressure exceeds strength AND random chance
    if (dynP > strength && Math.random() < 0.02) { // Only 2% chance per time step
      m *= 0.95; // Very gentle mass loss - only 5% per fragmentation event
    }

    step++;
    if (step >= maxSteps - 4) break;
  }

  // If we exited without crossing z<=0 but we have data and last alt > 0, attempt a fallback interpolation:
  const usedBeforeSlice = Math.min(step, xsArr.length);
  let used = usedBeforeSlice;
  if (used > 1) {
    const lastAlt = altsArr[used-1];
    if (lastAlt > 0) {
      // try to estimate impact using last state's vertical velocity if it's downward
      const lastWz = wzsArr[used-1];
      const lastX = xsArr[used-1];
      const lastV = vsArr[used-1];
      if (lastWz < -0.5) {
        const t_hit = lastAlt / (-lastWz);
        const estImpX = lastX + ( (vsArr[used-1] > 0) ? ( (xsArr[used-1] - (used>1?xsArr[used-2]:0)) / dt ) * t_hit : 0 );
        // fallback: simpler — linear extrap between last two points to hit zero
        const i1 = Math.max(0, used-2);
        const z0 = altsArr[i1];
        const z1 = altsArr[used-1];
        if (z1 < z0) {
          const frac = z0 / (z0 - z1);
          const impX = lerp(xsArr[i1], xsArr[used-1], frac);
          const [ilat, ilon] = destinationLatLon(lat0, lon0, bearing, impX);
          xsArr[used] = impX;
          latsArr[used] = ilat;
          lonsArr[used] = ilon;
          altsArr[used] = 0;
          vsArr[used] = vsArr[used-1];
          wzsArr[used] = lastWz;
          used = used + 1;
        } else {
          // as last resort, force last point altitude to 0 (hide small residual)
          altsArr[used-1] = 0;
        }
      } else {
        // lastWz not downward or too small -> force last sample alt to 0 (UI-level sanitization)
        altsArr[used-1] = 0;
      }
    }
  }

  const finalUsed = Math.min(used, xsArr.length);
  const xs = xsArr.slice(0, finalUsed).buffer;
  const lats = latsArr.slice(0, finalUsed).buffer;
  const lons = lonsArr.slice(0, finalUsed).buffer;
  const alts = altsArr.slice(0, finalUsed).buffer;
  const vs = vsArr.slice(0, finalUsed).buffer;
  const wzs = wzsArr.slice(0, finalUsed).buffer;

  // crater & seismic
  const Ek = 0.5 * mass0 * v0 * v0;
  const crater_km = 1.435e-5 * Math.pow(Ek, 1 / 3.4);
  const crater_m = crater_km * 1000;
  const k_eff = 5.7e-5;
  const Eseis = Ek * k_eff;
  const mw = (Math.log10(Eseis) - 5.24) / 1.44;
  const tsunami_m = crater_m * 5;

  return { xs, lats, lons, alts, vs, wzs, crater_km, mw, tsunami_m };
}
