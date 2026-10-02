use serde::{Deserialize, Serialize};

use crate::backend::{Airframe, Atmosphere, FConf, TOConf, VSpeeds, FT_TO_M, G, KT_TO_MS, MS_TO_KT};
use crate::runway::RunwayLimits;
use crate::Weight;

const DT: f32 = 0.05;
const REACTION_TIME: f32 = 1.5; // sec
const AEO_FACTOR: f32 = 1.15;
const SCREEN_HEIGHT_M: f32 = 35.0 * FT_TO_M;
const V2_MARGIN_KT: f32 = 4.0;
const MU_AQUA: f32 = 0.04;
const ROLLING_MU: f32 = 0.02;
const MAX_STEPS: u32 = 400_000;
const EPS: f32 = 1e-3;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum RnwCondition {
    Dry, // ICAO 6
    Wet, // ICAO 5
    CompactedSnow, // ICAO 4
    Snow, // ICAO 3
    Slush, // ICAO 2
    Ice, // ICAO 1
    WetIce, // ICAO 0
}

impl RnwCondition {
    fn base_braking_mu(&self) -> f32 {
        match self {
            RnwCondition::Dry => 0.65,
            RnwCondition::Wet => 0.50,
            RnwCondition::CompactedSnow => 0.38,
            RnwCondition::Snow => 0.30,
            RnwCondition::Slush => 0.21,
            RnwCondition::Ice => 0.12,
            RnwCondition::WetIce => 0.04,
        }
    }

    pub fn is_contaminated(&self) -> bool {
        !matches!(self, RnwCondition::Dry | RnwCondition::Wet)
    }

    pub fn reverse_credit(&self) -> bool {
        !matches!(self, RnwCondition::Dry)
    }

    pub fn braking_mu(&self, gs_kt: f32, aquaplaning_kt: f32) -> f32 {
        let base = self.base_braking_mu();
        if matches!(self, RnwCondition::Wet | RnwCondition::Slush) {
            let r = (gs_kt / aquaplaning_kt).clamp(0.0, 1.0);
            base - r * r * (base - MU_AQUA)
        } else {
            base
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Phase {
    Roll,
    Braking,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum Limit {
    Balanced,    // ASDR == TODR(OEI)
    Vmcg,        // V1 bottom limit (Vmcg)
    Vr,          // V1 to Vr
    BrakeEnergy, // V1 to (VMBE)
    Toda,        // TODA pushed V1 up
    Asda,        // ASDA pushed V1 down
}

#[derive(Debug, Clone, Serialize)]
pub struct TakeoffResult {
    pub vef_kt: f32,
    pub v1_kt: f32,
    pub vr_kt: f32,
    pub v2_kt: f32,
    pub v1_gs_kt: f32,
    pub todr_m: f32,
    pub todr_oei_m: f32,
    pub todr_aeo_m: f32,
    pub asdr_m: f32,
    pub toda_margin_m: f32,
    pub asda_margin_m: f32,
    pub limit: Limit,
    pub reverse_credited: bool,
    pub speeds: VSpeeds,
}

#[derive(Clone, Copy)]
struct Speeds {
    vmcg: f32,
    vr: f32,
    v2: f32,
}

struct Eval {
    v1: f32,
    todr_oei: f32,
    asdr: f32,
}

struct Ctx<'a> {
    a: &'a Airframe,
    atm: &'a Atmosphere,
    rating_temp_c: f32,
    w: &'a Weight,
    bleed: &'a TOConf,
    flap: &'a FConf,
    rnw: RnwCondition,
    slope: f32,
}

impl<'a> Ctx<'a> {
    fn engines(&self, aeo: bool) -> f32 {
        let n = self.a.n_engines as f32;
        if aeo { n } else { n - 1.0 }
    }
    fn failed(&self, aeo: bool) -> u8 {
        if aeo { 0 } else { 1 }
    }
    fn thrust(&self, aeo: bool, tas: f32) -> f32 {
        self.a.engine.takeoff_thrust(self.atm, self.rating_temp_c, tas, self.bleed) * self.engines(aeo)
    }
    fn rev_thrust(&self, aeo: bool, tas: f32) -> f32 {
        if !self.rnw.reverse_credit() {
            return 0.0;
        }
        -self.a.engine.reverse_thrust(self.atm, tas, self.bleed) * self.engines(aeo)
    }
    fn mass(&self) -> f32 {
        self.w.tow_kg as f32
    }
    fn weight_n(&self) -> f32 {
        self.mass() * G
    }
}

fn crossing<F>(mut lo: f32, mut hi: f32, mut g: F) -> Result<(f32, f32), String>
where
    F: FnMut(f32) -> Result<f32, String>,
{
    if g(lo)? >= 0.0 {
        return Ok((lo, lo));
    }
    if g(hi)? <= 0.0 {
        return Ok((hi, hi));
    }
    for _ in 0..40 {
        if hi - lo < 0.01 {
            break;
        }
        let mid = 0.5 * (lo + hi);
        if g(mid)? < 0.0 { lo = mid } else { hi = mid }
    }
    Ok((lo, hi))
}

pub struct PhyEng {
    time: f32,
    distance: f32,
    gs: f32,
    altitude: f32,
    headwind_kt: f32,
}

impl PhyEng {
    pub fn new() -> Self {
        Self { time: 0.0, distance: 0.0, gs: 0.0, altitude: 0.0, headwind_kt: 0.0 }
    }

    fn airspeed(&self) -> f32 {
        (self.gs + self.headwind_kt).max(0.0)
    }
    fn set_airspeed(&mut self, v: f32) {
        self.gs = (v - self.headwind_kt).max(0.0);
    }
    fn reset(&mut self) {
        self.time = 0.0;
        self.distance = 0.0;
        self.gs = 0.0;
        self.altitude = 0.0;
    }
    fn start_at(&mut self, v_air: f32) {
        self.reset();
        self.set_airspeed(v_air);
    }

    fn ground_accel(&self, ctx: &Ctx, thrust: f32, phase: Phase, aeo: bool) -> f32 {
        let tas = self.airspeed();
        let qs = ctx.atm.dynamic_pressure(tas) * ctx.a.aero.wing_area_m2;
        let (cl, cd) = ctx.a.ground_coeffs(ctx.flap, phase == Phase::Braking, ctx.failed(aeo));
        let weight = ctx.weight_n();
        let normal = (weight - qs * cl).max(0.0);
        let mu = match phase {
            Phase::Roll => ROLLING_MU,
            Phase::Braking => ctx.rnw.braking_mu(self.gs, ctx.a.gear.aquaplaning_speed_kt()),
        };
        let net = thrust - qs * cd - mu * normal - weight * ctx.slope;
        net / ctx.mass() * MS_TO_KT
    }

    fn dist_to_speed(&mut self, target: f32, ctx: &Ctx, aeo: bool) -> Result<(), String> {
        let mut steps = 0u32;
        while target - self.airspeed() > EPS {
            steps += 1;
            if steps > MAX_STEPS {
                return Err("Speed calculation couldn't ended.".to_string());
            }
            let acc = self.ground_accel(ctx, ctx.thrust(aeo, self.airspeed()), Phase::Roll, aeo);
            if acc <= 0.0 {
                return Err("Unsufficent acceleration.".to_string());
            }
            let remaining = target - self.airspeed();
            let dt = if acc * DT > remaining { remaining / acc } else { DT };
            let gs_old = self.gs;
            self.gs += acc * dt;
            self.distance += 0.5 * (gs_old + self.gs) * KT_TO_MS * dt;
            self.time += dt;
        }
        Ok(())
    }

    fn dist_to_time(&mut self, duration: f32, ctx: &Ctx, aeo: bool) -> Result<(), String> {
        let end = self.time + duration;
        let mut steps = 0u32;
        while end - self.time > 1e-4 {
            steps += 1;
            if steps > MAX_STEPS {
                return Err("Reaction dist can't calculated".to_string());
            }
            let acc = self.ground_accel(ctx, ctx.thrust(aeo, self.airspeed()), Phase::Roll, aeo);
            let dt = DT.min(end - self.time);
            let gs_old = self.gs;
            self.gs = (self.gs + acc * dt).max(0.0);
            self.distance += 0.5 * (gs_old + self.gs) * KT_TO_MS * dt;
            self.time += dt;
        }
        Ok(())
    }

    fn dist_to_stop(&mut self, ctx: &Ctx, aeo: bool) -> Result<(), String> {
        let mut steps = 0u32;
        while self.gs > EPS {
            steps += 1;
            if steps > MAX_STEPS {
                return Err("Err.".to_string());
            }
            let acc = self.ground_accel(ctx, ctx.rev_thrust(aeo, self.airspeed()), Phase::Braking, aeo);
            if acc >= 0.0 {
                return Err("No decel.".to_string());
            }
            let dt = if self.gs + acc * DT < 0.0 { self.gs / -acc } else { DT };
            let gs_old = self.gs;
            self.gs = (self.gs + acc * dt).max(0.0);
            self.distance += 0.5 * (gs_old + self.gs) * KT_TO_MS * dt;
            self.time += dt;
        }
        Ok(())
    }

    fn dist_to_alt(&mut self, target_m: f32, ctx: &Ctx, aeo: bool) -> Result<(), String> {
        let tas = self.airspeed();
        let v_ms = tas * KT_TO_MS;
        let qs = ctx.atm.dynamic_pressure(tas) * ctx.a.aero.wing_area_m2;
        let weight = ctx.weight_n();
        let cl = weight / qs;
        let drag = qs * ctx.a.airborne_cd(ctx.flap, cl, true, ctx.failed(aeo));
        let sin_g = ((ctx.thrust(aeo, tas) - drag) / weight).min(1.0);
        if sin_g <= 0.0 {
            return Err("Unsufficent thrust.".to_string());
        }
        let cos_g = (1.0 - sin_g * sin_g).sqrt();
        let climb_rate = v_ms * sin_g;
        let ground_speed = (v_ms * cos_g - self.headwind_kt * KT_TO_MS).max(0.0);
        let t = (target_m - self.altitude) / climb_rate;
        self.altitude = target_m;
        self.distance += ground_speed * t;
        self.time += t;
        Ok(())
    }

    fn v1_after_reaction(&mut self, ctx: &Ctx, vef: f32) -> Result<f32, String> {
        self.start_at(vef);
        self.dist_to_time(REACTION_TIME, ctx, false)?;
        Ok(self.airspeed())
    }

    fn eval(&mut self, ctx: &Ctx, sp: Speeds, vef: f32) -> Result<Eval, String> {
        self.reset();
        self.dist_to_speed(vef, ctx, true)?;
        let until_vef = self.distance;

        self.start_at(vef);
        self.dist_to_time(REACTION_TIME, ctx, false)?;
        let reaction_dist = self.distance;
        let v1 = self.airspeed();

        self.start_at(vef);
        self.dist_to_speed(sp.vr, ctx, false)?;
        let vef_to_vr = self.distance;
        self.start_at(sp.v2);
        self.dist_to_alt(SCREEN_HEIGHT_M, ctx, false)?;
        let todr_oei = until_vef + vef_to_vr + self.distance;

        self.start_at(v1);
        self.dist_to_stop(ctx, false)?;
        let asdr_oei = until_vef + reaction_dist + self.distance;

        self.reset();
        self.dist_to_speed(v1, ctx, true)?;
        let until_v1_aeo = self.distance;
        self.start_at(v1);
        self.dist_to_stop(ctx, true)?;
        let asdr_aeo = until_v1_aeo + self.distance;

        Ok(Eval { v1, todr_oei, asdr: asdr_oei.max(asdr_aeo) })
    }

    #[allow(clippy::too_many_arguments)]
    pub fn calculate_takeoff(
        &mut self,
        a: &Airframe,
        atm: &Atmosphere,
        rating_temp_c: f32,
        w: &Weight,
        bleed: &TOConf,
        flap: &FConf,
        rnw: RnwCondition,
        rwy: &RunwayLimits,
    ) -> Result<TakeoffResult, String> {
        self.headwind_kt = rwy.headwind_kt;
        let ctx = Ctx { a, atm, rating_temp_c, w, bleed, flap, rnw, slope: rwy.slope };

        let vs = a.calculate_speeds(w, atm, flap)?;
        let vr_cas = vs.vr_min.max(vs.v2_min - V2_MARGIN_KT);
        let v2_cas = vs.v2_min.max(vr_cas + V2_MARGIN_KT);
        let sp = Speeds {
            vmcg: atm.cas_to_tas(vs.vmcg),
            vr: atm.cas_to_tas(vr_cas),
            v2: atm.cas_to_tas(v2_cas),
        };
        if sp.vmcg >= sp.vr {
            return Err("Vmcg >= Vr: speeds inconsistent.".to_string());
        }

        if let Some(lim) = a.gear.tire_speed_limit_kt {
            let lof_gs = sp.v2 - self.headwind_kt;
            if lof_gs > lim {
                return Err(format!("Takeoff ground speed ({lof_gs:.0} kt) exceeds tire speed limit ({lim:.0} kt)."));
            }
        }

        self.reset();
        self.dist_to_speed(sp.vr, &ctx, true)?;
        let until_vr = self.distance;
        self.start_at(sp.v2);
        self.dist_to_alt(SCREEN_HEIGHT_M, &ctx, true)?;
        let todr_aeo = AEO_FACTOR * (until_vr + self.distance);
        if todr_aeo > rwy.toda_m {
            return Err(format!(
                "AEO TODR ×1.15 ({todr_aeo:.0} m) exceeds TODA ({:.0} m).",
                rwy.toda_m
            ));
        }

        let vmin = sp.vmcg;
        if self.v1_after_reaction(&ctx, vmin)? >= sp.vr {
            return Err("Vmcg + reaction time distance exceeds vr.".to_string());
        }
        let (mut vmax, _) = crossing(vmin, sp.vr, |x| Ok(self.v1_after_reaction(&ctx, x)? - sp.vr))?;

        let mut brake_limited = false;
        if let Some(vmbe) = a.vmbe_gs_kt(w) {
            let hw = self.headwind_kt;
            if self.v1_after_reaction(&ctx, vmin)? - hw > vmbe {
                return Err(format!("Braking energy: even the lowest V1 exceeds VMBE ({vmbe:.0} kt GS)."));
            }
            let (vb, _) = crossing(vmin, vmax, |x| Ok(self.v1_after_reaction(&ctx, x)? - hw - vmbe))?;
            if vb < vmax {
                vmax = vb;
                brake_limited = true;
            }
        }

        let e_min = self.eval(&ctx, sp, vmin)?;
        let e_max = self.eval(&ctx, sp, vmax)?;
        if e_max.todr_oei > rwy.toda_m {
            return Err("Even at the highest V1 level, OEI TODR surpasses TODA.".to_string());
        }
        if e_min.asdr > rwy.asda_m {
            return Err("Even at the lowest V1 level, ASDR exceeds ASDA.".to_string());
        }

        let (bl, bh) = crossing(vmin, vmax, |x| {
            let e = self.eval(&ctx, sp, x)?;
            Ok(e.asdr - e.todr_oei)
        })?;
        let vef_bal = 0.5 * (bl + bh);

        let lo_v = if e_min.todr_oei <= rwy.toda_m {
            vmin
        } else {
            crossing(vmin, vmax, |x| Ok(rwy.toda_m - self.eval(&ctx, sp, x)?.todr_oei))?.1
        };
        let hi_v = if e_max.asdr <= rwy.asda_m {
            vmax
        } else {
            crossing(vmin, vmax, |x| Ok(self.eval(&ctx, sp, x)?.asdr - rwy.asda_m))?.0
        };
        if lo_v > hi_v {
            return Err("No V1 offers both TODA and ASDA simultaneously.".to_string());
        }

        let (vef, limit) = if vef_bal < lo_v {
            (lo_v, Limit::Toda)
        } else if vef_bal > hi_v {
            (hi_v, Limit::Asda)
        } else if vef_bal <= vmin + 0.02 {
            (vef_bal, Limit::Vmcg)
        } else if vef_bal >= vmax - 0.02 {
            (vef_bal, if brake_limited { Limit::BrakeEnergy } else { Limit::Vr })
        } else {
            (vef_bal, Limit::Balanced)
        };

        let e = self.eval(&ctx, sp, vef)?;
        let todr = e.todr_oei.max(todr_aeo);
        Ok(TakeoffResult {
            vef_kt: atm.tas_to_cas(vef),
            v1_kt: atm.tas_to_cas(e.v1),
            vr_kt: vr_cas,
            v2_kt: v2_cas,
            v1_gs_kt: e.v1 - self.headwind_kt,
            todr_m: todr,
            todr_oei_m: e.todr_oei,
            todr_aeo_m: todr_aeo,
            asdr_m: e.asdr,
            toda_margin_m: rwy.toda_m - todr,
            asda_margin_m: rwy.asda_m - e.asdr,
            limit,
            reverse_credited: rnw.reverse_credit(),
            speeds: vs,
        })
    }
}
