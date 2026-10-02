
use serde::{Deserialize, Serialize};

use crate::backend::{Airframe, Atmosphere, ClimbSegment, FConf, TOConf, TrimSetting};
use crate::phyeng::{PhyEng, RnwCondition, TakeoffResult};
use crate::runway::RunwayLimits;
use crate::Weight;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum FlapStrategy {
    MaxFlex,
    MaxPerformance,
    Fixed(String),
}

pub struct PlanInput<'a> {
    pub airframe: &'a Airframe,
    pub atmosphere: &'a Atmosphere,
    pub weight: &'a Weight,
    pub toconf: &'a TOConf,
    pub rnw: RnwCondition,
    pub runway: &'a RunwayLimits,
    pub strategy: FlapStrategy,
    pub allow_flex: bool,
}

#[derive(Debug, Clone, Serialize)]
pub struct Candidate {
    pub flap: String,
    pub flex_temp_c: Option<i8>,
    pub thrust_ratio: f32,
    pub first_seg_gradient_pct: f32,
    pub oei_gradient_pct: f32,
    pub pitch_trim: Option<TrimSetting>,
    pub takeoff: TakeoffResult,
}

#[derive(Debug, Clone, Serialize)]
pub struct Rejected {
    pub flap: String,
    pub reason: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct PlanResult {
    pub best: Candidate,
    pub feasible: Vec<Candidate>,
    pub rejected: Vec<Rejected>,
}


fn evaluate(inp: &PlanInput, flap: &FConf, flex: Option<i8>) -> Result<Candidate, String> {
    let a = inp.airframe;
    let atm = inp.atmosphere;
    let eng = &a.engine;
    a.check_weight(inp.weight)?;
    let pitch_trim = a.pitch_trim(inp.weight, flap)?;

    let oat = atm.oat_c;
    let rating = flex.map_or(oat, |f| f as f32);
    let full = eng.takeoff_thrust(atm, oat, 0.0, inp.toconf);
    let used = eng.takeoff_thrust(atm, rating, 0.0, inp.toconf);
    let ratio = used / full;
    if ratio < 1.0 - eng.max_thrust_reduction {
        return Err(format!("Thrust reduction exceeds %{:.0} limit.", eng.max_thrust_reduction * 100.0));
    }

    let mut sim = PhyEng::new();
    let takeoff =
        sim.calculate_takeoff(a, atm, rating, inp.weight, inp.toconf, flap, inp.rnw, inp.runway)?;

    let g1 = a.climb_gradient_pct(inp.weight, atm, rating, inp.toconf, flap, takeoff.v2_kt, ClimbSegment::First);
    let g2 = a.climb_gradient_pct(inp.weight, atm, rating, inp.toconf, flap, takeoff.v2_kt, ClimbSegment::Second);
    let r1 = a.required_gradient_pct(ClimbSegment::First);
    let r2 = a.required_gradient_pct(ClimbSegment::Second);
    if g1 <= r1 {
        return Err(format!("first segment gradient %{g1:.2} ≤ %{r1}."));
    }
    if g2 < r2 {
        return Err(format!("second segment gradient %{g2:.2} < %{r2}."));
    }

    let flex_temp_c = if ratio >= 0.9999 { None } else { flex };
    Ok(Candidate {
        flap: flap.name.clone(),
        flex_temp_c,
        thrust_ratio: ratio,
        first_seg_gradient_pct: g1,
        oei_gradient_pct: g2,
        pitch_trim,
        takeoff,
    })
}

fn best_for_flap(inp: &PlanInput, flap: &FConf, use_flex: bool) -> Result<Candidate, String> {
    let toga = evaluate(inp, flap, None)?;
    if !use_flex {
        return Ok(toga);
    }

    let oat = inp.atmosphere.oat_c.ceil() as i16;
    let t_max = inp.airframe.engine.max_flex_temp_c.floor() as i16;
    if t_max <= oat {
        return Ok(toga);
    }

    if let Ok(c) = evaluate(inp, flap, Some(t_max as i8)) {
        return Ok(c);
    }

    let (mut lo, mut hi) = (oat, t_max);
    while hi - lo > 1 {
        let mid = lo + (hi - lo) / 2;
        if evaluate(inp, flap, Some(mid as i8)).is_ok() { lo = mid } else { hi = mid }
    }
    if lo == oat {
        return Ok(toga);
    }
    evaluate(inp, flap, Some(lo as i8))
}

pub fn plan(inp: &PlanInput) -> Result<PlanResult, String> {
    let a = inp.airframe;
    let rwy = inp.runway;

    a.check_weight(inp.weight)?;

    if let Some(max_xw) = a.limits.max_crosswind_kt {
        let xw = rwy.crosswind_kt.abs();
        if xw > max_xw {
            return Err(format!("Cross wind {xw:.0} kt, limits are {max_xw:.0} kt."));
        }
    }
    if let Some(max_tw) = a.limits.max_tailwind_kt {
        let tw = -rwy.headwind_raw_kt;
        if tw > max_tw {
            return Err(format!("Tail wind {tw:.0} kt, limits are {max_tw:.0} kt."));
        }
    }

    let use_flex = inp.allow_flex
        && !inp.rnw.is_contaminated()
        && !matches!(inp.strategy, FlapStrategy::MaxPerformance);

    let flaps: Vec<&FConf> = match &inp.strategy {
        FlapStrategy::Fixed(name) => vec![a.flap(name)?],
        _ => a.flaps.iter().collect(),
    };

    let mut feasible = Vec::new();
    let mut rejected = Vec::new();
    for f in flaps {
        match best_for_flap(inp, f, use_flex) {
            Ok(c) => feasible.push(c),
            Err(reason) => rejected.push(Rejected { flap: f.name.clone(), reason }),
        }
    }

    if feasible.is_empty() {
        let why: Vec<String> = rejected.iter().map(|r| format!("[{}] {}", r.flap, r.reason)).collect();
        return Err(format!("No flaps suitable: {}", why.join(" | ")));
    }

    let best = match &inp.strategy {
        FlapStrategy::MaxPerformance => feasible
            .iter()
            .min_by(|x, y| x.takeoff.todr_m.total_cmp(&y.takeoff.todr_m)),
        _ => feasible.iter().min_by(|x, y| {
            x.thrust_ratio
                .total_cmp(&y.thrust_ratio)
                .then(y.oei_gradient_pct.total_cmp(&x.oei_gradient_pct))
        }),
    }
    .cloned()
    .expect("feasible is not empty");

    Ok(PlanResult { best, feasible, rejected })
}
