use std::collections::HashSet;
use std::fmt;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::Weight;

pub const G: f32 = 9.80665;
pub const RHO_0: f32 = 1.225; // ISA density ASL (kg/m³)
pub const P_0: f32 = 101_325.0; // Pa
pub const T_0: f32 = 288.15; // K
pub const A_0: f32 = 340.294; // ISA mach 1 ASL (m/s)
pub const R_AIR: f32 = 287.053;
pub const GAMMA: f32 = 1.4;
pub const ISA_LAPSE: f32 = 0.0065; // K/m
pub const KT_TO_MS: f32 = 0.514444; // kts to m/s
pub const MS_TO_KT: f32 = 1.943844; // m/s to kts
pub const FT_TO_M: f32 = 0.3048; // feet to m

#[derive(Debug, Clone, PartialEq)]
pub enum PerfError {
    Io(String),
    Parse(String),
    InvalidData(Vec<String>),
    UnknownFlap(String),
    Overweight { tow_kg: f32, mtow_kg: f32 },
    CgOutOfRange { cg: f32, min: f32, max: f32 },
    AtmosphereOutOfRange(String),
    Trim(String),
}

impl fmt::Display for PerfError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PerfError::Io(s) => write!(f, "Couldn't read the JSON file: {s}"),
            PerfError::Parse(s) => write!(f, "JSON couldn't be resolved: {s}"),
            PerfError::InvalidData(v) => write!(f, "Invalid data: {}", v.join("; ")),
            PerfError::UnknownFlap(n) => write!(f, "Unknow flap config: {n}"),
            PerfError::Overweight { tow_kg, mtow_kg } => {
                write!(f, "TOW ({tow_kg:.0} kg) is out of MTOW ({mtow_kg:.0} kg) range.")
            }
            PerfError::CgOutOfRange { cg, min, max } => {
                write!(f, "CG %{cg:.1} MAC, out of range (%{min:.1}–%{max:.1}).")
            }
            PerfError::AtmosphereOutOfRange(s) => write!(f, "Atmosphere data is out of range: {s}"),
            PerfError::Trim(s) => write!(f, "Pitch trim: {s}"),
        }
    }
}

impl std::error::Error for PerfError {}

impl From<PerfError> for String {
    fn from(e: PerfError) -> Self {
        e.to_string()
    }
}

//Atmo

#[derive(Debug, Clone, Copy, Serialize)]
pub struct Atmosphere {
    pub qnh_hpa: f32,
    pub oat_c: f32,
    pub elevation_ft: f32,
    pub pressure_pa: f32,     // alan (istasyon) basıncı
    pub pressure_alt_ft: f32, // basınç irtifası
    pub isa_temp_c: f32,      // basınç irtifasındaki ISA sıcaklığı
    pub isa_dev_c: f32,       // OAT - ISA
    pub density_alt_ft: f32,
    pub rho: f32,
    pub sigma: f32, // ρ/ρ0
    pub delta: f32, // p/p0
    pub theta: f32, // T/T0
    pub speed_of_sound_ms: f32,
}

impl Atmosphere {
    pub fn new(qnh_hpa: f32, oat_c: f32, elevation_ft: f32) -> Result<Self, PerfError> {
        if !(870.0..=1090.0).contains(&qnh_hpa) {
            return Err(PerfError::AtmosphereOutOfRange(format!("QNH {qnh_hpa} hPa")));
        }
        if !(-60.0..=60.0).contains(&oat_c) {
            return Err(PerfError::AtmosphereOutOfRange(format!("OAT {oat_c} °C")));
        }
        if !(-1500.0..=16_000.0).contains(&elevation_ft) {
            return Err(PerfError::AtmosphereOutOfRange(format!("Elevation {elevation_ft} ft")));
        }

        let h_m = elevation_ft * FT_TO_M;
        let pressure_pa = qnh_hpa * 100.0 * (1.0 - ISA_LAPSE * h_m / T_0).powf(5.25588);
        let delta = pressure_pa / P_0;
        let pressure_alt_ft = T_0 / ISA_LAPSE * (1.0 - delta.powf(0.190263)) / FT_TO_M;

        let t_k = oat_c + 273.15;
        let rho = pressure_pa / (R_AIR * t_k);
        let sigma = rho / RHO_0;
        let isa_temp_c = 15.0 - ISA_LAPSE * pressure_alt_ft * FT_TO_M;
        let density_alt_ft = T_0 / ISA_LAPSE * (1.0 - sigma.powf(0.234969)) / FT_TO_M;

        Ok(Self {
            qnh_hpa,
            oat_c,
            elevation_ft,
            pressure_pa,
            pressure_alt_ft,
            isa_temp_c,
            isa_dev_c: oat_c - isa_temp_c,
            density_alt_ft,
            rho,
            sigma,
            delta,
            theta: t_k / T_0,
            speed_of_sound_ms: (GAMMA * R_AIR * t_k).sqrt(),
        })
    }

    // q = ½ρV² (Pa)
    pub fn dynamic_pressure(&self, tas_kt: f32) -> f32 {
        0.5 * self.rho * (tas_kt * KT_TO_MS).powi(2)
    }

    pub fn mach(&self, tas_kt: f32) -> f32 {
        tas_kt * KT_TO_MS / self.speed_of_sound_ms
    }

    //  TAS -> CAS
    pub fn tas_to_cas(&self, tas_kt: f32) -> f32 {
        let m = self.mach(tas_kt);
        let qc = self.pressure_pa * ((1.0 + 0.2 * m * m).powf(3.5) - 1.0);
        A_0 * (5.0 * ((qc / P_0 + 1.0).powf(2.0 / 7.0) - 1.0)).sqrt() * MS_TO_KT
    }

    // CAS -> TAS
    pub fn cas_to_tas(&self, cas_kt: f32) -> f32 {
        let c = cas_kt * KT_TO_MS / A_0;
        let qc = P_0 * ((1.0 + 0.2 * c * c).powf(3.5) - 1.0);
        let m = (5.0 * ((qc / self.pressure_pa + 1.0).powf(2.0 / 7.0) - 1.0)).sqrt();
        m * self.speed_of_sound_ms * MS_TO_KT
    }
}

#[allow(non_snake_case)]
#[derive(Debug, Clone, Copy, Default, Deserialize, Serialize)]
pub struct TOConf {
    pub Packs: bool,
    pub EngAI: bool,
    pub WingAI: bool,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(default)]
pub struct BleedPenalties {
    pub packs: f32,
    pub eng_ai: f32,
    pub wing_ai: f32,
}

impl Default for BleedPenalties {
    fn default() -> Self {
        Self { packs: 0.025, eng_ai: 0.007, wing_ai: 0.033 }
    }
}

impl BleedPenalties {
    pub fn factor(&self, c: &TOConf) -> f32 {
        let mut loss = 0.0;
        if c.Packs {
            loss += self.packs;
        }
        if c.EngAI {
            loss += self.eng_ai;
        }
        if c.WingAI {
            loss += self.wing_ai;
        }
        (1.0 - loss).max(0.5)
    }
}

fn d_temp_lapse() -> f32 { 0.009 }
fn d_pressure_exp() -> f32 { 0.9 }
fn d_mach_a() -> f32 { 0.45 }
fn d_mach_b() -> f32 { 0.35 }
fn d_max_flex() -> f32 { 68.0 }
fn d_max_red() -> f32 { 0.25 }
fn d_rev_ratio() -> f32 { 0.20 }
fn d_rev_cut() -> f32 { 70.0 }

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Engine {
    pub rated_thrust_n: f32,
    pub flat_rate_c: f32,
    #[serde(default = "d_temp_lapse")]
    pub temp_lapse_per_c: f32,
    #[serde(default = "d_pressure_exp")]
    pub pressure_exp: f32,
    #[serde(default = "d_mach_a")]
    pub mach_lapse_a: f32,
    #[serde(default = "d_mach_b")]
    pub mach_lapse_b: f32,
    #[serde(default = "d_max_flex")]
    pub max_flex_temp_c: f32,
    #[serde(default = "d_max_red")]
    pub max_thrust_reduction: f32,
    #[serde(default = "d_rev_ratio")]
    pub reverse_ratio: f32,
    #[serde(default = "d_rev_cut")]
    pub reverse_cutoff_kt: f32,
    #[serde(default)]
    pub bleed: BleedPenalties,
}

impl Engine {
    pub fn corner_temp_c(&self, atm: &Atmosphere) -> f32 {
        self.flat_rate_c - ISA_LAPSE * atm.pressure_alt_ft * FT_TO_M
    }

    pub fn takeoff_thrust(&self, atm: &Atmosphere, rating_temp_c: f32, tas_kt: f32, bleed: &TOConf) -> f32 {
        let t = rating_temp_c.max(atm.oat_c);
        let corner = self.corner_temp_c(atm);
        let temp_factor = if t <= corner {
            1.0
        } else {
            (1.0 - self.temp_lapse_per_c * (t - corner)).max(0.3)
        };
        let alt_factor = atm.delta.powf(self.pressure_exp);
        let m = atm.mach(tas_kt);
        let mach_factor = (1.0 - self.mach_lapse_a * m + self.mach_lapse_b * m * m).clamp(0.5, 1.0);
        self.rated_thrust_n * temp_factor * alt_factor * mach_factor * self.bleed.factor(bleed)
    }

    pub fn reverse_thrust(&self, atm: &Atmosphere, tas_kt: f32, bleed: &TOConf) -> f32 {
        if tas_kt < self.reverse_cutoff_kt {
            return 0.0;
        }
        self.takeoff_thrust(atm, atm.oat_c, 0.0, bleed) * self.reverse_ratio
    }
}


fn d_wing_h() -> f32 { 3.0 }
fn d_cd0() -> f32 { 0.022 }
fn d_cd_gear() -> f32 { 0.015 }
fn d_cd_wm() -> f32 { 0.0055 }
fn d_cd_trim() -> f32 { 0.0035 }
fn d_cd_spl() -> f32 { 0.05 }
fn d_cl_dump() -> f32 { 0.4 }
fn d_cg_ref() -> f32 { 25.0 }
fn d_cg_sens() -> f32 { 0.002 }

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Aero {
    pub wing_area_m2: f32,
    pub wingspan_m: f32,
    pub mac_m: f32,
    #[serde(default = "d_wing_h")]
    pub wing_height_m: f32,
    #[serde(default = "d_cd0")]
    pub cd0_clean: f32,
    #[serde(default = "d_cd_gear")]
    pub cd_gear: f32,
    #[serde(default = "d_cd_wm")]
    pub cd_windmill: f32,
    #[serde(default = "d_cd_trim")]
    pub cd_oei_trim: f32,
    #[serde(default = "d_cd_spl")]
    pub cd_spoilers: f32,
    #[serde(default = "d_cl_dump")]
    pub cl_lift_dump: f32,
    #[serde(default = "d_cg_ref")]
    pub cg_ref_pct: f32,
    #[serde(default = "d_cg_sens")]
    pub cg_tail_load_sens: f32,
}

impl Aero {
    pub fn ground_effect(&self, height_m: f32) -> f32 {
        let x = (16.0 * height_m / self.wingspan_m).powi(2);
        x / (1.0 + x)
    }
}

fn d_cl_roll() -> f32 { 0.3 }
fn d_one() -> f32 { 1.0 }

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct FConf {
    pub name: String,
    pub cl_max: f32,
    pub cl_unstick: f32,
    #[serde(default = "d_cl_roll")]
    pub cl_roll: f32,
    pub cd_flap: f32,
    pub k_induced: f32,
    #[serde(default = "d_one")]
    pub vmca_factor: f32,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Control {
    pub engine_arm_m: f32,
    pub fin_area_m2: f32,
    pub fin_arm_m: f32,
    pub fin_cy_max: f32,
    pub rudder_deriv: f32,
    pub rudder_max_rad: f32,
    pub weathercock_deriv: f32,
    pub max_sideslip_rad: f32,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Gear {
    pub tire_pressure_psi: f32,
    #[serde(default)]
    pub tire_speed_limit_kt: Option<f32>,
    #[serde(default)]
    pub max_brake_energy_mj: Option<f32>,
}

impl Gear {
    pub fn aquaplaning_speed_kt(&self) -> f32 {
        9.0 * self.tire_pressure_psi.sqrt()
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Limits {
    pub mtow_kg: f32,
    #[serde(default)]
    pub cg_min_pct: Option<f32>,
    #[serde(default)]
    pub cg_max_pct: Option<f32>,
    #[serde(default)]
    pub max_crosswind_kt: Option<f32>,
    #[serde(default)]
    pub max_tailwind_kt: Option<f32>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum TrimUnit {
    Deg,
    Units,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct TrimTable {
    #[serde(default)]
    pub flaps: Vec<String>,
    #[serde(default)]
    pub weight_kg: Option<f32>,
    /// [CG % MAC, trim]
    pub points: Vec<[f32; 2]>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct PitchTrim {
    pub unit: TrimUnit,
    /// allowed trim range
    #[serde(default)]
    pub min: Option<f32>,
    #[serde(default)]
    pub max: Option<f32>,
    pub tables: Vec<TrimTable>,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct TrimSetting {
    pub value: f32,
    pub unit: TrimUnit,
}

impl fmt::Display for TrimSetting {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.unit {
            TrimUnit::Deg => {
                let dir = if self.value >= 0.0 { "UP" } else { "DN" };
                write!(f, "{dir} {:.1}°", self.value.abs())
            }
            TrimUnit::Units => write!(f, "{:.2} units", self.value),
        }
    }
}

/// Linear interpolation with given trim values.
fn interp_linear(points: &[[f32; 2]], x: f32) -> Option<f32> {
    let first = points.first()?;
    let last = points.last()?;
    if x < first[0] || x > last[0] {
        return None;
    }
    for p in points.windows(2) {
        let (a, b) = (p[0], p[1]);
        if x <= b[0] {
            let t = (x - a[0]) / (b[0] - a[0]);
            return Some(a[1] + t * (b[1] - a[1]));
        }
    }
    Some(last[1])
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct VSpeeds {
    pub vsr: f32,
    pub vmu: f32,
    pub vmcg: f32,
    pub vmca: f32,
    pub vr_min: f32,
    pub v2_min: f32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ClimbSegment {
    /// VLOF
    First,
    /// V2 until 400 feet
    Second,
}

fn d_n_eng() -> u8 { 2 }

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Airframe {
    #[serde(default)]
    pub name: String,
    #[serde(default = "d_n_eng")]
    pub n_engines: u8,
    pub engine: Engine,
    pub aero: Aero,
    pub control: Control,
    pub gear: Gear,
    pub limits: Limits,
    pub flaps: Vec<FConf>,
    #[serde(default)]
    pub pitch_trim: Option<PitchTrim>,
}

fn check_pos(errs: &mut Vec<String>, name: &str, v: f32) {
    if !(v.is_finite() && v > 0.0) {
        errs.push(format!("{name} must be positive (value: {v})"));
    }
}

fn check_range(errs: &mut Vec<String>, name: &str, v: f32, lo: f32, hi: f32) {
    if !(v.is_finite() && v >= lo && v <= hi) {
        errs.push(format!("{name} must be in range of {lo}..{hi} (value: {v})"));
    }
}

impl Airframe {
    pub fn get_from_json<P: AsRef<Path>>(path: P) -> Result<Self, PerfError> {
        let file = File::open(path).map_err(|e| PerfError::Io(e.to_string()))?;
        let a: Self = serde_json::from_reader(BufReader::new(file)).map_err(|e| PerfError::Parse(e.to_string()))?;
        a.validate()?;
        Ok(a)
    }

    pub fn from_json_str(json: &str) -> Result<Self, PerfError> {
        let a: Self = serde_json::from_str(json).map_err(|e| PerfError::Parse(e.to_string()))?;
        a.validate()?;
        Ok(a)
    }

    pub fn validate(&self) -> Result<(), PerfError> {
        let mut e = Vec::new();

        if self.n_engines < 2 {
            e.push("n_engines(number of engines) must be at least 2 (for OEI gradient and to data.)".to_string());
        }

        let en = &self.engine;
        check_pos(&mut e, "engine.rated_thrust_n", en.rated_thrust_n);
        check_range(&mut e, "engine.flat_rate_c", en.flat_rate_c, -20.0, 60.0);
        check_range(&mut e, "engine.temp_lapse_per_c", en.temp_lapse_per_c, 0.0, 0.05);
        check_range(&mut e, "engine.max_thrust_reduction", en.max_thrust_reduction, 0.0, 0.5);
        check_range(&mut e, "engine.reverse_ratio", en.reverse_ratio, 0.0, 1.0);

        let ae = &self.aero;
        check_pos(&mut e, "aero.wing_area_m2", ae.wing_area_m2);
        check_pos(&mut e, "aero.wingspan_m", ae.wingspan_m);
        check_pos(&mut e, "aero.mac_m", ae.mac_m);
        check_pos(&mut e, "aero.wing_height_m", ae.wing_height_m);
        check_pos(&mut e, "aero.cd0_clean", ae.cd0_clean);

        let c = &self.control;
        check_pos(&mut e, "control.engine_arm_m", c.engine_arm_m);
        check_pos(&mut e, "control.fin_area_m2", c.fin_area_m2);
        check_pos(&mut e, "control.fin_arm_m", c.fin_arm_m);
        check_pos(&mut e, "control.fin_cy_max", c.fin_cy_max);
        check_pos(
            &mut e,
            "control yaw (rudder_deriv·rudder_max_rad + weathercock_deriv·max_sideslip_rad)",
            c.rudder_deriv * c.rudder_max_rad + c.weathercock_deriv * c.max_sideslip_rad,
        );

        check_pos(&mut e, "gear.tire_pressure_psi", self.gear.tire_pressure_psi);
        check_pos(&mut e, "limits.mtow_kg", self.limits.mtow_kg);
        if let (Some(lo), Some(hi)) = (self.limits.cg_min_pct, self.limits.cg_max_pct) {
            if lo >= hi {
                e.push(format!("limits.cg_min_pct ({lo}) must be lower than cg_max_pct'ten ({hi})."));
            }
        }

        if self.flaps.is_empty() {
            e.push("flap list is empty".to_string());
        }
        let mut names = HashSet::new();
        for f in &self.flaps {
            let p = format!("flap '{}'", f.name);
            if !names.insert(f.name.as_str()) {
                e.push(format!("{p}: occures more than one time."));
            }
            check_pos(&mut e, &format!("{p}.cl_max"), f.cl_max);
            check_pos(&mut e, &format!("{p}.cl_unstick"), f.cl_unstick);
            check_pos(&mut e, &format!("{p}.k_induced"), f.k_induced);
            check_range(&mut e, &format!("{p}.cd_flap"), f.cd_flap, 0.0, 0.2);
            check_range(&mut e, &format!("{p}.vmca_factor"), f.vmca_factor, 0.8, 1.2);
            if f.cl_unstick > f.cl_max {
                e.push(format!("{p}: cl_unstick, can't be higher than cl_max."));
            }
            if f.cl_roll >= f.cl_unstick {
                e.push(format!("{p}: cl_roll, must be lower than cl_unstick."));
            }
        }

        if let Some(pt) = &self.pitch_trim {
            if pt.tables.is_empty() {
                e.push("pitch_trim.tables is empty.".to_string());
            }
            if let (Some(lo), Some(hi)) = (pt.min, pt.max) {
                if lo >= hi {
                    e.push("pitch_trim.min, must be lower than pitch_trim.max.".to_string());
                }
            }
            for (i, t) in pt.tables.iter().enumerate() {
                if t.points.len() < 2 {
                    e.push(format!("pitch_trim.tables[{i}]: requires at least 2 points. "));
                }
                if t.points.windows(2).any(|p| p[1][0] <= p[0][0]) {
                    e.push(format!("pitch_trim.tables[{i}]: CG values must be incremental. "));
                }
                for n in &t.flaps {
                    if !self.flaps.iter().any(|f| &f.name == n) {
                        e.push(format!("pitch_trim.tables[{i}]: unknown flap '{n}'"));
                    }
                }
            }
        }

        if e.is_empty() { Ok(()) } else { Err(PerfError::InvalidData(e)) }
    }

    pub fn flap(&self, name: &str) -> Result<&FConf, PerfError> {
        self.flaps
            .iter()
            .find(|f| f.name == name)
            .ok_or_else(|| PerfError::UnknownFlap(name.to_string()))
    }

    pub fn check_weight(&self, w: &Weight) -> Result<(), PerfError> {
        let tow = w.tow_kg as f32;
        if tow > self.limits.mtow_kg {
            return Err(PerfError::Overweight { tow_kg: tow, mtow_kg: self.limits.mtow_kg });
        }
        let lo = self.limits.cg_min_pct.unwrap_or(f32::NEG_INFINITY);
        let hi = self.limits.cg_max_pct.unwrap_or(f32::INFINITY);
        if w.cg < lo || w.cg > hi {
            return Err(PerfError::CgOutOfRange { cg: w.cg, min: lo, max: hi });
        }
        Ok(())
    }

    fn effective_weight_n(&self, w: &Weight) -> f32 {
        w.tow_kg as f32 * G * (1.0 + (self.aero.cg_ref_pct - w.cg) * self.aero.cg_tail_load_sens)
    }

    fn effective_fin_arm(&self, w: &Weight) -> f32 {
        (self.control.fin_arm_m - (w.cg - self.aero.cg_ref_pct) / 100.0 * self.aero.mac_m).max(0.1)
    }

    fn vmc_thrust(&self, atm: &Atmosphere, tas_kt: f32) -> f32 {
        self.engine.takeoff_thrust(atm, atm.oat_c, tas_kt, &TOConf::default())
    }

    pub fn ground_coeffs(&self, flap: &FConf, spoilers: bool, engines_failed: u8) -> (f32, f32) {
        let a = &self.aero;
        let cl = if spoilers { (flap.cl_roll - a.cl_lift_dump).max(0.0) } else { flap.cl_roll };
        let phi = a.ground_effect(a.wing_height_m);
        let cd = a.cd0_clean + flap.cd_flap + a.cd_gear + if spoilers { a.cd_spoilers } else { 0.0 } + a.cd_windmill * engines_failed as f32 + flap.k_induced * phi * cl * cl;
        (cl, cd)
    }

    pub fn airborne_cd(&self, flap: &FConf, cl: f32, gear_down: bool, engines_failed: u8) -> f32 {
        let a = &self.aero;
        let oei = if engines_failed > 0 {
            a.cd_windmill * engines_failed as f32 + a.cd_oei_trim
        } else {
            0.0
        };
        a.cd0_clean + flap.cd_flap + if gear_down { a.cd_gear } else { 0.0 } + oei + flap.k_induced * cl * cl
    }


    pub fn pitch_trim(&self, w: &Weight, flap: &FConf) -> Result<Option<TrimSetting>, PerfError> {
        let Some(pt) = &self.pitch_trim else { return Ok(None) };

        let specific: Vec<&TrimTable> =
            pt.tables.iter().filter(|t| t.flaps.iter().any(|n| *n == flap.name)).collect();
        let set: Vec<&TrimTable> = if specific.is_empty() {
            pt.tables.iter().filter(|t| t.flaps.is_empty()).collect()
        } else {
            specific
        };
        if set.is_empty() {
            return Err(PerfError::Trim(format!("'{}' there is not any trim table to selected flap. ", flap.name)));
        }

        let mut rows = Vec::with_capacity(set.len());
        for t in &set {
            let v = interp_linear(&t.points, w.cg).ok_or_else(|| {
                PerfError::Trim(format!("CG %{:.1} MAC CG out of trim range. ", w.cg))
            })?;
            rows.push((t.weight_kg, v));
        }

        let value = if rows.len() == 1 {
            rows[0].1
        } else {
            let mut pts = Vec::with_capacity(rows.len());
            for (wt, v) in rows {
                let wt = wt.ok_or_else(|| {
                    PerfError::Trim("Weight_kg is missing. ".to_string())
                })?;
                pts.push([wt, v]);
            }
            pts.sort_by(|a, b| a[0].total_cmp(&b[0]));
            let tow = (w.tow_kg as f32).clamp(pts[0][0], pts[pts.len() - 1][0]);
            interp_linear(&pts, tow).expect("Weight compressed within the table range. ")
        };

        if pt.min.is_some_and(|lo| value < lo) || pt.max.is_some_and(|hi| value > hi) {
            return Err(PerfError::Trim(format!("Calculated trim ({value:.2}) out of green band(limits)")));
        }
        Ok(Some(TrimSetting { value, unit: pt.unit }))
    }

    // speeds

    pub fn calculate_speeds(&self, w: &Weight, atm: &Atmosphere, flap: &FConf) -> Result<VSpeeds, PerfError> {
        self.check_weight(w)?;
        let s = self.aero.wing_area_m2;
        let we = self.effective_weight_n(w);

        let vsr = (2.0 * we / (RHO_0 * s * flap.cl_max)).sqrt() * MS_TO_KT;
        let vmu = (2.0 * we / (RHO_0 * s * flap.cl_unstick)).sqrt() * MS_TO_KT;

        let c = &self.control;
        let arm = self.effective_fin_arm(w);
        let arm_ratio = arm / c.fin_arm_m;
        let yaw_authority = (c.rudder_deriv * c.rudder_max_rad + c.weathercock_deriv * c.max_sideslip_rad) * arm_ratio;

        let (mut vmcg, mut vmca) = (vsr, vsr);
        for _ in 0..8 {
            let tg = self.vmc_thrust(atm, atm.cas_to_tas(vmcg));
            vmcg = (2.0 * tg * c.engine_arm_m / (RHO_0 * c.fin_area_m2 * c.fin_cy_max * arm)).sqrt() * MS_TO_KT;

            let ta = self.vmc_thrust(atm, atm.cas_to_tas(vmca));
            vmca = (2.0 * ta * c.engine_arm_m / (RHO_0 * s * self.aero.wingspan_m * yaw_authority)).sqrt()
                * MS_TO_KT
                * flap.vmca_factor;
        }
        if !(vsr.is_finite() && vmu.is_finite() && vmcg.is_finite() && vmca.is_finite()) {
            return Err(PerfError::InvalidData(vec!["Speed calculation is not finite.".to_string()]));
        }

        Ok(VSpeeds {vsr,vmu,vmcg,vmca,vr_min: (1.05 * vmca).max(1.05 * vmu), v2_min: (1.13 * vsr).max(1.10 * vmca) })
    }

    pub fn vmbe_gs_kt(&self, w: &Weight) -> Option<f32> {
        self.gear
            .max_brake_energy_mj
            .map(|e| (2.0 * e * 1.0e6 / w.tow_kg as f32).sqrt() * MS_TO_KT)
    }

    pub fn required_gradient_pct(&self, seg: ClimbSegment) -> f32 {
        let i = match self.n_engines {
            0..=2 => 0,
            3 => 1,
            _ => 2,
        };
        match seg {
            ClimbSegment::First => [0.0, 0.3, 0.5][i],
            ClimbSegment::Second => [2.4, 2.7, 3.0][i],
        }
    }

    pub fn climb_gradient_pct(
        &self,
        w: &Weight,
        atm: &Atmosphere,
        rating_temp_c: f32,
        bleed: &TOConf,
        flap: &FConf,
        v2_cas: f32,
        seg: ClimbSegment,
    ) -> f32
    {
        let tas = atm.cas_to_tas(v2_cas);
        let qs = atm.dynamic_pressure(tas) * self.aero.wing_area_m2;
        let weight = w.tow_kg as f32 * G;
        let cl = weight / qs;
        let cd = self.airborne_cd(flap, cl, seg == ClimbSegment::First, 1);
        let thrust = self.engine.takeoff_thrust(atm, rating_temp_c, tas, bleed) * (self.n_engines - 1) as f32;
        (thrust - qs * cd) / weight * 100.0
    }
}