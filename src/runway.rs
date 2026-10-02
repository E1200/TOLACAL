use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Deserialize, Serialize)]
pub struct Wind {
    pub dir_deg: f32,
    pub speed_kt: f32,
}

impl Wind {
    pub fn components(&self, rwy_heading_deg: f32) -> (f32, f32) {
        let a = (self.dir_deg - rwy_heading_deg).to_radians();
        (self.speed_kt * a.cos(), self.speed_kt * a.sin())
    }
}

pub fn factored_headwind(hw_kt: f32) -> f32 {
    if hw_kt >= 0.0 { 0.5 * hw_kt } else { 1.5 * hw_kt }
}

#[derive(Debug, Clone, Copy, Deserialize, Serialize)]
pub struct Runway {
    pub heading_deg: f32,
    pub tora_m: f32,
    pub swy_m: f32,
    pub cwy_m: f32,
    pub elevation_ft: f32, 
    pub slope_pct: f32,    // up
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct RunwayLimits {
    pub tora_m: f32,
    pub toda_m: f32,
    pub asda_m: f32,
    pub slope: f32,          
    pub headwind_kt: f32,    
    pub headwind_raw_kt: f32, // limit check
    pub crosswind_kt: f32,   // limit check
}

impl Runway {
    pub fn asda_m(&self) -> f32 { self.tora_m + self.swy_m }
    pub fn toda_m(&self) -> f32 { self.tora_m + self.cwy_m }

    pub fn limits(&self, wind: &Wind) -> RunwayLimits {
        let (hw, xw) = wind.components(self.heading_deg);
        RunwayLimits {
            tora_m: self.tora_m,
            toda_m: self.toda_m(),
            asda_m: self.asda_m(),
            slope: self.slope_pct / 100.0,
            headwind_kt: factored_headwind(hw),
            headwind_raw_kt: hw,
            crosswind_kt: xw,
        }
    }
}
