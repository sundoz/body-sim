//! Прогон модели во времени: ткань + лечение + почасовая история метрик.

use crate::params::Params;
use crate::report::Metrics;
use crate::sim::{self, Scratch};
use crate::therapy::Therapy;
use crate::tissue::{Tissue, WoundShape};

pub struct Simulation {
    pub p: Params,
    pub tissue: Tissue,
    pub therapy: Therapy,
    pub steps: usize,
    /// Метрики каждый модельный час, начиная с 0.
    pub history: Vec<Metrics>,
    /// Когда рана впервые эпителизировалась (часы).
    pub closed_at: Option<f32>,
    scratch: Scratch,
}

impl Simulation {
    pub fn new(p: Params, w: usize, h: usize, shape: WoundShape, depth_mm: f32) -> Self {
        let mut tissue = Tissue::healthy(w, h, &p);
        tissue.injure(shape, depth_mm, &p);
        let mut s = Self {
            scratch: Scratch::new(w * h),
            p,
            tissue,
            therapy: Therapy::default(),
            steps: 0,
            history: Vec::new(),
            closed_at: None,
        };
        s.record();
        s
    }

    pub fn hours(&self) -> f32 {
        self.steps as f32 * self.p.dt
    }

    pub fn steps_per_hour(&self) -> usize {
        (1.0 / self.p.dt).round() as usize
    }

    pub fn latest(&self) -> &Metrics {
        self.history.last().expect("история заполняется в new()")
    }

    pub fn run_steps(&mut self, n: usize) {
        let per_hour = self.steps_per_hour();
        for _ in 0..n {
            let now = self.hours();
            self.therapy.tick(&mut self.tissue, &self.p, now);
            sim::step(&mut self.tissue, &self.p, &mut self.scratch);
            self.steps += 1;
            if self.steps % per_hour == 0 {
                self.record();
            }
        }
    }

    /// Новое повреждение (например, мышью): рана снова считается открытой.
    pub fn injure_disk(&mut self, cx: f32, cy: f32, r: f32, depth_mm: f32) {
        self.tissue.injure_disk(cx, cy, r, depth_mm, &self.p);
        self.closed_at = None;
    }

    fn record(&mut self) {
        let m = Metrics::measure(&self.tissue, self.hours(), self.p.cell_mm);
        if self.closed_at.is_none() && self.steps > 0 && m.open_fraction < 0.01 {
            self.closed_at = Some(self.hours());
        }
        self.history.push(m);
    }
}
