//! Прогон модели во времени: ткань + лечение + почасовая история метрик.

use std::sync::atomic::{AtomicU64, Ordering};

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
    /// Меняется при любом изменении ткани; уникальна между симуляциями —
    /// по ней визуализация понимает, что пора перерисоваться.
    pub revision: u64,
    scratch: Scratch,
}

static REVISION: AtomicU64 = AtomicU64::new(1);

fn next_revision() -> u64 {
    REVISION.fetch_add(1, Ordering::Relaxed)
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
            revision: next_revision(),
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

    /// Отметить, что ткань изменилась снаружи (например, процедурой).
    pub fn touch(&mut self) {
        self.revision = next_revision();
    }

    pub fn run_steps(&mut self, n: usize) {
        if n > 0 {
            self.touch();
        }
        let per_hour = self.steps_per_hour();
        for _ in 0..n {
            let now = self.hours();
            self.therapy.tick(&mut self.tissue, &self.p, now);
            sim::step(&mut self.tissue, &self.p, &mut self.scratch);
            self.steps += 1;
            if self.steps.is_multiple_of(per_hour) {
                self.record();
            }
        }
    }

    /// Новое повреждение (например, мышью): рана снова считается открытой.
    pub fn injure_disk(&mut self, cx: f32, cy: f32, r: f32, depth_mm: f32) {
        self.tissue.injure_disk(cx, cy, r, depth_mm, &self.p);
        self.touch();
        self.closed_at = None;
    }

    fn record(&mut self) {
        let mut m = Metrics::measure(&self.tissue, &self.p, self.hours());
        m.initial_wound_mm2 = self.history.first().map_or(m.wound_mm2, |f| f.initial_wound_mm2);
        if self.closed_at.is_none() && self.steps > 0 && m.open_fraction < 0.01 {
            self.closed_at = Some(self.hours());
        }
        self.history.push(m);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::Scenario;

    fn sim(depth: f32) -> Simulation {
        let p = Scenario::Healthy.params();
        Simulation::new(p, 96, 48, WoundShape::Circle { radius: 16.0 }, depth)
    }

    #[test]
    fn history_is_hourly() {
        let mut s = sim(2.5);
        assert_eq!(s.history.len(), 1);
        s.run_steps(25 * s.steps_per_hour());
        assert_eq!(s.history.len(), 26);
        for (k, m) in s.history.iter().enumerate() {
            assert!((m.hours - k as f32).abs() < 1e-3);
        }
        assert!((s.hours() - 25.0).abs() < 1e-3);
    }

    #[test]
    fn revision_tracks_every_change() {
        let mut s = sim(2.5);
        let r0 = s.revision;
        s.run_steps(0);
        assert_eq!(s.revision, r0, "без шагов ткань не менялась");
        s.run_steps(1);
        let r1 = s.revision;
        assert_ne!(r1, r0);
        s.injure_disk(10.0, 10.0, 3.0, 1.0);
        assert_ne!(s.revision, r1);
        let r2 = s.revision;
        s.touch();
        assert_ne!(s.revision, r2);
        assert_ne!(sim(2.5).revision, sim(2.5).revision, "у разных симуляций разные ревизии");
    }

    #[test]
    fn superficial_wound_closes_within_two_weeks_and_new_injury_reopens_it() {
        let mut s = sim(0.1);
        s.run_steps(14 * 24 * s.steps_per_hour());
        assert!(s.closed_at.is_some());
        s.injure_disk(48.0, 24.0, 4.0, 1.0);
        assert!(s.closed_at.is_none());
    }
}
