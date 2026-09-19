//! Простой покадровый профайлер: `--bench N` прогоняет N кадров и печатает среднее и максимум по участкам.

use std::time::Instant;

pub struct Prof {
    names: Vec<&'static str>,
    sums: Vec<f64>,
    maxs: Vec<f64>,
    frame: Vec<f64>,
    frames: u32,
}

impl Prof {
    pub fn new() -> Self {
        Self { names: Vec::new(), sums: Vec::new(), maxs: Vec::new(), frame: Vec::new(), frames: 0 }
    }

    fn slot(&mut self, name: &'static str) -> usize {
        match self.names.iter().position(|n| *n == name) {
            Some(i) => i,
            None => {
                self.names.push(name);
                self.sums.push(0.0);
                self.maxs.push(0.0);
                self.frame.push(0.0);
                self.names.len() - 1
            }
        }
    }

    /// Добавить к участку `name` время в миллисекундах.
    pub fn add(&mut self, name: &'static str, ms: f64) {
        let i = self.slot(name);
        self.frame[i] += ms;
    }

    pub fn time<R>(&mut self, name: &'static str, f: impl FnOnce() -> R) -> R {
        let t = Instant::now();
        let r = f();
        self.add(name, t.elapsed().as_secs_f64() * 1000.0);
        r
    }

    pub fn end_frame(&mut self) {
        for i in 0..self.names.len() {
            self.sums[i] += self.frame[i];
            self.maxs[i] = self.maxs[i].max(self.frame[i]);
            self.frame[i] = 0.0;
        }
        self.frames += 1;
    }

    /// Сбросить накопленное (после прогрева: первые кадры растеризуют шрифты и грузят текстуры).
    pub fn reset(&mut self) {
        *self = Self::new();
    }

    pub fn report(&self) -> String {
        let n = self.frames.max(1) as f64;
        let mut s = format!("{:<28} {:>10} {:>10}\n", "участок", "среднее мс", "макс мс");
        for i in 0..self.names.len() {
            s += &format!("{:<28} {:>10.2} {:>10.2}\n", self.names[i], self.sums[i] / n, self.maxs[i]);
        }
        s += &format!("кадров: {}", self.frames);
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn averages_and_maxima() {
        let mut p = Prof::new();
        p.add("a", 2.0);
        p.add("a", 1.0);
        p.add("b", 5.0);
        p.end_frame();
        p.add("a", 1.0);
        p.end_frame();
        assert_eq!(p.frames, 2);
        let a = p.names.iter().position(|n| *n == "a").unwrap();
        assert!((p.sums[a] / 2.0 - 2.0).abs() < 1e-9, "среднее (3 + 1) / 2");
        assert!((p.maxs[a] - 3.0).abs() < 1e-9);
        let r = p.report();
        assert!(r.contains("a") && r.contains("b") && r.contains("кадров: 2"));
        let v = p.time("c", || 42);
        assert_eq!(v, 42);
        p.reset();
        assert_eq!(p.frames, 0);
        assert!(p.names.is_empty());
    }
}
