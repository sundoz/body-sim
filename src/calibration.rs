//! Калибровка параметров по литературным данным.
//!
//! Каждый этап подбирает несколько параметров так, чтобы модель воспроизводила
//! набор наблюдений (пик клеток, сроки эпителизации, прочность рубца). Ошибка —
//! сумма квадратов отклонений, делённых на допуск; оптимизатор — Нелдер–Мид
//! в логарифмическом пространстве параметров с границами.
//!
//! Ориентиры:
//! - нейтрофилы преобладают первые 24–48 ч, макрофаги — на 2–3-й день,
//!   фибробласты — с 3-го дня с пиком на 1–2-й неделе (обзоры фаз заживления);
//! - поверхностная ссадина эпителизируется примерно за 1–1.5 недели, донорское место
//!   расщеплённого трансплантата — за 10–20 дней, небольшая полнослойная рана
//!   вторичным натяжением — за 2–4 недели;
//! - прочность рубца: минимальная в первую неделю, ~20% к 3-й неделе, 30–50% к 4–6-й,
//!   ~60% к полугоду, не выше ~80% нормальной кожи (Levenson и др., 1965).

use crate::params::Params;
use crate::report::Metrics;
use crate::simulation::Simulation;
use crate::tissue::WoundShape;

/// Наблюдение: что модель показала и что должно быть по литературе.
#[derive(Clone, Debug)]
pub struct Obs {
    pub name: &'static str,
    pub target: f32,
    pub tol: f32,
    pub value: f32,
}

impl Obs {
    pub fn z(&self) -> f32 {
        (self.value - self.target) / self.tol
    }
}

pub fn loss(obs: &[Obs]) -> f32 {
    obs.iter().map(|o| o.z() * o.z()).sum()
}

/// Настраиваемый параметр с границами поиска.
#[derive(Clone, Copy)]
pub struct Knob {
    pub name: &'static str,
    pub get: fn(&Params) -> f32,
    pub set: fn(&mut Params, f32),
    pub lo: f32,
    pub hi: f32,
}

impl Knob {
    /// Значение ↔ координата 0..1 (логарифмическая шкала между границами).
    fn to_unit(&self, v: f32) -> f32 {
        ((v / self.lo).ln() / (self.hi / self.lo).ln()).clamp(0.0, 1.0)
    }

    fn from_unit(&self, u: f32) -> f32 {
        self.lo * (self.hi / self.lo).powf(u.clamp(0.0, 1.0))
    }
}

macro_rules! knob {
    ($field:ident, $lo:expr, $hi:expr) => {
        Knob { name: stringify!($field), get: |p| p.$field, set: |p, v| p.$field = v, lo: $lo, hi: $hi }
    };
}

/// Этап калибровки: какие параметры двигаем и как считаем наблюдения.
pub struct Stage {
    pub name: &'static str,
    pub knobs: Vec<Knob>,
    pub observe: fn(&Params) -> Vec<Obs>,
}

pub fn stages() -> Vec<Stage> {
    vec![
        Stage {
            name: "Воспаление и пролиферация",
            knobs: vec![
                knob!(death_neut, 0.01, 0.1),
                knob!(recruit_mac, 0.01, 0.2),
                knob!(switch_m1_m2, 0.005, 0.1),
                knob!(death_m2, 0.002, 0.03),
                knob!(prolif_fib, 0.02, 0.2),
                knob!(fib_return, 0.005, 0.08),
            ],
            observe: observe_cells,
        },
        Stage {
            name: "Прочность рубца",
            knobs: vec![knob!(collagen_rate, 0.005, 0.05), knob!(maturation_rate, 0.0002, 0.005)],
            observe: observe_strength,
        },
        Stage {
            name: "Эпителизация",
            knobs: vec![
                knob!(epi_rate, 0.005, 0.4),
                knob!(d_epi, 0.005, 0.5),
                knob!(adnexal_rate, 0.00002, 0.02),
                knob!(adnexal_exponent, 1.0, 6.0),
                knob!(fill_rate, 0.01, 0.1),
            ],
            observe: observe_closure,
        },
    ]
}

/// Здоровый пациент, круглая рана Ø8 мм. Если `until_closed` — прогон обрывается после эпителизации.
fn run(p: &Params, depth_mm: f32, days: f32, until_closed: bool) -> Simulation {
    let r = 4.0 / p.cell_mm;
    let mut sim = Simulation::new(p.clone(), 96, 48, WoundShape::Circle { radius: r }, depth_mm);
    let per_day = 24 * sim.steps_per_hour();
    for _ in 0..days.ceil() as usize {
        sim.run_steps(per_day);
        if until_closed && sim.closed_at.is_some() {
            break;
        }
    }
    sim
}

fn peak_day(sim: &Simulation, f: impl Fn(&Metrics) -> f32) -> f32 {
    sim.history.iter().max_by(|a, b| f(a).total_cmp(&f(b))).map_or(0.0, |m| m.hours / 24.0)
}

fn at_day(sim: &Simulation, day: f32) -> &Metrics {
    let i = ((day * 24.0).round() as usize).min(sim.history.len() - 1);
    &sim.history[i]
}

fn observe_cells(p: &Params) -> Vec<Obs> {
    let s = run(p, 2.5, 16.0, false);
    let neut_peak = s.history.iter().map(|m| m.neutrophils).fold(0.0, f32::max).max(1e-6);
    vec![
        Obs { name: "пик нейтрофилов, день", target: 1.0, tol: 0.5, value: peak_day(&s, |m| m.neutrophils) },
        Obs { name: "пик макрофагов, день", target: 2.5, tol: 0.75, value: peak_day(&s, |m| m.m1 + m.m2) },
        Obs {
            name: "нейтрофилы на 7-й день / пик",
            target: 0.1,
            tol: 0.1,
            value: at_day(&s, 7.0).neutrophils / neut_peak,
        },
        Obs { name: "пик фибробластов, день", target: 10.0, tol: 3.0, value: peak_day(&s, |m| m.fibroblasts) },
    ]
}

/// День эпителизации раны глубиной `depth_mm`; если не закрылась за `limit` дней — штраф.
fn closure_day(p: &Params, depth_mm: f32, limit: f32) -> f32 {
    let s = run(p, depth_mm, limit, true);
    s.closed_at.map_or(limit + 10.0, |h| h / 24.0)
}

fn observe_closure(p: &Params) -> Vec<Obs> {
    vec![
        Obs { name: "эпителизация ссадины 0.1 мм, день", target: 8.0, tol: 3.0, value: closure_day(p, 0.1, 25.0) },
        Obs {
            name: "эпителизация донорского места 0.3 мм, день",
            target: 12.0,
            tol: 4.0,
            value: closure_day(p, 0.3, 30.0),
        },
        Obs {
            name: "эпителизация глубокой дермальной раны 1 мм, день",
            target: 21.0,
            tol: 6.0,
            value: closure_day(p, 1.0, 45.0),
        },
        Obs {
            name: "эпителизация полнослойной раны 2.5 мм, день",
            target: 24.0,
            tol: 6.0,
            value: closure_day(p, 2.5, 50.0),
        },
    ]
}

fn observe_strength(p: &Params) -> Vec<Obs> {
    let s = run(p, 2.5, 90.0, false);
    vec![
        Obs { name: "прочность на 7-й день", target: 0.03, tol: 0.05, value: at_day(&s, 7.0).strength },
        Obs { name: "прочность на 21-й день", target: 0.2, tol: 0.07, value: at_day(&s, 21.0).strength },
        Obs { name: "прочность на 42-й день", target: 0.4, tol: 0.1, value: at_day(&s, 42.0).strength },
        Obs { name: "прочность на 90-й день", target: 0.62, tol: 0.08, value: at_day(&s, 90.0).strength },
    ]
}

/// Минимизация Нелдера–Мида. `f` — функция от точки, `x0` — старт, `step` — размер начального симплекса.
/// `on_iter` получает номер итерации и лучшее значение.
pub fn nelder_mead(
    mut f: impl FnMut(&[f32]) -> f32,
    x0: &[f32],
    step: f32,
    max_iter: usize,
    tol: f32,
    mut on_iter: impl FnMut(usize, f32),
) -> (Vec<f32>, f32) {
    let n = x0.len();
    let mut simplex: Vec<(Vec<f32>, f32)> = Vec::with_capacity(n + 1);
    simplex.push((x0.to_vec(), f(x0)));
    for i in 0..n {
        let mut x = x0.to_vec();
        x[i] += if x[i] + step <= 1.0 { step } else { -step };
        let fx = f(&x);
        simplex.push((x, fx));
    }
    for it in 0..max_iter {
        simplex.sort_by(|a, b| a.1.total_cmp(&b.1));
        on_iter(it, simplex[0].1);
        if (simplex[n].1 - simplex[0].1).abs() <= tol * (1.0 + simplex[0].1.abs()) {
            break;
        }
        let centroid: Vec<f32> = (0..n).map(|j| simplex[..n].iter().map(|(x, _)| x[j]).sum::<f32>() / n as f32).collect();
        let along = |t: f32| -> Vec<f32> { (0..n).map(|j| centroid[j] + t * (simplex[n].0[j] - centroid[j])).collect() };
        let xr = along(-1.0);
        let fr = f(&xr);
        if fr < simplex[0].1 {
            let xe = along(-2.0);
            let fe = f(&xe);
            simplex[n] = if fe < fr { (xe, fe) } else { (xr, fr) };
        } else if fr < simplex[n - 1].1 {
            simplex[n] = (xr, fr);
        } else {
            let (xc, fc) = if fr < simplex[n].1 {
                let x = along(-0.5);
                let v = f(&x);
                (x, v)
            } else {
                let x = along(0.5);
                let v = f(&x);
                (x, v)
            };
            if fc < simplex[n].1.min(fr) {
                simplex[n] = (xc, fc);
            } else {
                // Сжатие всего симплекса к лучшей точке.
                let best = simplex[0].0.clone();
                for k in 1..=n {
                    let x: Vec<f32> = (0..n).map(|j| best[j] + 0.5 * (simplex[k].0[j] - best[j])).collect();
                    let fx = f(&x);
                    simplex[k] = (x, fx);
                }
            }
        }
    }
    simplex.sort_by(|a, b| a.1.total_cmp(&b.1));
    simplex.swap_remove(0)
}

/// Результат этапа: наблюдения до и после, найденные значения параметров.
pub struct StageResult {
    pub before: Vec<Obs>,
    pub after: Vec<Obs>,
    pub values: Vec<(&'static str, f32, f32)>,
}

/// Откалибровать один этап, меняя `p` на месте.
pub fn calibrate(stage: &Stage, p: &mut Params, max_iter: usize, mut progress: impl FnMut(usize, f32)) -> StageResult {
    let before = (stage.observe)(p);
    let old: Vec<f32> = stage.knobs.iter().map(|k| (k.get)(p)).collect();
    let x0: Vec<f32> = stage.knobs.iter().map(|k| k.to_unit((k.get)(p))).collect();
    let base = p.clone();
    let apply = |x: &[f32], p: &mut Params| {
        for (k, &u) in stage.knobs.iter().zip(x) {
            (k.set)(p, k.from_unit(u));
        }
    };
    let (best, _) = nelder_mead(
        |x| {
            let mut q = base.clone();
            apply(x, &mut q);
            loss(&(stage.observe)(&q))
        },
        &x0,
        0.12,
        max_iter,
        1e-3,
        &mut progress,
    );
    apply(&best, p);
    let after = (stage.observe)(p);
    let values = stage.knobs.iter().zip(old).map(|(k, o)| (k.name, o, (k.get)(p))).collect();
    StageResult { before, after, values }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nelder_mead_finds_minimum_of_a_bowl() {
        let (x, fx) = nelder_mead(
            |x| (x[0] - 0.3).powi(2) + 2.0 * (x[1] - 0.7).powi(2),
            &[0.5, 0.5],
            0.1,
            500,
            1e-9,
            |_, _| {},
        );
        assert!(fx < 1e-6, "{fx}");
        assert!((x[0] - 0.3).abs() < 1e-2 && (x[1] - 0.7).abs() < 1e-2, "{x:?}");
    }

    #[test]
    fn nelder_mead_handles_rosenbrock() {
        let (x, _) = nelder_mead(
            |x| {
                let (a, b) = (x[0] * 2.0, x[1] * 2.0);
                (1.0 - a).powi(2) + 100.0 * (b - a * a).powi(2)
            },
            &[0.1, 0.1],
            0.1,
            3000,
            1e-12,
            |_, _| {},
        );
        assert!((x[0] - 0.5).abs() < 0.05 && (x[1] - 0.5).abs() < 0.05, "{x:?}");
    }

    #[test]
    fn loss_is_zero_on_target_and_scaled_by_tolerance() {
        let o = |value| Obs { name: "x", target: 10.0, tol: 2.0, value };
        assert_eq!(loss(&[o(10.0)]), 0.0);
        assert!((loss(&[o(14.0)]) - 4.0).abs() < 1e-6, "два допуска → 4");
    }

    #[test]
    fn knob_log_mapping_round_trips() {
        let k = knob!(epi_rate, 0.01, 1.0);
        for v in [0.01, 0.05, 0.1, 0.5, 1.0] {
            assert!((k.from_unit(k.to_unit(v)) - v).abs() / v < 1e-4);
        }
        assert!((k.from_unit(0.5) - 0.1).abs() < 1e-5, "середина — среднее геометрическое");
        let mut p = Params::default();
        (k.set)(&mut p, 0.2);
        assert_eq!((k.get)(&p), 0.2);
    }

    #[test]
    fn every_stage_knob_starts_inside_its_bounds() {
        let p = Params::default();
        for s in stages() {
            for k in &s.knobs {
                let v = (k.get)(&p);
                assert!(k.lo <= v && v <= k.hi, "{}: {v} вне [{}, {}]", k.name, k.lo, k.hi);
            }
        }
    }
}
