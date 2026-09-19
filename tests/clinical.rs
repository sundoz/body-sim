//! Клиническая регрессия: модель должна воспроизводить известные исходы.
//! Пороги намеренно широкие — проверяется качественное поведение, а не точные числа.

use body_sim::params::Scenario;
use body_sim::report::{Condition, Metrics};
use body_sim::simulation::Simulation;
use body_sim::therapy::{Antibiotic, Antiseptic, Debrider};
use body_sim::tissue::WoundShape;

const W: usize = 96;
const H: usize = 48;

/// Круглая рана Ø8 мм заданной глубины.
fn wound(scenario: Scenario, depth_mm: f32) -> Simulation {
    let p = scenario.params();
    let r = 4.0 / p.cell_mm;
    Simulation::new(p, W, H, WoundShape::Circle { radius: r }, depth_mm)
}

#[derive(Default)]
struct Plan {
    antiseptic: Option<(Antiseptic, f32)>,
    antibiotics: Vec<(Antibiotic, f32)>,
    debriders: Vec<Debrider>,
    surgery_days: Vec<f32>,
    start_day: f32,
}

fn run(mut sim: Simulation, days: f32, plan: &Plan) -> Simulation {
    let per_hour = sim.steps_per_hour();
    let mut started = false;
    for hour in 0..(days * 24.0).round() as usize {
        let now = hour as f32;
        if !started && now >= plan.start_day * 24.0 {
            started = true;
            if let Some((a, every)) = plan.antiseptic {
                sim.therapy.start_antiseptic(a, every, now);
            }
            for (drug, every) in &plan.antibiotics {
                sim.therapy.start_antibiotic(*drug, *every, now);
            }
            for d in &plan.debriders {
                sim.tissue.debriders[d.index()] = true;
            }
        }
        if plan.surgery_days.iter().any(|d| (d * 24.0).round() as usize == hour) {
            sim.therapy.debride_now(&mut sim.tissue, &sim.p, now);
        }
        sim.run_steps(per_hour);
    }
    sim
}

fn untreated(scenario: Scenario, depth_mm: f32, days: f32) -> Simulation {
    run(wound(scenario, depth_mm), days, &Plan::default())
}

fn closed_day(sim: &Simulation) -> Option<f32> {
    sim.closed_at.map(|h| h / 24.0)
}

fn worst(sim: &Simulation) -> Condition {
    sim.history.iter().map(Metrics::condition).fold(Condition::Healed, |a, b| if b > a { b } else { a })
}

fn last(sim: &Simulation) -> &Metrics {
    sim.latest()
}

// ---------------------------------------------------------------- заживление и глубина

#[test]
fn healthy_full_thickness_wound_heals_in_three_to_five_weeks() {
    // Полнослойная рана Ø8 мм вторичным натяжением: литература — 2–5 недель.
    let s = untreated(Scenario::Healthy, 2.5, 40.0);
    let day = closed_day(&s).expect("здоровая рана должна зажить");
    assert!((18.0..36.0).contains(&day), "эпителизация на {day:.1} день");
    assert!(worst(&s) <= Condition::Healing, "без осложнений");
    assert_eq!(last(&s).bacteria, 0.0);
    let peak_neut_day = s.history.iter().max_by(|a, b| a.neutrophils.total_cmp(&b.neutrophils)).unwrap().hours / 24.0;
    assert!(peak_neut_day < 2.5, "пик нейтрофилов в первые сутки-двое, а не на {peak_neut_day:.1}");
}

#[test]
fn deeper_wounds_heal_slower() {
    let days: Vec<f32> = [0.1, 1.0, 2.5, 5.0]
        .iter()
        .map(|&d| closed_day(&untreated(Scenario::Healthy, d, 45.0)).expect("должна зажить"))
        .collect();
    assert!((5.0..13.0).contains(&days[0]), "ссадина — за 1–1.5 недели: {days:?}");
    assert!(days.windows(2).all(|w| w[0] < w[1]), "чем глубже, тем дольше: {days:?}");
}

#[test]
fn superficial_wound_heals_without_scar_but_deep_one_scars() {
    let shallow = untreated(Scenario::Healthy, 0.1, 20.0);
    assert!(last(&shallow).strength > 0.95, "ссадина без рубца");
    let deep = untreated(Scenario::Healthy, 2.5, 60.0);
    let s = last(&deep).strength;
    assert!((0.25..0.6).contains(&s), "полнослойный рубец слабее кожи: {s}");
}

#[test]
fn comorbidities_slow_healing() {
    let healthy = closed_day(&untreated(Scenario::Healthy, 2.5, 60.0)).unwrap();
    let elderly = closed_day(&untreated(Scenario::Elderly, 2.5, 60.0)).unwrap();
    let diabetic = closed_day(&untreated(Scenario::Diabetic, 2.5, 60.0)).unwrap_or(f32::INFINITY);
    assert!(elderly > healthy && diabetic > elderly, "{healthy} < {elderly} < {diabetic}");
}

// ---------------------------------------------------------------- глубокие слои

#[test]
fn small_muscle_loss_regenerates() {
    let s = untreated(Scenario::Healthy, 9.0, 21.0);
    let m = last(&s);
    assert!((m.lost_muscle_mm - 0.6).abs() < 1e-3);
    assert!(m.muscle_regen > 0.6, "к 3-й неделе большая часть волокон: {}", m.muscle_regen);
    assert!(m.muscle_fibrosis < 0.3, "фиброза немного: {}", m.muscle_fibrosis);
    let day7 = &s.history[7 * 24];
    assert!((0.3..0.7).contains(&day7.muscle_regen), "к 7-му дню — около половины: {}", day7.muscle_regen);
}

#[test]
fn large_muscle_loss_heals_with_fibrosis() {
    let s = untreated(Scenario::Healthy, 12.0, 42.0);
    let m = last(&s);
    assert!(m.lost_muscle_mm > 3.0);
    assert!(m.muscle_fibrosis > m.muscle_regen, "волокна {} против фиброза {}", m.muscle_regen, m.muscle_fibrosis);
}

#[test]
fn fat_barely_regenerates_and_leaves_scar() {
    let s = untreated(Scenario::Healthy, 5.0, 60.0);
    assert!(s.closed_at.is_some());
    let m = last(&s);
    assert!((m.lost_fat_mm - 2.9).abs() < 1e-3);
    assert!(m.fat_regen < 0.15, "жир почти не возвращается: {}", m.fat_regen);
}

// ---------------------------------------------------------------- некроз и инфекция

#[test]
fn ischemia_causes_dry_necrosis_before_any_infection() {
    // Ишемия убивает ткань сама по себе; инфекция если и приходит, то позже, на мёртвую ткань.
    let s = untreated(Scenario::Ischemic, 2.5, 20.0);
    let necrosis = first_day(&s, Condition::Necrosis).expect("ишемия должна дать некроз");
    assert!(necrosis < 5.0, "некроз на {necrosis:.1} день");
    let infection = s.history.iter().find(|m| m.infected_mm2 > 1.0).map_or(f32::INFINITY, |m| m.hours / 24.0);
    assert!(necrosis < infection, "некроз {necrosis:.1} раньше инфекции {infection:.1}");
}

#[test]
fn untreated_necrotizing_infection_becomes_septic() {
    let s = untreated(Scenario::Necrotizing, 2.5, 15.0);
    assert_eq!(worst(&s), Condition::Sepsis);
}

#[test]
fn antibiotics_alone_do_not_cure_necrotizing_infection() {
    let plan = Plan { antibiotics: vec![(Antibiotic::Cefazolin, 8.0)], start_day: 0.5, ..Default::default() };
    let s = run(wound(Scenario::Necrotizing, 2.5), 30.0, &plan);
    assert!(s.closed_at.is_none());
    assert!(last(&s).condition() >= Condition::Necrosis);
}

#[test]
fn clindamycin_added_to_beta_lactam_limits_necrosis() {
    let cef = Plan { antibiotics: vec![(Antibiotic::Cefazolin, 8.0)], start_day: 0.5, ..Default::default() };
    let both = Plan {
        antibiotics: vec![(Antibiotic::Cefazolin, 8.0), (Antibiotic::Clindamycin, 8.0)],
        start_day: 0.5,
        ..Default::default()
    };
    let a = run(wound(Scenario::Necrotizing, 2.5), 30.0, &cef);
    let b = run(wound(Scenario::Necrotizing, 2.5), 30.0, &both);
    assert!(last(&b).open_mm2 < last(&a).open_mm2, "{} < {}", last(&b).open_mm2, last(&a).open_mm2);
}

#[test]
fn early_surgery_with_antibiotics_cures_necrotizing_infection() {
    let plan = Plan {
        antibiotics: vec![(Antibiotic::Cefazolin, 8.0), (Antibiotic::Clindamycin, 8.0)],
        surgery_days: vec![0.5, 1.5, 3.0],
        start_day: 0.5,
        ..Default::default()
    };
    let s = run(wound(Scenario::Necrotizing, 2.5), 35.0, &plan);
    let day = closed_day(&s).expect("должна зажить");
    assert!(day < 30.0, "{day}");
    assert_eq!(last(&s).condition(), Condition::Healed);
    assert_eq!(last(&s).bacteria, 0.0);
}

#[test]
fn late_surgery_fails() {
    let plan = Plan {
        antibiotics: vec![(Antibiotic::Cefazolin, 8.0)],
        surgery_days: vec![3.0, 4.0],
        start_day: 0.5,
        ..Default::default()
    };
    let s = run(wound(Scenario::Necrotizing, 2.5), 35.0, &plan);
    assert!(s.closed_at.is_none(), "опоздавшая хирургия не спасает");
}

// ---------------------------------------------------------------- устойчивость

#[test]
fn long_beta_lactam_course_without_source_control_selects_resistance() {
    let plan = Plan { antibiotics: vec![(Antibiotic::Cefazolin, 8.0)], ..Default::default() };
    let s = run(wound(Scenario::DiabeticFoot, 2.5), 42.0, &plan);
    assert!(last(&s).resistant_fraction > 0.1, "{}", last(&s).resistant_fraction);
    let untreated = untreated(Scenario::DiabeticFoot, 2.5, 42.0);
    assert!(last(&untreated).resistant_fraction < 0.01, "без давления антибиотика устойчивые не отбираются");
}

#[test]
fn vancomycin_kills_the_resistant_strain_better_than_cefazolin() {
    let resistant_left = |drug: Antibiotic, every: f32| {
        let mut sim = wound(Scenario::Infected, 2.5);
        for i in 0..sim.tissue.len() {
            if sim.tissue.wound_mask[i] {
                sim.tissue.bacteria_res.data[i] = sim.tissue.bacteria.data[i];
                sim.tissue.bacteria.data[i] = 0.0;
            }
        }
        let plan = Plan { antibiotics: vec![(drug, every)], ..Default::default() };
        let s = run(sim, 4.0, &plan);
        s.tissue.bacteria_res.data.iter().sum::<f32>()
    };
    let cef = resistant_left(Antibiotic::Cefazolin, 8.0);
    let van = resistant_left(Antibiotic::Vancomycin, 12.0);
    assert!(van < cef, "ванкомицин {van} против цефазолина {cef}");
}

// ---------------------------------------------------------------- местное лечение

#[test]
fn gentle_antiseptics_heal_faster_than_cytotoxic_ones() {
    let heal = |a: Antiseptic, every: f32| {
        let plan = Plan { antiseptic: Some((a, every)), ..Default::default() };
        closed_day(&run(wound(Scenario::Infected, 2.5), 70.0, &plan)).unwrap_or(f32::INFINITY)
    };
    let octenidine = heal(Antiseptic::Octenidine, 24.0);
    let chlorhexidine = heal(Antiseptic::Chlorhexidine, 24.0);
    let untreated = closed_day(&untreated(Scenario::Infected, 2.5, 70.0)).unwrap_or(f32::INFINITY);
    assert!(octenidine < chlorhexidine, "{octenidine} < {chlorhexidine}");
    assert!(octenidine < untreated, "антисептик помогает инфицированной ране: {octenidine} < {untreated}");
}

#[test]
fn too_frequent_cytotoxic_antiseptic_delays_healing() {
    let heal = |every: f32| {
        let plan = Plan { antiseptic: Some((Antiseptic::PovidoneIodine, every)), ..Default::default() };
        closed_day(&run(wound(Scenario::Infected, 2.5), 60.0, &plan)).unwrap_or(f32::INFINITY)
    };
    let (daily, twice) = (heal(24.0), heal(12.0));
    assert!(daily.is_finite(), "раз в сутки рана заживает: {daily}");
    assert!(twice > daily, "{twice} > {daily}");
}

#[test]
fn surgical_debridement_clears_diabetic_foot_infection() {
    let plan = Plan {
        antiseptic: Some((Antiseptic::Octenidine, 24.0)),
        antibiotics: vec![(Antibiotic::Cefazolin, 8.0)],
        surgery_days: vec![2.0, 9.0],
        ..Default::default()
    };
    let s = run(wound(Scenario::DiabeticFoot, 2.5), 42.0, &plan);
    let m = last(&s);
    assert_eq!(m.bacteria, 0.0);
    assert_eq!(m.necrotic_mm2, 0.0);
    assert!(m.condition() < Condition::Necrosis);
}

#[test]
fn simulation_is_reproducible() {
    let a = untreated(Scenario::Infected, 2.5, 3.0);
    let b = untreated(Scenario::Infected, 2.5, 3.0);
    assert_eq!(a.tissue.bacteria.data, b.tissue.bacteria.data);
    assert_eq!(a.history.len(), b.history.len());
    assert_eq!(last(&a).csv_row(), last(&b).csv_row());
}

fn first_day(sim: &Simulation, cond: Condition) -> Option<f32> {
    sim.history.iter().find(|m| m.condition() >= cond).map(|m| m.hours / 24.0)
}

#[test]
fn necrotizing_infection_spreads_within_days_and_clindamycin_delays_it() {
    let untreated = untreated(Scenario::Necrotizing, 2.5, 12.0);
    let spread = first_day(&untreated, Condition::Spreading).expect("должна распространиться");
    assert!(spread < 4.0, "распространение на {spread:.1} день");
    let sepsis = first_day(&untreated, Condition::Sepsis).expect("должен развиться сепсис");
    assert!(sepsis < 10.0, "сепсис на {sepsis:.1} день");

    let plan = Plan {
        antibiotics: vec![(Antibiotic::Cefazolin, 8.0), (Antibiotic::Clindamycin, 8.0)],
        start_day: 0.5,
        ..Default::default()
    };
    let treated = run(wound(Scenario::Necrotizing, 2.5), 12.0, &plan);
    let delayed = first_day(&treated, Condition::Spreading).unwrap_or(f32::INFINITY);
    assert!(delayed > spread + 2.0, "клиндамицин тормозит токсины: {delayed:.1} против {spread:.1}");
}
