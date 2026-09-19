//! Калибровка параметров по литературным данным: `cargo run --release --bin calibrate [итераций] [с этапа]`.
//! Этапы идут по очереди, каждый стартует с результата предыдущего.
//! В конце печатается таблица «литература / было / стало» и новые значения параметров.

use std::time::Instant;

use body_sim::calibration::{calibrate, loss, stages, Obs};
use body_sim::params::Params;

fn table(before: &[Obs], after: &[Obs]) {
    println!("  {:<48} {:>10} {:>9} {:>9}", "наблюдение", "литература", "было", "стало");
    for (b, a) in before.iter().zip(after) {
        println!("  {:<48} {:>6.2}±{:<4.2} {:>9.2} {:>9.2}", b.name, b.target, b.tol, b.value, a.value);
    }
    println!("  ошибка: {:.2} → {:.2}", loss(before), loss(after));
}

fn main() {
    let max_iter: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(60);
    let first: usize = std::env::args().nth(2).and_then(|s| s.parse().ok()).unwrap_or(1);
    let mut p = Params::default();
    let start = Instant::now();
    let mut changed = Vec::new();
    for (k, stage) in stages().iter().enumerate().skip(first.saturating_sub(1)) {
        println!("\n=== Этап {}: {} ({} параметров) ===", k + 1, stage.name, stage.knobs.len());
        let t = Instant::now();
        let r = calibrate(stage, &mut p, max_iter, |it, best| {
            if it % 5 == 0 {
                println!("  итерация {it:>3}: ошибка {best:.3} ({:.0} с)", t.elapsed().as_secs_f32());
            }
        });
        table(&r.before, &r.after);
        for (name, old, new) in r.values {
            println!("  {name:<16} {old:.5} → {new:.5}");
            changed.push((name, new));
        }
    }
    println!("\nНовые значения для Params::default ({:.0} с):", start.elapsed().as_secs_f32());
    for (name, v) in changed {
        println!("            {name}: {v:.5},");
    }
}
