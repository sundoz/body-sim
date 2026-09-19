//! Параметры по умолчанию должны воспроизводить литературные ориентиры калибровки.
//! Если тест упал после правки модели — перезапустите `cargo run --release --bin calibrate`
//! и обновите значения в `Params::default`.

use body_sim::calibration::stages;
use body_sim::params::Params;

#[test]
fn default_parameters_match_literature() {
    let p = Params::default();
    let mut report = String::new();
    let mut failed = false;
    for stage in stages() {
        for o in (stage.observe)(&p) {
            let ok = o.z().abs() <= 1.25;
            failed |= !ok;
            report += &format!(
                "{} {}: {:.2} (литература {:.2} ± {:.2})\n",
                if ok { "  " } else { "✗ " },
                o.name,
                o.value,
                o.target,
                o.tol
            );
        }
    }
    assert!(!failed, "модель ушла от литературы:\n{report}");
}
