//! График показателей раны во времени с полосой фаз и отметками процедур.

use body_sim::report::{Condition, Metrics, Phase};
use body_sim::therapy::{Event, EventKind};
use macroquad::prelude::*;

use crate::ui::{self, Ui};

pub const N: usize = 14;

pub const NAMES: [&str; N] = [
    "Открытая рана",
    "Глубина",
    "Некроз",
    "Бактерии",
    "Устойчивые",
    "Биоплёнка",
    "Нейтрофилы",
    "Макрофаги M1",
    "Макрофаги M2",
    "Фибробласты",
    "Сосуды",
    "Коллаген",
    "Прочность",
    "Антибиотик",
];

pub const COLORS: [Color; N] = [
    Color::new(0.93, 0.93, 0.95, 1.0),
    Color::new(0.55, 0.62, 0.72, 1.0),
    Color::new(0.62, 0.47, 0.30, 1.0),
    Color::new(0.62, 0.84, 0.30, 1.0),
    Color::new(0.95, 0.30, 0.55, 1.0),
    Color::new(0.45, 0.62, 0.36, 1.0),
    Color::new(1.00, 0.76, 0.22, 1.0),
    Color::new(0.96, 0.52, 0.78, 1.0),
    Color::new(0.63, 0.52, 1.00, 1.0),
    Color::new(0.30, 0.80, 0.93, 1.0),
    Color::new(0.94, 0.33, 0.33, 1.0),
    Color::new(0.86, 0.72, 0.56, 1.0),
    Color::new(0.30, 0.58, 1.00, 1.0),
    Color::new(1.00, 0.60, 0.20, 1.0),
];

/// Какие серии видны по умолчанию.
pub const DEFAULT_VISIBLE: [bool; N] =
    [true, true, true, true, true, false, true, false, true, true, false, true, true, false];

/// Значение серии, нормированное на 0..1, и подпись реального значения.
pub fn value(i: usize, m: &Metrics) -> f32 {
    match i {
        0 => m.open_fraction,
        1 => m.depth / 5.0,
        2 => m.necrotic_mm2 / m.wound_mm2.max(1.0),
        3 => m.bacteria,
        4 => m.resistant_fraction,
        5 => m.biofilm,
        6 => m.neutrophils,
        7 => m.m1,
        8 => m.m2,
        9 => m.fibroblasts,
        10 => m.vessels,
        11 => m.collagen,
        12 => m.strength,
        _ => m.abx_plasma / 10.0,
    }
}

fn value_label(i: usize, m: &Metrics) -> String {
    match i {
        0 => format!("{:.0}%", m.open_fraction * 100.0),
        1 => format!("{:.1} мм", m.depth),
        2 => format!("{:.1} мм²", m.necrotic_mm2),
        4 => format!("{:.0}%", m.resistant_fraction * 100.0),
        12 => format!("{:.0}%", m.strength * 100.0),
        13 => format!("{:.1} МПК", m.abx_plasma),
        _ => format!("{:.2}", value(i, m)),
    }
}

pub fn phase_color(p: Phase) -> Color {
    match p {
        Phase::Hemostasis => Color::new(0.86, 0.24, 0.30, 1.0),
        Phase::Inflammation => Color::new(0.96, 0.58, 0.22, 1.0),
        Phase::Proliferation => Color::new(0.30, 0.78, 0.50, 1.0),
        Phase::Remodeling => Color::new(0.38, 0.58, 1.00, 1.0),
    }
}

pub fn condition_color(c: Condition) -> Color {
    match c {
        Condition::Healed => Color::new(0.30, 0.78, 0.50, 1.0),
        Condition::Healing => Color::new(0.55, 0.75, 0.95, 1.0),
        Condition::Chronic => Color::new(0.95, 0.72, 0.25, 1.0),
        Condition::Necrosis => Color::new(0.72, 0.52, 0.32, 1.0),
        Condition::Spreading => Color::new(0.95, 0.35, 0.30, 1.0),
        Condition::Sepsis => Color::new(0.80, 0.10, 0.20, 1.0),
    }
}

pub fn event_color(k: EventKind) -> Color {
    match k {
        EventKind::Antiseptic => Color::new(0.35, 0.80, 0.95, 1.0),
        EventKind::Dose(_) => Color::new(1.00, 0.60, 0.20, 1.0),
        EventKind::Debridement => Color::new(1.0, 1.0, 1.0, 1.0),
    }
}

pub fn draw(ui: &Ui, r: Rect, hist: &[Metrics], events: &[Event], visible: &mut [bool; N], now_h: f32) {
    ui::panel(r);
    ui.bold("Динамика в области раны", r.x + 16.0, r.y + 24.0, 15, ui::TEXT);

    // Легенда фаз справа от заголовка.
    let mut px = r.x + r.w - 16.0;
    for ph in [Phase::Remodeling, Phase::Proliferation, Phase::Inflammation, Phase::Hemostasis] {
        let w = ui.measure(ph.title(), 12, false);
        px -= w;
        ui.text(ph.title(), px, r.y + 23.0, 12, ui::MUTED);
        px -= 13.0;
        draw_rectangle(px, r.y + 15.0, 8.0, 8.0, phase_color(ph));
        px -= 12.0;
    }

    // Легенда серий — кликабельные «чипы», с переносом строк.
    let mut lx = r.x + 10.0;
    let mut ly = r.y + 32.0;
    for (i, name) in NAMES.iter().enumerate() {
        let cw = ui.measure(name, 12, false) + 24.0;
        if lx + cw > r.x + r.w - 10.0 {
            lx = r.x + 10.0;
            ly += 22.0;
        }
        draw_chip(ui, Rect::new(lx, ly, cw, 20.0), i, visible);
        lx += cw + 2.0;
    }

    let plot_top = ly + 34.0;
    let plot = Rect::new(r.x + 46.0, plot_top, r.w - 64.0, r.y + r.h - 28.0 - plot_top);
    let now_d = now_h / 24.0;
    let span = (((now_d / 7.0).ceil() * 7.0).max(28.0)) as usize;
    let xd = |d: f32| plot.x + plot.w * d / span as f32;
    let yv = |v: f32| plot.y + plot.h * (1.0 - v.clamp(0.0, 1.0));

    for k in 0..=4 {
        let v = k as f32 / 4.0;
        let y = yv(v);
        draw_line(plot.x, y, plot.x + plot.w, y, 1.0, ui::GRID);
        ui.text_right(&format!("{v:.2}"), plot.x - 8.0, y + 4.0, 11, ui::MUTED);
    }
    let tick = if span <= 42 {
        7
    } else if span <= 120 {
        14
    } else {
        30
    };
    for d in (0..=span).step_by(tick) {
        let x = xd(d as f32);
        draw_line(x, plot.y, x, plot.y + plot.h, 1.0, ui::GRID);
        let label = format!("{d} д");
        let w = ui.measure(&label, 11, false);
        ui.text(&label, x - w / 2.0, plot.y + plot.h + 16.0, 11, ui::MUTED);
    }

    // Полоса фаз и полоса состояния над графиком.
    for pair in hist.windows(2) {
        let x0 = xd(pair[0].hours / 24.0);
        let x1 = xd(pair[1].hours / 24.0);
        draw_rectangle(x0, plot.y - 16.0, x1 - x0 + 0.6, 4.0, phase_color(pair[0].phase()));
        draw_rectangle(x0, plot.y - 10.0, x1 - x0 + 0.6, 4.0, condition_color(pair[0].condition()));
    }

    // Процедуры — засечки по нижнему краю.
    for ev in events {
        let x = xd(ev.hours / 24.0);
        let (h, th) = match ev.kind {
            EventKind::Debridement => (plot.h, 1.5),
            EventKind::Antiseptic => (10.0, 2.0),
            EventKind::Dose(_) => (6.0, 1.0),
        };
        let mut c = event_color(ev.kind);
        if ev.kind == EventKind::Debridement {
            c.a = 0.5;
        }
        draw_line(x, plot.y + plot.h, x, plot.y + plot.h - h, th, c);
    }

    // Линии серий, с прореживанием до ~ширины графика.
    let stride = (hist.len() / plot.w as usize).max(1);
    for i in (0..N).filter(|&i| visible[i]) {
        let thick = if i == 12 { 3.0 } else { 2.0 };
        let mut prev: Option<Vec2> = None;
        let last = hist.len().saturating_sub(1);
        for (k, m) in hist.iter().enumerate() {
            if k % stride != 0 && k != last {
                continue;
            }
            let pt = vec2(xd(m.hours / 24.0), yv(value(i, m)));
            if let Some(p) = prev {
                draw_line(p.x, p.y, pt.x, pt.y, thick, COLORS[i]);
            }
            prev = Some(pt);
        }
    }

    let nx = xd(now_d);
    draw_line(nx, plot.y - 16.0, nx, plot.y + plot.h, 1.0, Color::new(1.0, 1.0, 1.0, 0.35));

    // Наведение: курсор и значения в выбранный день.
    let m = ui.mouse();
    if plot.contains(m) && !hist.is_empty() {
        let day = (m.x - plot.x) / plot.w * span as f32;
        let idx = ((day * 24.0).round() as usize).min(hist.len() - 1);
        let s = &hist[idx];
        let x = xd(s.hours / 24.0);
        draw_line(x, plot.y, x, plot.y + plot.h, 1.0, Color::new(1.0, 1.0, 1.0, 0.6));
        let rows: Vec<usize> = (0..N).filter(|&i| visible[i]).collect();
        let bw = 214.0;
        let bh = 48.0 + rows.len() as f32 * 17.0;
        let bx = if x + 12.0 + bw > plot.x + plot.w { x - 12.0 - bw } else { x + 12.0 };
        // Прижимаем к низу графика, но не даём уехать выше его верхнего края.
        let by = (plot.y + plot.h - bh).max(plot.y + 4.0);
        ui::fill_rounded(Rect::new(bx, by, bw, bh), 6.0, Color::new(0.05, 0.055, 0.07, 0.97));
        ui.bold(&format!("День {:.1} · {}", s.hours / 24.0, s.phase().title()), bx + 10.0, by + 19.0, 13, ui::TEXT);
        ui.text(s.condition().title(), bx + 10.0, by + 36.0, 12, condition_color(s.condition()));
        for (row, &i) in rows.iter().enumerate() {
            let y = by + 54.0 + row as f32 * 17.0;
            draw_rectangle(bx + 10.0, y - 8.0, 8.0, 8.0, COLORS[i]);
            ui.text(NAMES[i], bx + 24.0, y, 12, ui::MUTED);
            ui.text_right(&value_label(i, s), bx + bw - 10.0, y, 12, ui::TEXT);
        }
    }
}

fn draw_chip(ui: &Ui, chip: Rect, i: usize, visible: &mut [bool; N]) {
    if ui.hovered(chip) {
        ui::fill_rounded(chip, 5.0, ui::BTN_HOVER);
    }
    let on = visible[i];
    let sw = if on { COLORS[i] } else { ui::BORDER };
    draw_rectangle(chip.x + 7.0, chip.y + 6.0, 8.0, 8.0, sw);
    ui.text(NAMES[i], chip.x + 19.0, chip.y + 14.0, 12, if on { ui::TEXT } else { ui::MUTED });
    if ui.clicked_in(chip) {
        visible[i] = !visible[i];
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn series_are_normalized() {
        let m = Metrics {
            open_fraction: 0.5,
            depth: 2.5,
            necrotic_mm2: 10.0,
            wound_mm2: 50.0,
            abx_plasma: 5.0,
            strength: 0.4,
            ..Default::default()
        };
        assert_eq!(value(0, &m), 0.5);
        assert!((value(1, &m) - 0.5).abs() < 1e-6, "глубина 2.5 из 5 мм");
        assert!((value(2, &m) - 0.2).abs() < 1e-6, "некроз — доля площади раны");
        assert!((value(13, &m) - 0.5).abs() < 1e-6);
        assert_eq!(value_label(12, &m), "40%");
        assert_eq!(NAMES.len(), N);
        assert_eq!(COLORS.len(), N);
        assert_eq!(DEFAULT_VISIBLE.len(), N);
    }

    #[test]
    fn every_phase_and_condition_has_a_distinct_color() {
        let phases = [Phase::Hemostasis, Phase::Inflammation, Phase::Proliferation, Phase::Remodeling];
        for (i, a) in phases.iter().enumerate() {
            for b in &phases[i + 1..] {
                assert_ne!(phase_color(*a), phase_color(*b));
            }
        }
        assert_ne!(condition_color(Condition::Healed), condition_color(Condition::Sepsis));
        assert_eq!(
            event_color(EventKind::Dose(body_sim::therapy::Antibiotic::Cefazolin)),
            event_color(EventKind::Dose(body_sim::therapy::Antibiotic::Vancomycin))
        );
    }
}
