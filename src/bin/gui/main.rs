//! Графический интерфейс симулятора: живая ткань (сверху или в разрезе),
//! пациент, рана, лечение, графики.
//!
//! ЛКМ по ткани — нанести рану выбранной глубины, ПКМ — провести линию разреза,
//! колесо — размер кисти, Пробел — пауза, R — начать заново, → — шаг на 1 час при паузе,
//! Tab — сверху / разрез.

mod chart;
mod paint;
mod prof;
mod render3d;
mod ui;

use body_sim::body::BodySite;
use body_sim::params::Scenario;
use body_sim::report::{Metrics, View};
use body_sim::simulation::Simulation;
use body_sim::therapy::{Antibiotic, Antiseptic, Debrider, EventKind};
use body_sim::tissue::WoundShape;
use macroquad::prelude::*;
use render3d::{OrbitCam, Scene3D};
use ui::Ui;

/// Размеры моделируемого участка: ширина и высота в мм и сторона клетки.
/// Крупный участок берётся более грубой сеткой, иначе шаг модели и меш поверхности
/// перестают укладываться в кадр.
const PATCHES: [(f32, f32, f32, &str); 3] =
    [(24.0, 12.0, 0.25, "24×12 мм"), (48.0, 24.0, 0.25, "48×24"), (96.0, 48.0, 0.5, "96×48")];
/// Размер окна просмотра ткани в пикселях — от размера сетки не зависит.
const VIEW_PX: (f32, f32) = (960.0, 480.0);
const MAX_STEPS_PER_FRAME: usize = 600;
const SPEEDS: [(f32, &str); 6] = [(0.25, "¼"), (0.5, "½"), (1.0, "1"), (2.0, "2"), (5.0, "5"), (10.0, "10")];
const DEPTHS: [(f32, &str, &str); 5] = [
    (0.1, "0.1", "эпидермис: заживает из фолликулов, без рубца"),
    (1.0, "1", "дерма: островки эпителия из придатков, мягкий рубец"),
    (2.5, "2.5", "полнослойная: полость заполняют грануляции, рубец"),
    (5.0, "5", "до клетчатки: дно раны плохо кровоснабжается"),
    (9.0, "9", "до мышцы: долгое заполнение глубокой полости"),
];
/// Глубина колонки ткани (`Params::max_depth_mm`) — предел для отладочного `--depth`.
const MAX_DEPTH_MM: f32 = 12.0;
const DRESSINGS: [(Option<f32>, &str); 3] = [(None, "нет"), (Some(12.0), "12 ч"), (Some(24.0), "24 ч")];
const ABX: [(Option<f32>, &str); 4] = [(None, "нет"), (Some(8.0), "8 ч"), (Some(12.0), "12 ч"), (Some(24.0), "24 ч")];

fn view_rect() -> Rect {
    Rect::new(20.0, 72.0, VIEW_PX.0, VIEW_PX.1)
}

#[derive(Clone, Copy, PartialEq)]
enum Shape {
    Circle,
    Cut,
}

#[derive(Clone, Copy, PartialEq)]
enum Tab {
    Body,
    Wound,
    Treatment,
}

struct App {
    scenario: Scenario,
    site: BodySite,
    /// Индекс в `PATCHES`.
    patch: usize,
    shape: Shape,
    size_mm: f32,
    depth_mm: f32,
    sim: Simulation,
    playing: bool,
    speed: f32,
    pending_steps: f32,
    view: Option<View>,
    /// Срезать ли блок через рану и по какой строке сетки.
    cut: bool,
    cut_row: usize,
    orbiting: bool,
    last_mouse: Vec2,
    brush_mm: f32,
    visible: [bool; chart::N],
    tab: Tab,
    agent: Antiseptic,
    dressing: Option<f32>,
    /// Интервал введения каждого антибиотика (по `Antibiotic::index`).
    abx: [Option<f32>; 3],
}

impl App {
    fn new(scenario: Scenario) -> Self {
        let mut app = Self {
            scenario,
            site: BodySite::Forearm,
            patch: 0,
            shape: Shape::Circle,
            size_mm: 4.0,
            depth_mm: 2.5,
            sim: Self::make_sim(scenario, BodySite::Forearm, 0, Shape::Circle, 4.0, 2.5),
            playing: true,
            speed: 1.0,
            pending_steps: 0.0,
            view: None,
            cut: true,
            cut_row: 24,
            orbiting: false,
            last_mouse: Vec2::ZERO,
            brush_mm: 1.0,
            visible: chart::DEFAULT_VISIBLE,
            tab: Tab::Body,
            agent: Antiseptic::Octenidine,
            dressing: None,
            abx: [None; 3],
        };
        app.reset();
        app
    }

    fn make_sim(
        scenario: Scenario,
        site: BodySite,
        patch: usize,
        shape: Shape,
        size_mm: f32,
        depth_mm: f32,
    ) -> Simulation {
        let (aw, ah, cell, _) = PATCHES[patch];
        let mut p = scenario.params_at(site);
        p.cell_mm = cell;
        p.fit_dt();
        let r = size_mm / p.cell_mm;
        let shape = match shape {
            Shape::Circle => WoundShape::Circle { radius: r },
            Shape::Cut => WoundShape::Cut { half_length: 2.0 * r, half_width: 1.0 / p.cell_mm },
        };
        let (w, h) = ((aw / cell).round() as usize, (ah / cell).round() as usize);
        let depth_mm = depth_mm.min(p.max_depth_mm);
        Simulation::new(p, w, h, shape, depth_mm)
    }

    fn grid(&self) -> (usize, usize) {
        (self.sim.tissue.w, self.sim.tissue.h)
    }

    fn reset(&mut self) {
        self.sim = Self::make_sim(self.scenario, self.site, self.patch, self.shape, self.size_mm, self.depth_mm);
        self.cut_row = self.sim.tissue.h / 2;
        self.pending_steps = 0.0;
        self.dressing = None;
        self.abx = [None; 3];
    }

    fn update(&mut self, frame_dt: f32) {
        if !self.playing {
            return;
        }
        self.pending_steps += frame_dt.min(0.1) * self.speed * 24.0 / self.sim.p.dt;
        let n = self.pending_steps.floor() as usize;
        let n = if n > MAX_STEPS_PER_FRAME {
            self.pending_steps = 0.0;
            MAX_STEPS_PER_FRAME
        } else {
            self.pending_steps -= n as f32;
            n
        };
        self.sim.run_steps(n);
    }

    fn set_dressing(&mut self, every: Option<f32>) {
        self.dressing = every;
        let now = self.sim.hours();
        match every {
            Some(h) => self.sim.therapy.start_antiseptic(self.agent, h, now),
            None => self.sim.therapy.stop_antiseptic(),
        }
    }

    fn set_antibiotic(&mut self, drug: Antibiotic, every: Option<f32>) {
        self.abx[drug.index()] = every;
        let now = self.sim.hours();
        match every {
            Some(h) => self.sim.therapy.start_antibiotic(drug, h, now),
            None => self.sim.therapy.stop_antibiotic(drug),
        }
    }

    fn apply_antiseptic(&mut self) {
        let now = self.sim.hours();
        let s = &mut self.sim;
        s.therapy.apply_antiseptic_now(&mut s.tissue, self.agent, now);
        s.touch();
    }

    fn debride(&mut self) {
        let now = self.sim.hours();
        let s = &mut self.sim;
        s.therapy.debride_now(&mut s.tissue, &s.p, now);
        s.touch();
    }
}

// ---------------------------------------------------------------- шапка

fn badge(ui: &Ui, right: f32, y: f32, label: &str, color: Color) -> f32 {
    let bw = ui.measure(label, 14, true) + 24.0;
    let r = Rect::new(right - bw, y, bw, 26.0);
    ui::fill_rounded(r, 13.0, color);
    ui.bold(label, r.x + 12.0, r.y + 18.0, 14, Color::new(0.06, 0.06, 0.08, 1.0));
    r.x
}

fn draw_header(ui: &Ui, app: &App) {
    ui.bold("Регенерация кожи", 20.0, 34.0, 24, ui::TEXT);
    let wound = match app.shape {
        Shape::Circle => format!("круглая рана Ø{:.0} мм", 2.0 * app.size_mm),
        Shape::Cut => format!("разрез {:.0} мм", 4.0 * app.size_mm),
    };
    let sub = format!(
        "{} · {}, глубина {} мм ({})",
        app.scenario.title(),
        wound,
        app.depth_mm,
        app.sim.p.layer_title(app.depth_mm)
    );
    ui.text(&sub, 20.0, 56.0, 14, ui::MUTED);

    let m = app.sim.latest();
    let right = view_rect().x + view_rect().w;
    let cond = m.condition();
    let x = badge(ui, right, 16.0, cond.title(), chart::condition_color(cond));
    let x = badge(ui, x - 8.0, 16.0, m.phase().title(), chart::phase_color(m.phase()));

    let day = format!("День {:.1}", app.sim.hours() / 24.0);
    let dw = ui.measure(&day, 24, true);
    ui.bold(&day, x - 16.0 - dw, 38.0, 24, ui::TEXT);
    if let Some(h) = app.sim.closed_at {
        ui.text_right(&format!("эпителизация за {:.1} дн.", h / 24.0), right, 60.0, 13, ui::MUTED);
    }
}

// ---------------------------------------------------------------- ткань (3D)

fn shadow_text(ui: &Ui, s: &str, x: f32, y: f32, size: u16, color: Color) {
    ui.text(s, x + 1.0, y + 1.0, size, Color::new(0.0, 0.0, 0.0, 0.7));
    ui.text(s, x, y, size, color);
}

/// Строка сетки, через которую проходит передний срез блока (None — блок целиком).
fn cut_row(app: &App) -> Option<usize> {
    app.cut.then_some(app.cut_row)
}

/// До какой строки (в непрерывных координатах сетки) простирается видимая поверхность.
fn gy_limit(app: &App) -> f32 {
    if app.cut {
        app.cut_row as f32
    } else {
        app.sim.tissue.h as f32 - 0.5
    }
}

/// Управление камерой: ПКМ — вращать, колесо — приблизить, Shift+колесо — кисть, ↑/↓ — срез.
fn camera_input(ui: &Ui, app: &mut App, scene: &mut Scene3D) {
    let r = view_rect();
    let m = ui.mouse();
    if ui.right_clicked() && (r.contains(m) || app.orbiting) {
        if app.orbiting {
            let d = m - app.last_mouse;
            scene.cam.yaw -= d.x * 0.008;
            scene.cam.pitch = (scene.cam.pitch + d.y * 0.006).clamp(0.08, 1.52);
        }
        app.orbiting = true;
    } else {
        app.orbiting = false;
    }
    app.last_mouse = m;

    let wheel = if ui.inert { 0.0 } else { mouse_wheel().1 };
    if wheel != 0.0 && r.contains(m) {
        if is_key_down(KeyCode::LeftShift) || is_key_down(KeyCode::RightShift) {
            app.brush_mm = (app.brush_mm + wheel.signum() * 0.25).clamp(0.5, 4.0);
        } else {
            scene.cam.dist = (scene.cam.dist * if wheel > 0.0 { 0.9 } else { 1.11 }).clamp(14.0, 90.0);
        }
    }
    if is_key_pressed(KeyCode::Up) && app.cut_row > 0 {
        app.cut_row -= 1;
    }
    if is_key_pressed(KeyCode::Down) && app.cut_row + 1 < app.sim.tissue.h {
        app.cut_row += 1;
    }
}

/// 3D-вид ткани. Возвращает клетку под курсором.
fn draw_tissue_3d(ui: &Ui, app: &mut App, scene: &mut Scene3D, prof: &mut prof::Prof) -> Option<usize> {
    let r = view_rect();
    // Фон — мягкий вертикальный градиент, как в студии.
    let bands = 32;
    for k in 0..bands {
        let t = k as f32 / (bands - 1) as f32;
        let y = r.y + r.h * k as f32 / bands as f32;
        let c = Color::new(0.13 - 0.07 * t, 0.135 - 0.07 * t, 0.16 - 0.075 * t, 1.0);
        draw_rectangle(r.x, y, r.w, r.h / bands as f32 + 1.0, c);
    }

    scene.update(&app.sim.tissue, &app.sim.p, app.sim.revision, app.view, cut_row(app));
    prof.add("сцена: вид клеток", scene.timings[0]);
    prof.add("сцена: меши", scene.timings[1]);
    prof.add("сцена: текстуры срезов", scene.timings[2]);
    let pick = if app.orbiting {
        None
    } else {
        prof.time("выбор точки лучом", || scene.pick(r, ui.mouse(), gy_limit(app)))
    };
    if let Some(p) = &pick {
        if ui.mouse_down() {
            let rad = app.brush_mm / app.sim.p.cell_mm;
            app.sim.injure_disk(p.gx, p.gy, rad, app.depth_mm);
        }
    }
    prof.time("3D: отправка на GPU", || scene.draw(r, pick.as_ref(), app.brush_mm));

    // Подписи слоёв у ребра среза — как выноски в атласе.
    if app.view.is_none() {
        let p = &app.sim.p;
        let layers = [
            ("эпидермис и дерма", p.epidermis_mm + p.dermis_mm / 2.0),
            ("жировая клетчатка", p.skin_bottom_mm() + p.fat_mm / 2.0),
            ("фасция", p.fat_bottom_mm() + 0.15),
            ("мышца", p.fat_bottom_mm() + 2.0),
        ];
        let gy = gy_limit(app);
        for (name, depth) in layers {
            let left = scene.project(r, scene.front_edge_point(false, gy, depth));
            let right = scene.project(r, scene.front_edge_point(true, gy, depth));
            let (Some(a), Some(b)) = (left, right) else { continue };
            let w = ui.measure(name, 13, false);
            let (anchor, tx) = if a.x - 40.0 - w > r.x + 6.0 { (a, a.x - 40.0 - w) } else { (b, b.x + 40.0) };
            if tx < r.x || tx + w > r.x + r.w || anchor.y < r.y || anchor.y > r.y + r.h {
                continue;
            }
            let lx = if tx < anchor.x { tx + w + 4.0 } else { tx - 4.0 };
            draw_line(anchor.x, anchor.y, lx, anchor.y, 1.0, Color::new(1.0, 1.0, 1.0, 0.6));
            draw_circle(anchor.x, anchor.y, 2.5, WHITE);
            shadow_text(ui, name, tx, anchor.y + 4.0, 13, WHITE);
        }
    }

    let cut_note = if app.cut {
        format!("срез через y = {:.1} мм (↑/↓)", (app.cut_row as f32 + 0.5) * app.sim.p.cell_mm)
    } else {
        "блок целиком".to_string()
    };
    draw_rectangle(r.x, r.y, r.w, 30.0, Color::new(0.0, 0.0, 0.0, 0.35));
    shadow_text(
        ui,
        &format!(
            "ЛКМ — рана · ПКМ — вращать · колесо — масштаб · Shift+колесо — кисть {:.1} мм · {cut_note}",
            app.brush_mm
        ),
        r.x + 12.0,
        r.y + 20.0,
        12,
        Color::new(1.0, 1.0, 1.0, 0.85),
    );

    pick.map(|p| {
        let (gw, gh) = (app.sim.tissue.w, app.sim.tissue.h);
        let ix = (p.gx.round().max(0.0) as usize).min(gw - 1);
        let iy = (p.gy.round().max(0.0) as usize).min(gh - 1);
        iy * gw + ix
    })
}

fn draw_legend_overlay(ui: &Ui, view: Option<View>) {
    let r = view_rect();
    match view {
        None => {
            let cols = 3;
            let cw = 124.0;
            let rows = paint::LEGEND.len().div_ceil(cols);
            let bw = cols as f32 * cw + 12.0;
            let bh = rows as f32 * 18.0 + 12.0;
            let (bx, by) = (r.x + r.w - bw - 8.0, r.y + r.h - bh - 8.0);
            draw_rectangle(bx, by, bw, bh, Color::new(0.0, 0.0, 0.0, 0.5));
            for (k, (name, c)) in paint::LEGEND.iter().enumerate() {
                let x = bx + 8.0 + (k % cols) as f32 * cw;
                let y = by + 18.0 + (k / cols) as f32 * 18.0;
                draw_rectangle(x, y - 10.0, 11.0, 11.0, Color::new(c[0], c[1], c[2], 1.0));
                ui.text(name, x + 16.0, y, 12, WHITE);
            }
        }
        Some(v) => {
            let w = 200.0;
            let (bx, by) = (r.x + r.w - w - 150.0, r.y + r.h - 40.0);
            draw_rectangle(bx - 30.0, by - 18.0, w + 170.0, 34.0, Color::new(0.0, 0.0, 0.0, 0.55));
            ui.text("0", bx - 20.0, by + 4.0, 12, WHITE);
            for k in 0..w as usize {
                let c = paint::colormap(k as f32 / w);
                draw_rectangle(bx + k as f32, by - 7.0, 1.0, 12.0, Color::new(c[0], c[1], c[2], 1.0));
            }
            ui.text(v.scale_label(), bx + w + 6.0, by + 4.0, 12, WHITE);
            ui.bold(v.title(), bx + w + 50.0, by + 4.0, 12, WHITE);
        }
    }
}

fn draw_view_chips(ui: &Ui, app: &mut App, scene: &mut Scene3D) {
    let r = view_rect();
    let mut x = r.x;
    let mut y = r.y + r.h + 8.0;
    let h = 24.0;
    let chip = |x: f32, y: f32, label: &str, active: bool| -> (bool, f32) {
        let w = ui.measure(label, 13, active) + 20.0;
        (ui.button_sized(Rect::new(x, y, w, h), label, active, 13), w)
    };
    let top = scene.cam.pitch > 1.4;
    let (hit, w) = chip(x, y, "3D", !top);
    if hit {
        scene.cam = OrbitCam::atlas();
    }
    x += w + 4.0;
    let (hit, w) = chip(x, y, "Сверху", top);
    if hit {
        scene.cam = OrbitCam::top();
    }
    x += w + 4.0;
    let (hit, w) = chip(x, y, "Срез", app.cut);
    if hit {
        app.cut = !app.cut;
    }
    x += w + 16.0;
    let views: Vec<Option<View>> = std::iter::once(None).chain(View::ALL.iter().copied().map(Some)).collect();
    for v in views {
        let label = v.map_or("Ткань", |v| v.title());
        let w = ui.measure(label, 13, false) + 20.0;
        if x + w > r.x + r.w {
            x = r.x;
            y += h + 4.0;
        }
        let (hit, w) = chip(x, y, label, app.view == v);
        if hit {
            app.view = v;
        }
        x += w + 4.0;
    }
}

fn draw_tooltip(ui: &Ui, app: &App, i: usize) {
    let t = &app.sim.tissue;
    let p = &app.sim.p;
    let x_mm = (i % t.w) as f32 * p.cell_mm;
    let y_mm = (i / t.w) as f32 * p.cell_mm;
    let rows = [
        ("Глубина полости / макс.", format!("{:.1} / {:.1} мм", t.depth.data[i], t.depth_max.data[i])),
        ("Слой дна", p.layer_title(t.depth_max.data[i]).to_string()),
        ("Мёртвая ткань", format!("{:.2} мм", t.slough.data[i])),
        ("Эпителий", format!("{:.2}", t.epithelium.data[i])),
        ("Коллаген / зрелость", format!("{:.2} / {:.2}", t.collagen.data[i], t.maturity.data[i])),
        ("Сосуды / кислород", format!("{:.2} / {:.2}", t.vessels.data[i], t.oxygen.data[i])),
        ("Бактерии (устойч.)", format!("{:.3} ({:.3})", t.bacteria_total(i), t.bacteria_res.data[i])),
        ("Биоплёнка", format!("{:.2}", t.biofilm.data[i])),
        ("Нейтрофилы", format!("{:.2}", t.neutrophils.data[i])),
        ("Макрофаги M1 / M2", format!("{:.2} / {:.2}", t.m1.data[i], t.m2.data[i])),
        ("Фибробласты", format!("{:.2}", t.fibroblasts.data[i])),
        ("Антисептик / антибиотики", format!("{:.2} / {:.1} МПК", t.antiseptic.data[i], t.antibiotic.data[i])),
        ("Мышца: волокна / фиброз", format!("{:.0}% / {:.0}%", t.myo.data[i] * 100.0, t.muscle_scar.data[i] * 100.0)),
        ("Клетчатка: жир вернулся", format!("{:.0}%", t.fat_new.data[i] * 100.0)),
        ("Прочность", format!("{:.0}%", t.strength(i) * 100.0)),
    ];
    let w = 272.0;
    let h = 36.0 + rows.len() as f32 * 18.0;
    let m = ui.mouse();
    let mut x = m.x + 18.0;
    let mut y = m.y + 18.0;
    if x + w > screen_width() - 8.0 {
        x = m.x - 18.0 - w;
    }
    if y + h > screen_height() - 8.0 {
        y = m.y - 18.0 - h;
    }
    ui::fill_rounded(Rect::new(x - 1.0, y - 1.0, w + 2.0, h + 2.0), 8.0, ui::BORDER);
    ui::fill_rounded(Rect::new(x, y, w, h), 7.0, Color::new(0.05, 0.055, 0.07, 0.97));
    ui.bold(&format!("x {x_mm:.1} мм, y {y_mm:.1} мм"), x + 12.0, y + 22.0, 13, ui::TEXT);
    for (k, (name, val)) in rows.iter().enumerate() {
        let ry = y + 42.0 + k as f32 * 18.0;
        ui.text(name, x + 12.0, ry, 13, ui::MUTED);
        ui.text_right(val, x + w - 12.0, ry, 13, ui::TEXT);
    }
}

// ---------------------------------------------------------------- панель

/// Строка сводки: подпись, значение справа и полоска заполнения под ними.
struct MetricRow<'a> {
    name: &'a str,
    value: &'a str,
    /// Насколько заполнена полоска, 0..1.
    frac: f32,
    color: Color,
}

fn metric_row(ui: &Ui, x: f32, y: f32, w: f32, row: &MetricRow) {
    ui.text(row.name, x, y + 12.0, 13, ui::TEXT);
    ui.text_right(row.value, x + w, y + 12.0, 13, ui::MUTED);
    ui::fill_rounded(Rect::new(x, y + 17.0, w, 3.0), 1.5, ui::BTN);
    if row.frac > 0.002 {
        ui::fill_rounded(Rect::new(x, y + 17.0, w * row.frac.clamp(0.0, 1.0), 3.0), 1.5, row.color);
    }
}

fn button_row<T: Copy + PartialEq>(ui: &Ui, x: f32, y: f32, w: f32, items: &[(T, &str)], current: T) -> Option<T> {
    let gap = 6.0;
    let bw = (w - gap * (items.len() - 1) as f32) / items.len() as f32;
    let mut hit = None;
    for (k, (v, label)) in items.iter().enumerate() {
        if ui.button_sized(Rect::new(x + k as f32 * (bw + gap), y, bw, 30.0), label, *v == current, 14) {
            hit = Some(*v);
        }
    }
    hit
}

fn dots(n: f32) -> String {
    let k = (n.clamp(0.0, 1.0) * 4.0).round() as usize;
    "●".repeat(k) + &"○".repeat(4 - k)
}

fn draw_body_tab(ui: &Ui, app: &mut App, x: f32, mut y: f32, iw: f32) {
    ui.section("МЕСТО НА ТЕЛЕ", x, y);
    y += 20.0;
    let quarter = (iw - 18.0) / 4.0;
    for (k, site) in BodySite::ALL.iter().enumerate() {
        let r = Rect::new(x + (k % 4) as f32 * (quarter + 6.0), y + (k / 4) as f32 * 32.0, quarter, 28.0);
        if ui.button_sized(r, site.short(), app.site == *site, 13) && app.site != *site {
            app.site = *site;
            app.reset();
        }
    }
    y += 68.0;
    ui.text(app.site.note(), x, y + 11.0, 12, ui::MUTED);
    y += 20.0;
    let p = &app.sim.p;
    let anatomy = format!(
        "эпидермис {:.2} · дерма {:.1} · клетчатка {:.0} мм · кровоток ×{:.2}",
        p.epidermis_mm, p.dermis_mm, p.fat_mm, p.site_perfusion
    );
    ui.text(&anatomy, x, y + 11.0, 12, ui::MUTED);
    y += 26.0;

    ui.section("УЧАСТОК КОЖИ", x, y);
    y += 20.0;
    let patches: Vec<(usize, &str)> = PATCHES.iter().enumerate().map(|(k, (_, _, _, l))| (k, *l)).collect();
    if let Some(k) = button_row(ui, x, y, iw, &patches, app.patch) {
        app.patch = k;
        app.reset();
    }
    y += 34.0;
    let (gw, gh) = app.grid();
    let region = app.site.region_cm2();
    ui.text(
        &format!("{gw}×{gh} клеток по {:.2} мм · вся область {region:.0} см²", app.sim.p.cell_mm),
        x,
        y + 11.0,
        12,
        ui::MUTED,
    );
    y += 26.0;

    ui.section("СЦЕНАРИЙ", x, y);
    y += 20.0;
    let half = (iw - 6.0) / 2.0;
    for (k, s) in Scenario::ALL.iter().enumerate() {
        let r = Rect::new(x + (k % 2) as f32 * (half + 6.0), y + (k / 2) as f32 * 34.0, half, 29.0);
        if ui.button_sized(r, s.short(), app.scenario == *s, 14) && app.scenario != *s {
            app.scenario = *s;
            app.reset();
        }
    }
}

fn draw_wound_tab(ui: &Ui, app: &mut App, x: f32, mut y: f32, iw: f32) {
    ui.section("РАНА", x, y);
    y += 20.0;
    if let Some(s) = button_row(ui, x, y, iw, &[(Shape::Circle, "Круглая"), (Shape::Cut, "Разрез")], app.shape)
    {
        app.shape = s;
        app.reset();
    }
    y += 36.0;
    if ui.button(Rect::new(x, y, 30.0, 30.0), "−", false) && app.size_mm > 1.0 {
        app.size_mm -= 1.0;
        app.reset();
    }
    let size_label = format!("{:.0} мм", app.size_mm);
    let sw = ui.measure(&size_label, 14, false);
    ui.text(&size_label, x + 30.0 + (60.0 - sw) / 2.0, y + 20.0, 14, ui::TEXT);
    if ui.button(Rect::new(x + 90.0, y, 30.0, 30.0), "+", false) && app.size_mm < 5.0 {
        app.size_mm += 1.0;
        app.reset();
    }
    if ui.button_sized(Rect::new(x + 132.0, y, iw - 132.0, 30.0), "Заново  (R)", false, 14) {
        app.reset();
    }
    y += 42.0;

    ui.section("ГЛУБИНА, ММ", x, y);
    y += 20.0;
    let items: Vec<(f32, &str)> = DEPTHS.iter().map(|(d, l, _)| (*d, *l)).collect();
    if let Some(d) = button_row(ui, x, y, iw, &items, app.depth_mm) {
        app.depth_mm = d;
        app.reset();
    }
    y += 36.0;
    if let Some((_, _, note)) = DEPTHS.iter().find(|(d, _, _)| *d == app.depth_mm) {
        ui.text(note, x, y + 12.0, 13, ui::MUTED);
    }
    y += 30.0;
    ui.text("ЛКМ по ткани — рана выбранной глубины", x, y + 12.0, 12, ui::MUTED);
    ui.text(&format!("колесо — кисть {:.1} мм · ПКМ — линия разреза", app.brush_mm), x, y + 28.0, 12, ui::MUTED);
}

fn draw_treatment_tab(ui: &Ui, app: &mut App, x: f32, mut y: f32, iw: f32) {
    ui.section("АНТИСЕПТИК (МЕСТНО)", x, y);
    y += 20.0;
    let third = (iw - 12.0) / 3.0;
    for (k, a) in Antiseptic::ALL.iter().enumerate() {
        let r = Rect::new(x + (k % 3) as f32 * (third + 6.0), y + (k / 3) as f32 * 32.0, third, 27.0);
        if ui.button_sized(r, a.title(), app.agent == *a, 13) && app.agent != *a {
            app.agent = *a;
            if app.dressing.is_some() {
                app.set_dressing(app.dressing);
            }
        }
    }
    y += 64.0;
    let pr = app.agent.props();
    let info = format!(
        "бактерицидность {} · в биоплёнку {} · в некроз {}",
        dots(pr.kill / 3.0),
        dots(pr.biofilm_pen),
        dots(pr.necro_pen)
    );
    ui.text(&info, x, y + 11.0, 12, ui::MUTED);
    let tox = format!("токсичность: эпителий {} · фибробласты {}", dots(pr.tox_epi / 0.06), dots(pr.tox_fib / 0.06));
    ui.text(&tox, x, y + 27.0, 12, ui::MUTED);
    y += 34.0;
    if ui.button_sized(Rect::new(x, y, 130.0, 28.0), "Обработать", false, 14) {
        app.apply_antiseptic();
    }
    ui.text("перевязки:", x + 142.0, y + 19.0, 13, ui::MUTED);
    if let Some(d) = button_row(ui, x + 216.0, y, iw - 216.0, &DRESSINGS, app.dressing) {
        app.set_dressing(d);
    }
    y += 38.0;

    ui.section("АНТИБИОТИКИ (СИСТЕМНО), ДОЗА КАЖДЫЕ", x, y);
    y += 20.0;
    for drug in Antibiotic::ALL {
        ui.text(drug.title(), x, y + 18.0, 13, ui::TEXT);
        if let Some(a) = button_row(ui, x + 104.0, y, iw - 104.0, &ABX, app.abx[drug.index()]) {
            app.set_antibiotic(drug, a);
        }
        y += 31.0;
    }
    let plasma = app.sim.tissue.abx_plasma;
    let levels = format!(
        "в плазме, МПК: цефазолин {:.1} · клиндамицин {:.1} · ванкомицин {:.1}",
        plasma[0], plasma[1], plasma[2]
    );
    ui.text(&levels, x, y + 11.0, 12, ui::MUTED);
    y += 22.0;

    ui.section("ОЧИЩЕНИЕ ОТ НЕКРОЗА", x, y);
    y += 20.0;
    let quarter = (iw - 18.0) / 4.0;
    if ui.button_sized(Rect::new(x, y, quarter, 28.0), "Хирургия", false, 13) {
        app.debride();
    }
    for (k, d) in Debrider::ALL.iter().enumerate() {
        let r = Rect::new(x + (k + 1) as f32 * (quarter + 6.0), y, quarter, 28.0);
        let on = app.sim.tissue.debriders[d.index()];
        if ui.button_sized(r, d.title(), on, 13) {
            app.sim.tissue.debriders[d.index()] = !on;
            app.sim.touch();
        }
    }
    y += 34.0;
    ui.text("хирургия — сразу · личинки — 1–2 нед. · ферменты и гель — дольше", x, y + 11.0, 12, ui::MUTED);
    y += 18.0;
    let last = |kind: EventKind| app.sim.therapy.events.iter().rev().find(|e| e.kind == kind).map(|e| e.hours / 24.0);
    let fmt = |d: Option<f32>| d.map_or("не было".to_string(), |d| format!("день {d:.1}"));
    let text = format!(
        "последняя хирургия: {} · антисептик: {}",
        fmt(last(EventKind::Debridement)),
        fmt(last(EventKind::Antiseptic))
    );
    ui.text(&text, x, y + 11.0, 12, ui::MUTED);
}

fn draw_metrics(ui: &Ui, m: &Metrics, x: f32, mut y: f32, iw: f32) {
    ui.section("ПОКАЗАТЕЛИ В ОБЛАСТИ РАНЫ", x, y);
    y += 22.0;
    let c = |i: usize| chart::COLORS[i];
    let muscle = if m.lost_muscle_mm > 0.0 {
        format!(
            "−{:.1} мм: волокна {:.0}% · фиброз {:.0}%",
            m.lost_muscle_mm,
            m.muscle_regen * 100.0,
            m.muscle_fibrosis * 100.0
        )
    } else {
        "не задета".to_string()
    };
    let fat = if m.lost_fat_mm > 0.0 {
        format!("−{:.1} мм: жир вернулся {:.0}%", m.lost_fat_mm, m.fat_regen * 100.0)
    } else {
        "не задета".to_string()
    };
    let rows: [(&str, String, f32, Color); 13] = [
        ("Открытая площадь", format!("{:.1} из {:.1} мм²", m.open_mm2, m.wound_mm2), m.open_fraction, c(0)),
        ("Глубина полости", format!("{:.1} мм (макс. {:.1})", m.depth, m.depth_max), m.depth / 5.0, c(1)),
        ("Некроз", format!("{:.1} мм²", m.necrotic_mm2), m.necrotic_mm2 / m.wound_mm2.max(1.0), c(2)),
        ("Бактерии", format!("{:.2} · устойч. {:.0}%", m.bacteria, m.resistant_fraction * 100.0), m.bacteria, c(3)),
        ("Биоплёнка", format!("{:.2}", m.biofilm), m.biofilm, c(5)),
        ("Нейтрофилы · M1 · M2", format!("{:.2} · {:.2} · {:.2}", m.neutrophils, m.m1, m.m2), m.neutrophils, c(6)),
        ("Фибробласты", format!("{:.2}", m.fibroblasts), m.fibroblasts, c(9)),
        ("Сосуды · кислород", format!("{:.2} · {:.2}", m.vessels, m.oxygen), m.vessels, c(10)),
        ("Коллаген", format!("{:.2}", m.collagen), m.collagen, c(11)),
        ("Мышца", muscle, m.muscle_regen, Color::new(0.80, 0.30, 0.34, 1.0)),
        ("Клетчатка", fat, m.fat_regen, Color::new(0.96, 0.80, 0.40, 1.0)),
        ("Прочность рубца", format!("{:.0}%", m.strength * 100.0), m.strength, c(12)),
        (
            "Доля поверхности тела",
            format!("{:.2}% · инф. {:.2}%", m.tbsa_percent, m.infected_tbsa_percent),
            m.tbsa_percent / 10.0,
            c(4),
        ),
    ];
    for (name, val, frac, color) in rows.iter() {
        metric_row(ui, x, y, iw, &MetricRow { name, value: val, frac: *frac, color: *color });
        y += 24.0;
    }
}

fn draw_panel(ui: &Ui, app: &mut App) {
    let pr = Rect::new(1000.0, 72.0, 420.0, 808.0);
    ui::panel(pr);
    let x = pr.x + 16.0;
    let iw = pr.w - 32.0;

    let tw = iw / 3.0;
    for (k, (tab, label)) in [(Tab::Body, "Тело"), (Tab::Wound, "Рана"), (Tab::Treatment, "Лечение")].iter().enumerate()
    {
        if ui.tab(Rect::new(x + k as f32 * tw, pr.y + 6.0, tw, 34.0), label, app.tab == *tab) {
            app.tab = *tab;
        }
    }

    // Время — общее для всех вкладок.
    let mut y = pr.y + 52.0;
    let play_label = if app.playing { "Пауза" } else { "Пуск" };
    if ui.button_sized(Rect::new(x, y, 84.0, 30.0), play_label, !app.playing, 14) {
        app.playing = !app.playing;
    }
    if let Some(v) = button_row(ui, x + 92.0, y, iw - 92.0 - 64.0, &SPEEDS, app.speed) {
        app.speed = v;
    }
    ui.text("дн/сек", x + iw - 54.0, y + 20.0, 12, ui::MUTED);
    y += 44.0;

    match app.tab {
        Tab::Body => draw_body_tab(ui, app, x, y, iw),
        Tab::Wound => draw_wound_tab(ui, app, x, y, iw),
        Tab::Treatment => draw_treatment_tab(ui, app, x, y, iw),
    }

    let m = app.sim.latest().clone();
    draw_metrics(ui, &m, x, pr.y + pr.h - 22.0 - 13.0 * 24.0 - 14.0, iw);
}

// ---------------------------------------------------------------- запуск

struct Opts {
    scenario: Scenario,
    site: Option<BodySite>,
    depth: f32,
    view: Option<View>,
    cut: bool,
    top: bool,
    cam: Option<(f32, f32, f32)>,
    treat_tab: bool,
    days: Option<f32>,
    antiseptic: Option<Antiseptic>,
    antibiotics: Vec<(Antibiotic, f32)>,
    debriders: Vec<Debrider>,
    debride: Vec<f32>,
    screenshot: Option<String>,
    bench: Option<u32>,
}

/// Для отладки и скриншотов: `--scenario`, `--depth`, `--view`, `--no-cut`, `--top`,
/// `--cam yaw,pitch,dist`, `--treat`, `--days` (промотать), `--antiseptic`,
/// `--antibiotic-every`, `--debride 1,3`, `--screenshot out.png`.
fn parse_opts() -> Opts {
    let mut o = Opts {
        scenario: Scenario::Healthy,
        site: None,
        depth: 2.5,
        view: None,
        cut: true,
        top: false,
        cam: None,
        treat_tab: false,
        days: None,
        antiseptic: None,
        antibiotics: Vec::new(),
        debriders: Vec::new(),
        debride: Vec::new(),
        screenshot: None,
        bench: None,
    };
    let mut args = std::env::args().skip(1);
    while let Some(a) = args.next() {
        match a.as_str() {
            "--no-cut" => o.cut = false,
            "--top" => o.top = true,
            "--treat" => o.treat_tab = true,
            _ => {
                let v = args.next().unwrap_or_default();
                match a.as_str() {
                    "--scenario" => o.scenario = Scenario::parse(&v).unwrap_or(Scenario::Healthy),
                    "--site" => o.site = BodySite::parse(&v),
                    // Глубже колонки ткани нельзя: модель всё равно обрежет полость по max_depth_mm.
                    "--depth" => o.depth = v.parse().map_or(2.5, |d: f32| d.clamp(0.0, MAX_DEPTH_MM)),
                    "--view" => o.view = View::parse(&v),
                    "--cam" => {
                        let n: Vec<f32> = v.split(',').filter_map(|s| s.trim().parse().ok()).collect();
                        if n.len() == 3 {
                            o.cam = Some((n[0], n[1], n[2]));
                        }
                    }
                    "--days" => o.days = v.parse().ok(),
                    "--antiseptic" => o.antiseptic = Antiseptic::parse(&v),
                    "--antibiotic-every" => o.antibiotics.extend(v.parse().ok().map(|h| (Antibiotic::Cefazolin, h))),
                    "--antibiotics" => {
                        for item in v.split(',') {
                            if let Some((name, every)) = item.split_once(':') {
                                if let (Some(d), Ok(h)) = (Antibiotic::parse(name.trim()), every.trim().parse()) {
                                    o.antibiotics.push((d, h));
                                }
                            }
                        }
                    }
                    "--debriders" => o.debriders = v.split(',').filter_map(|s| Debrider::parse(s.trim())).collect(),
                    "--debride" => o.debride = v.split(',').filter_map(|s| s.trim().parse().ok()).collect(),
                    "--screenshot" => o.screenshot = Some(v),
                    "--bench" => o.bench = v.parse().ok(),
                    _ => eprintln!("неизвестная опция: {a}"),
                }
            }
        }
    }
    o
}

fn window_conf() -> macroquad::conf::Conf {
    macroquad::conf::Conf {
        miniquad_conf: Conf {
            window_title: "Симулятор регенерации кожи".to_owned(),
            window_width: 1440,
            window_height: 900,
            window_resizable: false,
            sample_count: 4,
            ..Default::default()
        },
        // Меш поверхности — ~42 тыс. вершин за один вызов (индексы u16, поэтому < 65536).
        draw_call_vertex_capacity: 60_000,
        draw_call_index_capacity: 300_000,
        ..Default::default()
    }
}

#[macroquad::main(window_conf)]
async fn main() {
    let opts = parse_opts();
    let mut ui = Ui::load();
    ui.inert = opts.screenshot.is_some() || opts.bench.is_some();
    let mut prof = prof::Prof::new();
    let mut app = App::new(opts.scenario);
    app.depth_mm = opts.depth;
    app.reset();
    app.view = opts.view;
    app.cut = opts.cut;
    if opts.treat_tab {
        app.tab = Tab::Treatment;
    }
    if let Some(site) = opts.site {
        app.site = site;
        app.reset();
    }
    if let Some(a) = opts.antiseptic {
        app.agent = a;
        app.set_dressing(Some(24.0));
    }
    for (drug, h) in &opts.antibiotics {
        app.set_antibiotic(*drug, Some(*h));
    }
    for d in &opts.debriders {
        app.sim.tissue.debriders[d.index()] = true;
    }
    if let Some(days) = opts.days {
        let per_hour = app.sim.steps_per_hour();
        for hour in 0..(days * 24.0).round() as usize {
            if opts.debride.iter().any(|d| (d * 24.0).round() as usize == hour) {
                app.debride();
            }
            app.sim.run_steps(per_hour);
        }
        app.playing = false;
    }

    let (gw, gh) = app.grid();
    let mut scene = Scene3D::new(&app.sim.p, gw, gh);
    if opts.top {
        scene.cam = OrbitCam::top();
    }
    if let Some((yaw, pitch, dist)) = opts.cam {
        scene.cam = OrbitCam { yaw, pitch, dist };
    }

    if opts.bench.is_some() {
        // Худший случай: максимальная скорость и курсор над раной (работает выбор точки).
        app.speed = 10.0;
        app.playing = true;
        let r = view_rect();
        ui.fake_mouse = Some(vec2(r.x + r.w * 0.45, r.y + r.h * 0.35));
    }

    let mut frame = 0u32;
    let mut frame_start = std::time::Instant::now();
    loop {
        ui.begin();
        if is_key_pressed(KeyCode::Space) {
            app.playing = !app.playing;
        }
        if is_key_pressed(KeyCode::R) {
            app.reset();
        }
        if is_key_pressed(KeyCode::Tab) {
            app.cut = !app.cut;
        }
        if is_key_pressed(KeyCode::Right) && !app.playing {
            let n = app.sim.steps_per_hour();
            app.sim.run_steps(n);
        }
        // Смена места на теле или размера участка меняет геометрию блока — пересобираем сцену,
        // сохраняя положение камеры.
        let (gw, gh) = app.grid();
        if !scene.matches(&app.sim.p, gw, gh) {
            let cam = scene.cam;
            scene = Scene3D::new(&app.sim.p, gw, gh);
            scene.cam = cam;
            app.cut_row = app.cut_row.min(gh - 1);
        }
        camera_input(&ui, &mut app, &mut scene);
        let dt = if opts.bench.is_some() { 1.0 / 60.0 } else { get_frame_time() };
        prof.time("модель", || app.update(dt));

        clear_background(ui::BG);
        draw_header(&ui, &app);
        let r = view_rect();
        ui::fill_rounded(Rect::new(r.x - 1.0, r.y - 1.0, r.w + 2.0, r.h + 2.0), 3.0, ui::BORDER);
        let hovered = draw_tissue_3d(&ui, &mut app, &mut scene, &mut prof);
        let t_ui = std::time::Instant::now();
        draw_legend_overlay(&ui, app.view);
        draw_view_chips(&ui, &mut app, &mut scene);
        let now = app.sim.hours();
        chart::draw(
            &ui,
            Rect::new(20.0, 628.0, 960.0, 252.0),
            &app.sim.history,
            &app.sim.therapy.events,
            &mut app.visible,
            now,
        );
        draw_panel(&ui, &mut app);
        if let Some(i) = hovered {
            draw_tooltip(&ui, &app, i);
        }
        prof.add("интерфейс и график", t_ui.elapsed().as_secs_f64() * 1000.0);

        if let Some(path) = &opts.screenshot {
            if frame == 3 {
                get_screen_data().export_png(path);
                break;
            }
        }
        if let Some(n) = opts.bench {
            if frame == n {
                println!("{}", prof.report());
                println!("модельных дней: {:.1}", app.sim.hours() / 24.0);
                break;
            }
        }
        frame += 1;
        next_frame().await;
        prof.add("кадр целиком (с vsync)", frame_start.elapsed().as_secs_f64() * 1000.0);
        prof.end_frame();
        if opts.bench.is_some() && frame == 30 {
            prof.reset();
        }
        frame_start = std::time::Instant::now();
    }
}
