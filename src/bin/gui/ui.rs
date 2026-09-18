//! Минимальный immediate-mode UI поверх примитивов macroquad.

use macroquad::prelude::*;

pub const BG: Color = Color::new(0.067, 0.071, 0.086, 1.0);
pub const PANEL: Color = Color::new(0.102, 0.110, 0.130, 1.0);
pub const BORDER: Color = Color::new(0.190, 0.204, 0.240, 1.0);
pub const GRID: Color = Color::new(0.160, 0.172, 0.200, 1.0);
pub const TEXT: Color = Color::new(0.900, 0.910, 0.930, 1.0);
pub const MUTED: Color = Color::new(0.560, 0.590, 0.650, 1.0);
pub const ACCENT: Color = Color::new(0.290, 0.560, 0.980, 1.0);
pub const BTN: Color = Color::new(0.155, 0.166, 0.195, 1.0);
pub const BTN_HOVER: Color = Color::new(0.205, 0.220, 0.258, 1.0);

pub struct Ui {
    regular: Option<Font>,
    bold: Option<Font>,
    mouse: Vec2,
    clicked: bool,
    down: bool,
    right: bool,
    /// Игнорировать мышь (режим скриншота).
    pub inert: bool,
}

fn load_font(paths: &[&str]) -> Option<Font> {
    paths
        .iter()
        .find_map(|p| std::fs::read(p).ok())
        .and_then(|bytes| load_ttf_font_from_bytes(&bytes).ok())
}

impl Ui {
    /// Встроенный шрифт macroquad не умеет кириллицу, поэтому берём системный.
    pub fn load() -> Self {
        let regular = load_font(&[
            "C:/Windows/Fonts/segoeui.ttf",
            "C:/Windows/Fonts/arial.ttf",
            "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
            "/System/Library/Fonts/Supplemental/Arial.ttf",
        ]);
        let bold = load_font(&[
            "C:/Windows/Fonts/seguisb.ttf",
            "C:/Windows/Fonts/segoeuib.ttf",
            "C:/Windows/Fonts/arialbd.ttf",
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
            "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
        ])
        .or_else(|| regular.clone());
        Self { regular, bold, mouse: Vec2::ZERO, clicked: false, down: false, right: false, inert: false }
    }

    pub fn begin(&mut self) {
        if self.inert {
            self.mouse = vec2(-1e4, -1e4);
            return;
        }
        self.mouse = mouse_position().into();
        self.clicked = is_mouse_button_pressed(MouseButton::Left);
        self.down = is_mouse_button_down(MouseButton::Left);
        self.right = is_mouse_button_down(MouseButton::Right);
    }

    pub fn mouse_down(&self) -> bool {
        self.down
    }

    pub fn mouse(&self) -> Vec2 {
        self.mouse
    }

    pub fn hovered(&self, r: Rect) -> bool {
        r.contains(self.mouse)
    }

    pub fn clicked_in(&self, r: Rect) -> bool {
        self.clicked && self.hovered(r)
    }

    fn font(&self, bold: bool) -> Option<&Font> {
        if bold { self.bold.as_ref() } else { self.regular.as_ref() }
    }

    /// Текст с базовой линией на `y`. Возвращает ширину.
    pub fn text(&self, s: &str, x: f32, y: f32, size: u16, color: Color) -> f32 {
        self.draw(s, x, y, size, color, false)
    }

    pub fn bold(&self, s: &str, x: f32, y: f32, size: u16, color: Color) -> f32 {
        self.draw(s, x, y, size, color, true)
    }

    pub fn text_right(&self, s: &str, right: f32, y: f32, size: u16, color: Color) {
        let w = self.measure(s, size, false);
        self.draw(s, right - w, y, size, color, false);
    }

    pub fn measure(&self, s: &str, size: u16, bold: bool) -> f32 {
        measure_text(s, self.font(bold), size, 1.0).width
    }

    fn draw(&self, s: &str, x: f32, y: f32, size: u16, color: Color, bold: bool) -> f32 {
        let params = TextParams { font: self.font(bold), font_size: size, color, ..Default::default() };
        draw_text_ex(s, x.round(), y.round(), params).width
    }

    pub fn section(&self, title: &str, x: f32, y: f32) {
        self.bold(title, x, y + 12.0, 12, MUTED);
    }

    pub fn button(&self, r: Rect, label: &str, active: bool) -> bool {
        self.button_sized(r, label, active, 15)
    }

    pub fn button_sized(&self, r: Rect, label: &str, active: bool, size: u16) -> bool {
        let hover = self.hovered(r);
        let bg = if active {
            ACCENT
        } else if hover {
            BTN_HOVER
        } else {
            BTN
        };
        fill_rounded(r, 6.0, bg);
        let w = self.measure(label, size, active);
        let color = if active { WHITE } else { TEXT };
        let baseline = r.y + r.h / 2.0 + size as f32 * 0.33;
        self.draw(label, r.x + (r.w - w) / 2.0, baseline, size, color, active);
        hover && self.clicked
    }

    /// Вкладка: подчёркнутый текст.
    pub fn tab(&self, r: Rect, label: &str, active: bool) -> bool {
        let hover = self.hovered(r);
        let color = if active { TEXT } else if hover { TEXT } else { MUTED };
        let w = self.measure(label, 15, true);
        self.bold(label, r.x + (r.w - w) / 2.0, r.y + r.h / 2.0 + 5.0, 15, color);
        let line = if active { ACCENT } else { BORDER };
        draw_rectangle(r.x, r.y + r.h - 2.0, r.w, 2.0, line);
        hover && self.clicked
    }

    pub fn right_clicked(&self) -> bool {
        self.right
    }
}

pub fn fill_rounded(r: Rect, radius: f32, c: Color) {
    let rad = radius.min(r.w / 2.0).min(r.h / 2.0);
    draw_rectangle(r.x + rad, r.y, r.w - 2.0 * rad, r.h, c);
    draw_rectangle(r.x, r.y + rad, rad, r.h - 2.0 * rad, c);
    draw_rectangle(r.x + r.w - rad, r.y + rad, rad, r.h - 2.0 * rad, c);
    for (cx, cy) in [
        (r.x + rad, r.y + rad),
        (r.x + r.w - rad, r.y + rad),
        (r.x + rad, r.y + r.h - rad),
        (r.x + r.w - rad, r.y + r.h - rad),
    ] {
        draw_circle(cx, cy, rad, c);
    }
}

pub fn panel(r: Rect) {
    fill_rounded(Rect::new(r.x - 1.0, r.y - 1.0, r.w + 2.0, r.h + 2.0), 11.0, BORDER);
    fill_rounded(r, 10.0, PANEL);
}
