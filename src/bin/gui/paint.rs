//! Как выглядит ткань: цвет, влажность и рельеф поверхности для 3D-сцены
//! и «гистологические» текстуры граней блока (срезов через слои кожи).

use body_sim::params::Params;
use body_sim::report::View;
use body_sim::therapy::Debrider;
use body_sim::tissue::Tissue;
use rayon::prelude::*;

pub type Rgb = [f32; 3];

const fn rgb(r: u8, g: u8, b: u8) -> Rgb {
    [r as f32 / 255.0, g as f32 / 255.0, b as f32 / 255.0]
}

pub const SKIN: Rgb = rgb(226, 178, 156);
pub const BLOOD: Rgb = rgb(120, 6, 16);
pub const CLOT: Rgb = rgb(70, 10, 14);
pub const RAW: Rgb = rgb(150, 36, 42);
pub const GRANULATION: Rgb = rgb(205, 52, 58);
pub const EPI_NEW: Rgb = rgb(238, 160, 160);
pub const SCAR_YOUNG: Rgb = rgb(220, 120, 126);
pub const SCAR_OLD: Rgb = rgb(238, 212, 200);
pub const EXUDATE: Rgb = rgb(212, 194, 104);
pub const SLOUGH: Rgb = rgb(204, 184, 112);
pub const ESCHAR: Rgb = rgb(34, 26, 22);
pub const FAT: Rgb = rgb(244, 206, 110);
pub const MUSCLE: Rgb = rgb(152, 38, 44);
pub const ERYTHEMA: Rgb = rgb(214, 78, 70);
pub const BIOFILM: Rgb = rgb(146, 166, 90);
pub const IODINE: Rgb = rgb(160, 76, 20);

const EPIDERMIS: Rgb = rgb(186, 132, 104);
const CORNEUM: Rgb = rgb(222, 196, 176);
const DERMIS: Rgb = rgb(236, 190, 178);
const PAPILLARY: Rgb = rgb(236, 170, 164);
const SEPTUM: Rgb = rgb(232, 196, 178);
const FASCIA: Rgb = rgb(236, 232, 224);
const FOLLICLE: Rgb = rgb(150, 100, 78);
const HAIR: Rgb = rgb(52, 36, 28);
const GLAND: Rgb = rgb(214, 160, 170);
const CAPILLARY: Rgb = rgb(170, 20, 30);
const BACTERIA: Rgb = rgb(118, 172, 48);
const NEUTROPHIL: Rgb = rgb(110, 64, 160);
const AIR: Rgb = rgb(14, 14, 18);
const LARVA: Rgb = rgb(236, 226, 196);
const FIBROSIS: Rgb = rgb(226, 206, 196);

pub const LEGEND: [(&str, Rgb); 12] = [
    ("Кожа", SKIN),
    ("Кровь", BLOOD),
    ("Грануляции", GRANULATION),
    ("Новый эпителий", EPI_NEW),
    ("Молодой рубец", SCAR_YOUNG),
    ("Зрелый рубец", SCAR_OLD),
    ("Клетчатка", FAT),
    ("Мышца", MUSCLE),
    ("Слаф", SLOUGH),
    ("Сухой струп", ESCHAR),
    ("Бактерии", BACTERIA),
    ("Нейтрофилы", NEUTROPHIL),
];

fn mix(a: Rgb, b: Rgb, t: f32) -> Rgb {
    let t = t.clamp(0.0, 1.0);
    [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t]
}

fn scale(a: Rgb, k: f32) -> Rgb {
    [a[0] * k, a[1] * k, a[2] * k]
}

fn sat(x: f32, half: f32) -> f32 {
    x / (x + half)
}

fn smoothstep(e0: f32, e1: f32, x: f32) -> f32 {
    let t = ((x - e0) / (e1 - e0)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// Перцептивно ровная палитра в духе inferno.
pub fn colormap(v: f32) -> Rgb {
    const STOPS: [Rgb; 5] = [
        rgb(4, 5, 18),
        rgb(70, 16, 105),
        rgb(186, 54, 85),
        rgb(249, 142, 9),
        rgb(252, 255, 164),
    ];
    let x = v.clamp(0.0, 1.0) * (STOPS.len() - 1) as f32;
    let k = (x.floor() as usize).min(STOPS.len() - 2);
    mix(STOPS[k], STOPS[k + 1], x - k as f32)
}

// ---------------------------------------------------------------- шум

fn hash_u(x: i32, y: i32, seed: u32) -> u32 {
    let mut h = (x as u32).wrapping_mul(374_761_393) ^ (y as u32).wrapping_mul(668_265_263) ^ seed.wrapping_mul(2_246_822_519);
    h = (h ^ (h >> 13)).wrapping_mul(1_274_126_177);
    h ^ (h >> 16)
}

fn hash(x: i32, y: i32, seed: u32) -> f32 {
    hash_u(x, y, seed) as f32 / u32::MAX as f32
}

fn vnoise(x: f32, y: f32, seed: u32) -> f32 {
    let (ix, iy) = (x.floor() as i32, y.floor() as i32);
    let (fx, fy) = (x - ix as f32, y - iy as f32);
    let (ux, uy) = (fx * fx * (3.0 - 2.0 * fx), fy * fy * (3.0 - 2.0 * fy));
    let a = hash(ix, iy, seed);
    let b = hash(ix + 1, iy, seed);
    let c = hash(ix, iy + 1, seed);
    let d = hash(ix + 1, iy + 1, seed);
    let top = a + (b - a) * ux;
    let bot = c + (d - c) * ux;
    top + (bot - top) * uy
}

fn fbm(x: f32, y: f32, seed: u32) -> f32 {
    0.5 * vnoise(x, y, seed) + 0.3 * vnoise(2.03 * x, 2.03 * y, seed + 1) + 0.2 * vnoise(4.1 * x, 4.1 * y, seed + 2)
}

/// Шум Уорли: расстояния до ближайшей и второй по близости «клеточной» точки (в единицах ячейки).
fn worley(x: f32, y: f32, seed: u32) -> (f32, f32) {
    let (ix, iy) = (x.floor() as i32, y.floor() as i32);
    let (mut d1, mut d2) = (9.0f32, 9.0f32);
    for dy in -1..=1 {
        for dx in -1..=1 {
            let (cx, cy) = (ix + dx, iy + dy);
            let px = cx as f32 + hash(cx, cy, seed);
            let py = cy as f32 + hash(cx, cy, seed + 17);
            let d = ((px - x).powi(2) + (py - y).powi(2)).sqrt();
            if d < d1 {
                d2 = d1;
                d1 = d;
            } else if d < d2 {
                d2 = d;
            }
        }
    }
    (d1, d2)
}

/// Есть ли в точке «клетка» размером `r` (в ячейках сетки шага 1) с вероятностью `density`.
fn dots(x: f32, y: f32, density: f32, r: f32, seed: u32) -> bool {
    if density <= 0.001 {
        return false;
    }
    let (ix, iy) = (x.floor() as i32, y.floor() as i32);
    for dy in -1..=1 {
        for dx in -1..=1 {
            let (cx, cy) = (ix + dx, iy + dy);
            if hash(cx, cy, seed + 5) > density {
                continue;
            }
            let px = cx as f32 + 0.2 + 0.6 * hash(cx, cy, seed);
            let py = cy as f32 + 0.2 + 0.6 * hash(cx, cy, seed + 9);
            if (px - x).powi(2) + (py - y).powi(2) < r * r {
                return true;
            }
        }
    }
    false
}

// ---------------------------------------------------------------- поверхность

/// Внешний вид одной клетки сетки на поверхности.
#[derive(Clone, Copy, Default)]
pub struct CellLook {
    pub albedo: Rgb,
    /// 0 — сухо, 1 — мокро и блестит.
    pub wet: f32,
    /// 1 — неповреждённая кожа с порами.
    pub skin: f32,
    /// Затенение полости (ambient occlusion).
    pub ao: f32,
}

/// Глубина видимой поверхности (мм вниз от уровня кожи; отрицательная — отёк).
/// На месте невосстановленной клетчатки и фиброзной мышцы зажившая рана остаётся вдавленной.
pub fn surface_z(t: &Tissue, p: &Params, i: usize) -> f32 {
    let e = t.epithelium.data[i];
    let open = (t.depth.data[i] - t.slough.data[i]).max(0.0);
    // Свежая рана залита кровью, потом сгустком.
    let fill = ((t.bleeding.data[i] + t.clot.data[i] * (1.0 - e)) * 0.9).min(0.9);
    let inflamed = (0.8 * t.neutrophils.data[i] + t.bacteria_total(i)).min(1.0) * e;
    let dmax = t.depth_max.data[i];
    let healed = (1.0 - t.depth.data[i] / 0.3).clamp(0.0, 1.0) * t.collagen.data[i];
    let missing = p.lost_fat_mm(dmax) * (1.0 - t.fat_new.data[i]) + 0.5 * p.lost_muscle_mm(dmax) * t.muscle_scar.data[i];
    open * (1.0 - fill) - 0.25 * inflamed + p.atrophy * missing * healed
}

fn bed_color(p: &Params, dmax: f32) -> Rgb {
    if dmax <= p.skin_bottom_mm() {
        RAW
    } else if dmax <= p.fat_bottom_mm() {
        mix(RAW, FAT, ((dmax - p.skin_bottom_mm()) / 0.8).min(1.0))
    } else {
        MUSCLE
    }
}

/// Мёртвая ткань: влажный жёлтый слаф при инфекции, сухой чёрный струп при ишемии.
/// Под гидрогелем струп размягчается и светлеет.
fn necrosis_color(slough: f32, bacteria: f32, moist: bool) -> Rgb {
    mix(SLOUGH, ESCHAR, necrosis_dryness(slough, bacteria, moist))
}

fn necrosis_dryness(slough: f32, bacteria: f32, moist: bool) -> f32 {
    (slough / 1.2).min(1.0) * (1.0 - 0.7 * bacteria.min(1.0)) * if moist { 0.3 } else { 1.0 }
}

pub fn cell_look(t: &Tissue, p: &Params, i: usize, view: Option<View>) -> CellLook {
    let depth_vis = surface_z(t, p, i).max(0.0);
    let moist = t.debriders[Debrider::Hydrogel.index()];
    let ao = 1.0 - 0.06 * depth_vis.min(7.0);
    if let Some(v) = view {
        return CellLook { albedo: colormap(v.value(t, i)), wet: 0.25, skin: 0.0, ao };
    }
    let c = t.collagen.data[i];
    let q = t.maturity.data[i];
    let v = t.vessels.data[i];
    let f = t.fibroblasts.data[i];
    let e = t.epithelium.data[i];
    let cf = t.clot.data[i];
    let bl = t.bleeding.data[i];
    let bt = t.bacteria_total(i);
    let n = t.neutrophils.data[i];
    let dmax = t.depth_max.data[i];
    let slough = t.slough.data[i];
    let bf = t.biofilm.data[i];
    let asc = t.antiseptic.data[i];

    // Дно: оголённый слой → сочные грануляции по мере прорастания сосудов.
    let mut col = mix(bed_color(p, dmax), GRANULATION, 0.6 * v + 0.6 * f);
    let ripeness = (q / 0.8).clamp(0.0, 1.0);
    let mature = mix(SCAR_OLD, SKIN, (q - 0.8) / 0.2);
    col = mix(col, mix(SCAR_YOUNG, mature, ripeness), c);
    col = mix(col, mix(col, EPI_NEW, 0.6), e * (1.0 - c * ripeness));
    let inflamed = (0.8 * n + 0.3 * t.signal.data[i].min(1.0) + bt).min(1.0) * e;
    col = mix(col, ERYTHEMA, 0.55 * inflamed);
    let open = 1.0 - e;
    col = mix(col, CLOT, cf * open * 0.85 * (1.0 - sat(slough, 0.3)));
    col = mix(col, BLOOD, bl);
    let necro = sat(slough, 0.25) * 1.1;
    col = mix(col, necrosis_color(slough, bt, moist), necro);
    col = mix(col, BIOFILM, 0.5 * bf * open);
    let exudate = ((2.0 * bt).min(1.0) * 0.6 + 0.1 * n) * open;
    col = mix(col, EXUDATE, exudate.min(0.7));
    if t.antiseptic_agent == body_sim::therapy::Antiseptic::PovidoneIodine {
        col = mix(col, IODINE, 0.55 * asc);
    }

    // Влажность: открытая рана мокрая, кровь и биоплёнка блестят, сухой струп — матовый.
    let dry = necrosis_dryness(slough, bt, moist) * necro.min(1.0);
    let mut wet = 0.12 * e + open * 0.75;
    wet = wet.max(bl).max(bf * open).max(asc * 0.9).max(0.35 * e * (1.0 - c * ripeness));
    wet *= 1.0 - 0.9 * dry;
    // Поры и микрорельеф есть только у нормальной кожи — у рубца их нет.
    let skin = e * ((q - 0.8) / 0.2).clamp(0.0, 1.0) * (1.0 - bl);
    CellLook { albedo: col, wet: wet.clamp(0.0, 1.0), skin, ao }
}

// ---------------------------------------------------------------- срезы

/// Геометрия текстуры грани: по горизонтали — вдоль линии клеток, по вертикали — глубина.
pub struct SectionGeom {
    pub w: usize,
    pub h: usize,
    /// Сколько миллиметров над поверхностью кожи попадает в текстуру.
    pub top_mm: f32,
    pub px_per_mm: f32,
}

impl SectionGeom {
    pub fn v_of_z(&self, z_mm: f32) -> f32 {
        (z_mm + self.top_mm) * self.px_per_mm / self.h as f32
    }
}

fn put(buf: &mut [u8], o: usize, col: Rgb) {
    for k in 0..3 {
        buf[o + k] = (col[k].clamp(0.0, 1.0) * 255.0) as u8;
    }
    buf[o + 3] = 255;
}

/// Ткань исходного слоя на глубине `z` (мм) в точке `u` вдоль среза.
fn native_layer(p: &Params, u: f32, z: f32, v: f32, seed: u32) -> Rgb {
    let skin = p.skin_bottom_mm();
    let fat = p.fat_bottom_mm();
    let grain = fbm(u * 3.0, z * 3.0, seed);
    if z < p.epidermis_mm {
        return if z < 0.025 { CORNEUM } else if z > 0.08 { scale(EPIDERMIS, 0.85) } else { EPIDERMIS };
    }
    if z < skin {
        // Волосяной фолликул: наклонная трубка с волосом и луковицей.
        let bin = (u / 1.6).floor() as i32;
        let u0 = bin as f32 * 1.6 + 0.3 + 0.9 * hash(bin, 0, seed + 31);
        let fz = u0 + 0.35 * z;
        let len = 1.7 + 0.5 * hash(bin, 1, seed + 31);
        if z < len && (u - fz).abs() < 0.07 {
            return if (u - fz).abs() < 0.025 { HAIR } else { FOLLICLE };
        }
        if ((u - (u0 + 0.35 * len)).powi(2) / 0.02 + (z - len).powi(2) / 0.012) < 1.0 {
            return scale(FOLLICLE, 0.9);
        }
        // Потовая железа — клубок у границы с клетчаткой.
        if z > skin - 0.45 && dots(u * 9.0, z * 9.0, 0.35, 0.28, seed + 40) && hash(bin, 2, seed) < 0.6 {
            return GLAND;
        }
        // Капилляры.
        if dots(u * 7.0, z * 7.0, 0.1 + 0.25 * v, 0.14, seed + 50) {
            return CAPILLARY;
        }
        // Волокна коллагена.
        let fiber = 0.5 + 0.5 * (u * 14.0 + 5.0 * fbm(u * 1.5, z * 4.0, seed + 3)).sin();
        let base = if z < p.epidermis_mm + 0.35 { PAPILLARY } else { DERMIS };
        return scale(base, 0.93 + 0.07 * fiber + 0.05 * grain);
    }
    if z < fat {
        // Дольки жира с соединительнотканными перегородками и сосудами в них.
        let (d1, d2) = worley(u / 0.75, z / 0.75, seed + 60);
        if d2 - d1 < 0.07 {
            return if dots(u * 6.0, z * 6.0, 0.25, 0.25, seed + 61) { CAPILLARY } else { SEPTUM };
        }
        return scale(FAT, 0.9 + 0.12 * (1.0 - d1).clamp(0.0, 1.0) + 0.04 * grain);
    }
    if z < fat + 0.3 {
        let line = 0.5 + 0.5 * (z * 90.0 + 3.0 * grain).sin();
        return scale(FASCIA, 0.93 + 0.07 * line);
    }
    // Мышца: пучки волокон (перимизий) и поперечная исчерченность.
    let (d1, d2) = worley(u / 1.6, z / 0.5, seed + 70);
    if d2 - d1 < 0.06 {
        return scale(SEPTUM, 0.85);
    }
    let stria = 0.5 + 0.5 * (u * 45.0).sin();
    let fiber = 0.5 + 0.5 * (z * 55.0 + 2.0 * grain).sin();
    scale(MUSCLE, 0.82 + 0.1 * fiber + 0.06 * stria + 0.06 * grain)
}

/// Состояние ткани в одном вертикальном столбце среза (интерполяция между клетками).
#[derive(Clone, Copy)]
struct Column {
    u: f32,
    depth: f32,
    dmax: f32,
    slough: f32,
    e: f32,
    c: f32,
    q: f32,
    v: f32,
    f: f32,
    cf: f32,
    bl: f32,
    nt: f32,
    bf: f32,
    bt: f32,
    top: f32,
    myo: f32,
    mscar: f32,
    fat_new: f32,
    moist: bool,
    larvae: bool,
}

impl Column {
    /// Нетронутая кожа — такой столбец берётся из готовой анатомической подложки.
    fn intact(u: f32) -> Self {
        Self { u, depth: 0.0, dmax: 0.0, slough: 0.0, e: 1.0, c: 1.0, q: 1.0, v: 1.0, f: 0.1, cf: 0.0, bl: 0.0, nt: 0.0, bf: 0.0, bt: 0.0, top: 0.0, myo: 0.0, mscar: 0.0, fat_new: 0.0, moist: false, larvae: false }
    }

    fn is_quiet(&self) -> bool {
        self.dmax < 0.01
            && self.slough < 0.01
            && self.nt < 0.02
            && self.bt < 0.005
            && self.e > 0.98
            && self.top.abs() < 0.02
            && self.v > 0.95
    }

    fn sample(t: &Tissue, p: &Params, i0: usize, i1: usize, fs: f32, u: f32) -> Self {
        let at = |fld: &body_sim::grid::Field| fld.data[i0] + (fld.data[i1] - fld.data[i0]) * fs;
        Self {
            u,
            depth: at(&t.depth),
            dmax: at(&t.depth_max),
            slough: at(&t.slough),
            e: at(&t.epithelium),
            c: at(&t.collagen),
            q: at(&t.maturity),
            v: at(&t.vessels),
            f: at(&t.fibroblasts),
            cf: at(&t.clot),
            bl: at(&t.bleeding),
            nt: at(&t.neutrophils),
            bf: at(&t.biofilm),
            bt: t.bacteria_total(i0) + (t.bacteria_total(i1) - t.bacteria_total(i0)) * fs,
            top: surface_z(t, p, i0) + (surface_z(t, p, i1) - surface_z(t, p, i0)) * fs,
            myo: at(&t.myo),
            mscar: at(&t.muscle_scar),
            fat_new: at(&t.fat_new),
            moist: t.debriders[Debrider::Hydrogel.index()],
            larvae: t.debriders[Debrider::Larvae.index()],
        }
    }
}

/// Цвет пикселя среза на глубине `z` (мм) в столбце `k`.
fn section_pixel(k: &Column, p: &Params, z: f32, seed: u32) -> Rgb {
    let u = k.u;
    let open_top = (k.depth - k.slough).max(0.0);
    let grain = fbm(u * 5.0, z * 5.0, seed + 7);
    let col = if z < k.top - 0.02 {
        AIR
    } else if z < open_top {
        // Кровь и сгусток в полости: нити фибрина.
        let fib = smoothstep(0.55, 0.75, fbm(u * 6.0, z * 2.0, seed + 11));
        mix(mix(CLOT, BLOOD, k.bl), scale(CLOT, 1.6), fib * k.cf)
    } else if z < k.depth {
        // Мёртвая ткань с волокнистой структурой; сверху — биоплёнка.
        if z - open_top < 0.25 * k.bf {
            scale(BIOFILM, 0.9 + 0.2 * grain)
        } else if k.larvae && dots(u * 7.0, z * 12.0, 0.35, 0.3, seed + 15) {
            LARVA
        } else {
            let fibrous = 0.8 + 0.3 * fbm(u * 10.0, z * 3.0, seed + 13);
            scale(necrosis_color(k.slough, k.bt, k.moist), fibrous)
        }
    } else if z < k.dmax && z >= p.fat_bottom_mm() + 0.3 && k.myo + k.mscar > 0.01 {
        // Там, где была мышца: дольки новых волокон, фиброз и ещё не созревшие грануляции.
        let sel = vnoise(u / 1.6, z / 0.5, seed + 91);
        if sel < k.myo {
            scale(native_layer(p, u, z, 1.0, seed), 1.08)
        } else if sel < k.myo + k.mscar {
            scale(FIBROSIS, 0.92 + 0.08 * (0.5 + 0.5 * (z * 60.0 + 2.0 * grain).sin()))
        } else {
            new_scar_tissue(k, u, z, grain, seed)
        }
    } else if z < k.dmax && z >= p.skin_bottom_mm() && z < p.fat_bottom_mm() && k.fat_new > 0.01 {
        // Там, где была клетчатка: вернулись лишь отдельные жировые дольки, остальное — рубец.
        if vnoise(u / 0.75, z / 0.75, seed + 90) < k.fat_new {
            native_layer(p, u, z, 1.0, seed)
        } else {
            new_scar_tissue(k, u, z, grain, seed)
        }
    } else if z < k.dmax {
        new_scar_tissue(k, u, z, grain, seed)
    } else {
        native_layer(p, u, z, k.v, seed)
    };
    overlay_cells(k, z, col, seed)
}

/// Новая ткань на месте дефекта: грануляции с петлями капилляров → рубец с параллельными волокнами.
fn new_scar_tissue(k: &Column, u: f32, z: f32, grain: f32, seed: u32) -> Rgb {
    let gran = (1.0 - k.c).clamp(0.0, 1.0);
    if gran > 0.3 && dots(u * 8.0, z * 8.0, 0.25 * (k.v + k.f).min(1.0) * gran, 0.18, seed + 20) {
        CAPILLARY
    } else {
        let ripeness = (k.q / 0.8).clamp(0.0, 1.0);
        let scar = mix(SCAR_YOUNG, SCAR_OLD, ripeness);
        let lines = 0.5 + 0.5 * (z * 70.0 + 1.5 * grain).sin();
        let base = mix(GRANULATION, scar, k.c);
        scale(base, 0.9 + 0.08 * lines * k.c + 0.08 * grain)
    }
}

/// Поверх ткани: новый эпителий, воспаление, нейтрофилы и бактерии у поверхности.
fn overlay_cells(k: &Column, z: f32, mut col: Rgb, seed: u32) -> Rgb {
    let u = k.u;
    let open_top = (k.depth - k.slough).max(0.0);
    let live_top = k.depth.max(k.top);
    if z >= live_top - 0.01 {
        // Новый эпителий поверх заполненного дефекта.
        if k.dmax > 0.05 && k.e > 0.05 && z - live_top < 0.1 {
            col = mix(col, EPI_NEW, k.e);
        }
        // Воспаление у поверхности, нейтрофилы и прорастающие бактерии.
        let near = (1.0 - (z - live_top) / 1.8).clamp(0.0, 1.0);
        col = mix(col, ERYTHEMA, 0.3 * (k.nt + k.bt).min(1.0) * near);
        if dots(u * 14.0, z * 14.0, 0.5 * k.nt * near, 0.25, seed + 30) {
            col = NEUTROPHIL;
        }
        if dots(u * 22.0, z * 22.0, 0.8 * k.bt * near, 0.3, seed + 33) {
            col = BACTERIA;
        }
    } else if z >= open_top && dots(u * 22.0, z * 22.0, 0.9 * k.bt, 0.3, seed + 34) {
        col = BACTERIA;
    }
    col
}

/// Анатомическая подложка грани: нетронутые слои кожи. Считается один раз.
pub fn paint_base(p: &Params, n_cells: usize, g: &SectionGeom, out: &mut [u8], seed: u32) {
    let len_mm = n_cells as f32 * p.cell_mm;
    out.par_chunks_mut(g.w * 4).enumerate().for_each(|(py, row)| {
        let z = (py as f32 + 0.5) / g.px_per_mm - g.top_mm;
        for px in 0..g.w {
            let k = Column::intact((px as f32 + 0.5) / g.w as f32 * len_mm);
            put(row, px * 4, section_pixel(&k, p, z, seed));
        }
    });
}

/// Грань блока: срез по линии клеток `cells` (соседние клетки сетки, шаг `cell_mm`).
/// Столбцы с нетронутой тканью копируются из подложки `base`, остальные рисуются заново;
/// строки считаются параллельно.
pub fn paint_section(t: &Tissue, p: &Params, cells: &[usize], g: &SectionGeom, base: &[u8], out: &mut [u8], seed: u32) {
    let n = cells.len();
    let len_mm = n as f32 * p.cell_mm;
    let columns: Vec<Column> = (0..g.w)
        .map(|px| {
            let s = ((px as f32 + 0.5) / g.w as f32 * n as f32 - 0.5).clamp(0.0, (n - 1) as f32);
            let k0 = s.floor() as usize;
            let k1 = (k0 + 1).min(n - 1);
            let u = (px as f32 + 0.5) / g.w as f32 * len_mm;
            Column::sample(t, p, cells[k0], cells[k1], s - k0 as f32, u)
        })
        .collect();
    let quiet: Vec<bool> = columns.iter().map(Column::is_quiet).collect();
    let row_bytes = g.w * 4;
    out.par_chunks_mut(row_bytes).zip(base.par_chunks(row_bytes)).enumerate().for_each(|(py, (row, base_row))| {
        let z = (py as f32 + 0.5) / g.px_per_mm - g.top_mm;
        for px in 0..g.w {
            let o = px * 4;
            if quiet[px] {
                row[o..o + 4].copy_from_slice(&base_row[o..o + 4]);
            } else {
                put(row, o, section_pixel(&columns[px], p, z, seed));
            }
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use body_sim::tissue::WoundShape;

    const W: usize = 48;
    const H: usize = 24;

    fn geom(cells: usize) -> SectionGeom {
        SectionGeom { w: cells * 4, h: 120, top_mm: 1.0, px_per_mm: 9.0 }
    }

    fn row(t: &Tissue, y: usize) -> Vec<usize> {
        (0..t.w).map(|x| y * t.w + x).collect()
    }

    #[test]
    fn colormap_spans_dark_to_light() {
        let lum = |c: Rgb| c[0] + c[1] + c[2];
        assert!(lum(colormap(0.0)) < 0.2);
        assert!(lum(colormap(1.0)) > 2.4);
        assert_eq!(colormap(-5.0), colormap(0.0));
        assert_eq!(colormap(5.0), colormap(1.0));
        assert!(lum(colormap(0.3)) < lum(colormap(0.7)));
    }

    #[test]
    fn noise_is_bounded_and_deterministic() {
        for k in 0..200 {
            let (x, y) = (k as f32 * 0.37, k as f32 * 0.91);
            let v = fbm(x, y, 3);
            assert!((0.0..=1.0).contains(&v));
            assert_eq!(v, fbm(x, y, 3));
            let (d1, d2) = worley(x, y, 5);
            assert!(d1 <= d2);
        }
    }

    #[test]
    fn surface_of_healthy_skin_is_flat() {
        let p = Params::default();
        let t = Tissue::healthy(W, H, &p);
        assert!(surface_z(&t, &p, 0).abs() < 1e-6);
    }

    #[test]
    fn fresh_wound_is_filled_with_blood_and_healed_fat_defect_is_depressed() {
        let p = Params::default();
        let mut t = Tissue::healthy(W, H, &p);
        t.injure(WoundShape::Circle { radius: 6.0 }, 5.0, &p);
        let c = (H / 2) * W + W / 2;
        let fresh = surface_z(&t, &p, c);
        assert!(fresh > 0.0 && fresh < 1.0, "полость залита кровью: {fresh}");
        // «Зажившая» клетка: полость заполнена, коллаген есть, жир не вернулся.
        t.depth.data[c] = 0.0;
        t.bleeding.data[c] = 0.0;
        t.clot.data[c] = 0.0;
        t.collagen.data[c] = 1.0;
        t.epithelium.data[c] = 1.0;
        t.bacteria.data[c] = 0.0; // иначе поверхность чуть приподнимет воспалительный отёк
        let dip = surface_z(&t, &p, c);
        assert!((dip - p.atrophy * p.lost_fat_mm(5.0)).abs() < 1e-4, "вдавленный рубец: {dip}");
        t.fat_new.data[c] = 1.0;
        assert!(surface_z(&t, &p, c).abs() < 1e-4, "жир вернулся — вдавления нет");
    }

    #[test]
    fn looks_of_skin_wound_and_eschar() {
        let p = Params::default();
        let mut t = Tissue::healthy(W, H, &p);
        let skin = cell_look(&t, &p, 0, None);
        assert!(skin.skin > 0.99 && skin.wet < 0.2);
        t.injure(WoundShape::Circle { radius: 6.0 }, 2.5, &p);
        let c = (H / 2) * W + W / 2;
        let wound = cell_look(&t, &p, c, None);
        assert!(wound.wet > 0.9 && wound.skin == 0.0, "свежая рана мокрая, без пор");
        t.bleeding.data[c] = 0.0;
        t.slough.data[c] = 2.0;
        let eschar = cell_look(&t, &p, c, None);
        assert!(eschar.wet < 0.3, "сухой струп матовый");
        t.debriders[Debrider::Hydrogel.index()] = true;
        assert!(cell_look(&t, &p, c, None).wet > eschar.wet, "гидрогель размягчает струп");
        let heat = cell_look(&t, &p, c, Some(View::Necrosis));
        assert_eq!(heat.albedo, colormap(1.0));
    }

    #[test]
    fn section_of_healthy_skin_is_the_anatomy_base() {
        let p = Params::default();
        let t = Tissue::healthy(W, H, &p);
        let g = geom(W);
        let mut base = vec![0u8; g.w * g.h * 4];
        paint_base(&p, W, &g, &mut base, 1);
        let mut out = vec![0u8; base.len()];
        paint_section(&t, &p, &row(&t, H / 2), &g, &base, &mut out, 1);
        assert_eq!(out, base, "нетронутая кожа копируется из подложки");
        // Слои по глубине: воздух над кожей, жир в клетчатке, мышца внизу.
        let px = |z_mm: f32| {
            let py = ((z_mm + g.top_mm) * g.px_per_mm) as usize;
            let o = (py * g.w + g.w / 2) * 4;
            [base[o], base[o + 1], base[o + 2]]
        };
        assert!(px(-0.5).iter().all(|&v| v < 30), "над кожей — воздух");
        assert!(px(5.0)[1] > 150, "клетчатка жёлтая");
        let muscle = px(10.5);
        assert!(muscle[0] > muscle[1] && muscle[0] > muscle[2], "мышца красная");
    }

    #[test]
    fn wound_columns_are_repainted() {
        let p = Params::default();
        let mut t = Tissue::healthy(W, H, &p);
        t.injure(WoundShape::Circle { radius: 6.0 }, 5.0, &p);
        let g = geom(W);
        let mut base = vec![0u8; g.w * g.h * 4];
        paint_base(&p, W, &g, &mut base, 1);
        let mut out = vec![0u8; base.len()];
        paint_section(&t, &p, &row(&t, H / 2), &g, &base, &mut out, 1);
        let col = |buf: &[u8], x: usize| (0..g.h).map(|y| buf[(y * g.w + x) * 4]).collect::<Vec<_>>();
        assert_eq!(col(&out, 0), col(&base, 0), "край среза — здоровая кожа");
        assert_ne!(col(&out, g.w / 2), col(&base, g.w / 2), "в центре — рана");
    }
}