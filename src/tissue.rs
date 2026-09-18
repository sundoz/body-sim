use crate::grid::Field;
use crate::params::Params;
use crate::therapy::Antiseptic;

/// Участок кожи: набор полей плотностей клеток, матрикса и сигналов.
/// Каждая клетка сетки — колонка ткани «2.5D»: эпидермис → дерма → клетчатка → мышца,
/// из которой сверху может быть вынут дефект глубиной `depth`.
pub struct Tissue {
    pub w: usize,
    pub h: usize,

    // Гемостаз и повреждение
    pub bleeding: Field,
    pub clot: Field,
    pub debris: Field,
    /// Чувствительный к антибиотику штамм.
    pub bacteria: Field,
    /// Устойчивый штамм.
    pub bacteria_res: Field,
    pub biofilm: Field,

    // Глубина (мм)
    /// Текущая глубина полости (ещё не заполненной грануляциями).
    pub depth: Field,
    /// Самый глубокий уровень, до которого дошло поражение (определяет слой и рубец).
    pub depth_max: Field,
    /// Толщина мёртвой ткани (струп/слаф) на дне полости.
    pub slough: Field,

    // Клетки
    pub neutrophils: Field,
    pub m1: Field,
    pub m2: Field,
    pub fibroblasts: Field,

    // Структура ткани
    pub collagen: Field,
    /// Доля зрелого коллагена (тип I, упорядоченные волокна) в `collagen`.
    pub maturity: Field,
    pub vessels: Field,
    pub epithelium: Field,

    // Растворимые сигналы и лекарства
    pub oxygen: Field,
    pub signal: Field,
    pub growth_factor: Field,
    pub vegf: Field,
    /// Местная концентрация антисептика (0..1).
    pub antiseptic: Field,
    /// Концентрация антибиотика в ткани (в единицах МПК чувствительного штамма).
    pub antibiotic: Field,

    pub antiseptic_agent: Antiseptic,
    /// Концентрация антибиотика в плазме.
    pub abx_plasma: f32,

    /// Клетки сетки, относящиеся к ране (растёт при некрозе).
    pub wound_mask: Vec<bool>,
}

#[derive(Clone, Copy, Debug)]
pub enum WoundShape {
    /// Круглый дефект (как после биопсии-панча).
    Circle { radius: f32 },
    /// Резаная рана — вытянутый эллипс.
    Cut { half_length: f32, half_width: f32 },
}

impl Tissue {
    pub fn healthy(w: usize, h: usize, p: &Params) -> Self {
        let f = |v: f32| Field::new(w, h, v);
        Self {
            w,
            h,
            bleeding: f(0.0),
            clot: f(0.0),
            debris: f(0.0),
            bacteria: f(0.0),
            bacteria_res: f(0.0),
            biofilm: f(0.0),
            depth: f(0.0),
            depth_max: f(0.0),
            slough: f(0.0),
            neutrophils: f(0.0),
            m1: f(0.0),
            m2: f(0.0),
            fibroblasts: f(p.fib_baseline),
            collagen: f(1.0),
            maturity: f(1.0),
            vessels: f(1.0),
            epithelium: f(1.0),
            oxygen: f(0.9),
            signal: f(0.0),
            growth_factor: f(0.0),
            vegf: f(0.0),
            antiseptic: f(0.0),
            antibiotic: f(0.0),
            antiseptic_agent: Antiseptic::Octenidine,
            abx_plasma: 0.0,
            wound_mask: vec![false; w * h],
        }
    }

    pub fn len(&self) -> usize {
        self.w * self.h
    }

    /// Наносит повреждение глубиной `depth_mm` в центре участка.
    pub fn injure(&mut self, shape: WoundShape, depth_mm: f32, p: &Params) {
        let cx = (self.w as f32 - 1.0) / 2.0;
        let cy = (self.h as f32 - 1.0) / 2.0;
        for y in 0..self.h {
            for x in 0..self.w {
                let dx = x as f32 - cx;
                let dy = y as f32 - cy;
                // Расстояние (в клетках) внутрь от края раны; ≤ 0 — снаружи.
                let inside = match shape {
                    WoundShape::Circle { radius } => radius - (dx * dx + dy * dy).sqrt(),
                    WoundShape::Cut { half_length, half_width } => {
                        let d = ((dx / half_length).powi(2) + (dy / half_width).powi(2)).sqrt();
                        (1.0 - d) * half_width.min(half_length)
                    }
                };
                self.damage_cell(y * self.w + x, inside + 0.5, depth_mm, p);
            }
        }
    }

    /// Повреждение кругом радиуса `r` (в клетках) с центром (cx, cy) — например, мышью.
    pub fn injure_disk(&mut self, cx: f32, cy: f32, r: f32, depth_mm: f32, p: &Params) {
        let x0 = (cx - r - 1.0).floor().max(0.0) as usize;
        let y0 = (cy - r - 1.0).floor().max(0.0) as usize;
        let x1 = ((cx + r + 1.0).ceil() as usize).min(self.w.saturating_sub(1));
        let y1 = ((cy + r + 1.0).ceil() as usize).min(self.h.saturating_sub(1));
        for y in y0..=y1 {
            for x in x0..=x1 {
                let (dx, dy) = (x as f32 - cx, y as f32 - cy);
                let inside = r - (dx * dx + dy * dy).sqrt();
                self.damage_cell(y * self.w + x, inside + 0.5, depth_mm, p);
            }
        }
    }

    /// Хирургически иссечь клетку на глубину `depth_mm`.
    pub fn excise(&mut self, i: usize, depth_mm: f32, p: &Params) {
        self.damage_cell(i, 1.0, depth_mm, p);
    }

    /// Повреждает клетку на долю `amount` (0..1): клетки на краю раны
    /// повреждены частично, поэтому граница не ступенчатая.
    /// В неглубоких ранах часть дермы (с сосудами, фибробластами и фолликулами) уцелевает.
    fn damage_cell(&mut self, i: usize, amount: f32, depth_mm: f32, p: &Params) {
        let a = amount.clamp(0.0, 1.0);
        if a <= 0.0 {
            return;
        }
        if a >= 0.5 {
            self.wound_mask[i] = true;
        }
        let d = a * depth_mm;
        self.depth.data[i] = self.depth.data[i].max(d);
        self.depth_max.data[i] = self.depth_max.data[i].max(d);

        let lost = a * (1.0 - p.residual_dermis(depth_mm));
        let keep = 1.0 - lost;
        // Поверхностные ссадины почти не кровоточат: сосуды лежат в дерме.
        let bleed = a * (depth_mm / 0.5).min(1.0);
        self.bleeding.data[i] = self.bleeding.data[i].max(bleed);
        self.debris.data[i] = self.debris.data[i].max(a * (0.3 + depth_mm / 2.0).min(1.0));
        self.bacteria.data[i] = self.bacteria.data[i].max(p.initial_bacteria * a);
        self.clot.data[i] *= 1.0 - a;
        self.epithelium.data[i] *= 1.0 - a;
        self.fibroblasts.data[i] *= keep;
        self.collagen.data[i] *= keep;
        self.vessels.data[i] *= keep;
    }

    pub fn bacteria_total(&self, i: usize) -> f32 {
        self.bacteria.data[i] + self.bacteria_res.data[i]
    }

    /// Прочность на разрыв относительно здоровой кожи (0..1).
    /// Рубец упирается в ~0.84 из-за потолка зрелости коллагена.
    pub fn strength(&self, i: usize) -> f32 {
        let dead = self.slough.data[i] / (self.slough.data[i] + 0.3);
        self.collagen.data[i] * (0.2 + 0.8 * self.maturity.data[i]) * (1.0 - dead)
    }

    /// Интегральная «целостность» ткани: барьер + матрикс + кровоснабжение,
    /// минус незаполненная полость.
    pub fn integrity(&self, i: usize) -> f32 {
        let cavity = 1.0 / (1.0 + self.depth.data[i]);
        (0.35 * self.epithelium.data[i] + 0.45 * self.strength(i) + 0.2 * self.vessels.data[i]) * cavity
    }
}
