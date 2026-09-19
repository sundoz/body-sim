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
    /// Бактериальные экзотоксины.
    pub toxin: Field,

    // Что выросло на месте утраченных глубоких слоёв (доли от утраченного объёма слоя)
    /// Новые мышечные волокна.
    pub myo: Field,
    /// Фиброз (рубец) в мышце.
    pub muscle_scar: Field,
    /// Вернувшаяся жировая ткань; остальное — фиброзный рубец.
    pub fat_new: Field,

    pub antiseptic_agent: Antiseptic,
    /// Концентрации антибиотиков в плазме (по Antibiotic::index).
    pub abx_plasma: [f32; 3],
    /// Включённые методы очищения от некроза (по Debrider::index).
    pub debriders: [bool; 3],

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
            toxin: f(0.0),
            myo: f(0.0),
            muscle_scar: f(0.0),
            fat_new: f(0.0),
            antiseptic_agent: Antiseptic::Octenidine,
            abx_plasma: [0.0; 3],
            debriders: [false; 3],
            wound_mask: vec![false; w * h],
        }
    }

    pub fn len(&self) -> usize {
        self.w * self.h
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
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
                        // Вырожденный разрез (нулевая полуось) иначе даёт деление на ноль и NaN.
                        let (hl, hw) = (half_length.max(1e-3), half_width.max(1e-3));
                        let d = ((dx / hl).powi(2) + (dy / hw).powi(2)).sqrt();
                        (1.0 - d) * hw.min(hl)
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
    ///
    /// Глубина обрезается по `max_depth_mm`: шаг модели всё равно держит полость в этих
    /// пределах, а незажатая `depth_max` потом даёт нефизичные «утрачено N мм мышцы».
    /// Нефинитная доля (вырожденная форма раны) игнорируется, иначе NaN расползётся по полям.
    fn damage_cell(&mut self, i: usize, amount: f32, depth_mm: f32, p: &Params) {
        if !amount.is_finite() || !depth_mm.is_finite() {
            return;
        }
        let depth_mm = depth_mm.min(p.max_depth_mm);
        let a = amount.clamp(0.0, 1.0);
        if a <= 0.0 {
            return;
        }
        if a >= 0.5 {
            self.wound_mask[i] = true;
        }
        let d = a * depth_mm;
        self.depth.data[i] = self.depth.data[i].max(d);
        let old = self.depth_max.data[i];
        let new = old.max(d);
        self.depth_max.data[i] = new;
        // Восстановленная часть глубоких слоёв теперь приходится на больший объём дефекта.
        let (mo, mn) = (p.lost_muscle_mm(old), p.lost_muscle_mm(new));
        if mn > mo {
            self.myo.data[i] *= mo / mn;
            self.muscle_scar.data[i] *= mo / mn;
        }
        let (fo, fnew) = (p.lost_fat_mm(old), p.lost_fat_mm(new));
        if fnew > fo {
            self.fat_new.data[i] *= fo / fnew;
        }

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

#[cfg(test)]
mod tests {
    use super::*;

    const W: usize = 96;
    const H: usize = 48;

    fn center() -> usize {
        (H / 2) * W + W / 2
    }

    #[test]
    fn healthy_tissue_is_intact() {
        let p = Params::default();
        let t = Tissue::healthy(W, H, &p);
        assert!(t.wound_mask.iter().all(|m| !m));
        for i in [0, center(), t.len() - 1] {
            assert!((t.integrity(i) - 1.0).abs() < 1e-5);
            assert!((t.strength(i) - 1.0).abs() < 1e-5);
            assert_eq!(t.bacteria_total(i), 0.0);
        }
    }

    #[test]
    fn circular_wound_has_expected_area_and_depth() {
        let p = Params::default();
        let mut t = Tissue::healthy(W, H, &p);
        let r = 16.0;
        t.injure(WoundShape::Circle { radius: r }, 2.5, &p);
        let cells = t.wound_mask.iter().filter(|m| **m).count() as f32;
        let expected = std::f32::consts::PI * r * r;
        assert!((cells - expected).abs() / expected < 0.05, "{cells} против {expected}");
        let c = center();
        assert!((t.depth.data[c] - 2.5).abs() < 1e-5);
        assert_eq!(t.epithelium.data[c], 0.0);
        assert_eq!(t.bleeding.data[c], 1.0);
        assert_eq!(t.collagen.data[c], 0.0, "полнослойная рана уносит всю дерму");
        assert_eq!(t.depth.data[0], 0.0, "угол участка не задет");
    }

    #[test]
    fn partial_thickness_wound_keeps_part_of_dermis() {
        let p = Params::default();
        let mut t = Tissue::healthy(W, H, &p);
        t.injure(WoundShape::Circle { radius: 10.0 }, 1.1, &p);
        let c = center();
        assert!((t.collagen.data[c] - 0.5).abs() < 1e-4);
        assert!((t.vessels.data[c] - 0.5).abs() < 1e-4);
    }

    #[test]
    fn superficial_abrasion_barely_bleeds() {
        let p = Params::default();
        let mut t = Tissue::healthy(W, H, &p);
        t.injure(WoundShape::Circle { radius: 10.0 }, 0.1, &p);
        assert!(t.bleeding.data[center()] < 0.25);
    }

    #[test]
    fn deeper_injury_rescales_regenerated_fractions() {
        let p = Params::default();
        let mut t = Tissue::healthy(W, H, &p);
        t.injure(WoundShape::Circle { radius: 10.0 }, 9.0, &p);
        let c = center();
        t.myo.data[c] = 0.8;
        t.fat_new.data[c] = 0.5;
        t.injure_disk((W / 2) as f32, (H / 2) as f32, 3.0, 12.0, &p);
        let lost_before = p.lost_muscle_mm(9.0);
        let lost_after = p.lost_muscle_mm(12.0);
        assert!((t.myo.data[c] - 0.8 * lost_before / lost_after).abs() < 1e-5);
        assert!((t.fat_new.data[c] - 0.5).abs() < 1e-5, "вся клетчатка уже была утрачена");
    }

    #[test]
    fn injury_deeper_than_the_block_is_clamped() {
        let p = Params::default();
        let mut t = Tissue::healthy(W, H, &p);
        t.injure(WoundShape::Circle { radius: 10.0 }, 50.0, &p);
        let c = center();
        assert!((t.depth.data[c] - p.max_depth_mm).abs() < 1e-5);
        assert!((t.depth_max.data[c] - p.max_depth_mm).abs() < 1e-5);
        // Иначе метрики отрапортовали бы утрату мышцы толще самой колонки ткани.
        assert!(p.lost_muscle_mm(t.depth_max.data[c]) < p.max_depth_mm);
    }

    #[test]
    fn degenerate_wound_shapes_do_not_poison_the_fields() {
        let p = Params::default();
        let mut t = Tissue::healthy(W, H, &p);
        t.injure(WoundShape::Cut { half_length: 0.0, half_width: 0.0 }, 2.5, &p);
        t.injure_disk(10.0, 10.0, 3.0, f32::NAN, &p);
        for f in [&t.depth, &t.collagen, &t.epithelium, &t.vessels, &t.clot, &t.fibroblasts] {
            assert!(f.data.iter().all(|v| v.is_finite()), "NaN расползся по полям");
        }
    }

    #[test]
    fn excision_turns_cell_into_wound() {
        let p = Params::default();
        let mut t = Tissue::healthy(W, H, &p);
        t.excise(5, 2.6, &p);
        assert!(t.wound_mask[5]);
        assert!((t.depth.data[5] - 2.6).abs() < 1e-5);
        assert_eq!(t.epithelium.data[5], 0.0);
    }
}
