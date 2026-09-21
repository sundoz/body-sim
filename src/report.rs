use crate::body::{BodySite, BODY_SURFACE_CM2};
use crate::params::Params;
use crate::tissue::Tissue;

/// Показатели раны: средние по её области плюс площади по всему участку.
#[derive(Clone, Debug, Default)]
pub struct Metrics {
    pub hours: f32,
    pub wound_mm2: f32,
    pub open_fraction: f32,
    pub open_mm2: f32,
    pub bleeding: f32,
    pub clot: f32,
    pub debris: f32,
    pub bacteria: f32,
    pub resistant_fraction: f32,
    pub biofilm: f32,
    pub depth: f32,
    pub depth_max: f32,
    pub necrotic_mm2: f32,
    pub infected_mm2: f32,
    pub neutrophils: f32,
    pub m1: f32,
    pub m2: f32,
    pub fibroblasts: f32,
    pub collagen: f32,
    pub maturity: f32,
    pub vessels: f32,
    pub oxygen: f32,
    pub epithelium: f32,
    pub strength: f32,
    pub integrity: f32,
    /// Суммарная концентрация антибиотиков в плазме (в МПК чувствительного штамма).
    pub abx_plasma: f32,
    /// Площадь раны в начале наблюдения (заполняет `Simulation`) — чтобы видеть, что рана разрастается.
    pub initial_wound_mm2: f32,
    /// По каждому антибиотику (по `Antibiotic::index`).
    pub abx: [f32; 3],
    /// Утраченная толщина мышцы и клетчатки (максимум по ране), мм.
    pub lost_muscle_mm: f32,
    pub lost_fat_mm: f32,
    /// Какая доля утраченной мышцы стала новыми волокнами / рубцом (взвешено по объёму).
    pub muscle_regen: f32,
    pub muscle_fibrosis: f32,
    /// Какая доля утраченной клетчатки вернулась жиром.
    pub fat_regen: f32,

    // Привязка к телу
    pub site: BodySite,
    /// Площадь моделируемого участка кожи, мм².
    pub patch_mm2: f32,
    /// Площадь всего поражения на теле, см² (участок — его представительный кусок).
    pub lesion_cm2: f32,
    /// Какая доля этой области тела поражена, %.
    pub lesion_percent: f32,
    /// Какую долю всей поверхности тела занимает поражение, %.
    pub tbsa_percent: f32,
    /// Какую долю поверхности тела занимает инфицированная ткань, %.
    pub infected_tbsa_percent: f32,
}

/// Инфекция захватила больше половины моделируемого участка — она вышла из-под
/// местного контроля. Доля, а не абсолютная площадь: иначе порог зависел бы от
/// размера сетки. При участке 24×12 мм это прежние 150 мм².
const SEPSIS_PATCH_FRACTION: f32 = 150.0 / 288.0;

/// Обширное поражение даёт системную реакцию независимо от местного контроля.
/// Порог грубый: модель организма сюда не входит, это оценка по площади.
const SEPSIS_TBSA_PERCENT: f32 = 10.0;

/// Порог бактериальной нагрузки, выше которого ткань считается инфицированной.
const INFECTED: f32 = 0.2;
/// Толщина мёртвой ткани, выше которой клетка считается некротизированной (мм).
const NECROTIC: f32 = 0.2;

impl Metrics {
    pub fn measure(t: &Tissue, p: &Params, hours: f32) -> Self {
        let cell_mm = p.cell_mm;
        let mut m = Metrics {
            hours,
            abx: t.abx_plasma,
            abx_plasma: t.abx_plasma.iter().sum(),
            site: p.site,
            patch_mm2: t.len() as f32 * cell_mm * cell_mm,
            ..Default::default()
        };
        let (mut mus_w, mut fat_w) = (0.0f32, 0.0f32);
        let area = cell_mm * cell_mm;
        let mut count = 0usize;
        let mut open = 0usize;
        let (mut b_all, mut br_all) = (0.0f32, 0.0f32);
        for i in 0..t.len() {
            let bt = t.bacteria_total(i);
            b_all += bt;
            br_all += t.bacteria_res.data[i];
            if bt > INFECTED {
                m.infected_mm2 += area;
            }
            if t.slough.data[i] > NECROTIC {
                m.necrotic_mm2 += area;
            }
            let dmax = t.depth_max.data[i];
            m.depth_max = m.depth_max.max(dmax);
            let (lm, lf) = (p.lost_muscle_mm(dmax), p.lost_fat_mm(dmax));
            if lm > 0.0 {
                m.lost_muscle_mm = m.lost_muscle_mm.max(lm);
                m.muscle_regen += t.myo.data[i] * lm;
                m.muscle_fibrosis += t.muscle_scar.data[i] * lm;
                mus_w += lm;
            }
            if lf > 0.0 {
                m.lost_fat_mm = m.lost_fat_mm.max(lf);
                m.fat_regen += t.fat_new.data[i] * lf;
                fat_w += lf;
            }
            if !t.wound_mask[i] {
                continue;
            }
            count += 1;
            if t.epithelium.data[i] < 0.5 {
                open += 1;
            }
            m.bleeding += t.bleeding.data[i];
            m.clot += t.clot.data[i];
            m.debris += t.debris.data[i];
            m.bacteria += bt;
            m.biofilm += t.biofilm.data[i];
            m.depth += t.depth.data[i];
            m.neutrophils += t.neutrophils.data[i];
            m.m1 += t.m1.data[i];
            m.m2 += t.m2.data[i];
            m.fibroblasts += t.fibroblasts.data[i];
            m.collagen += t.collagen.data[i];
            m.maturity += t.maturity.data[i];
            m.vessels += t.vessels.data[i];
            m.oxygen += t.oxygen.data[i];
            m.epithelium += t.epithelium.data[i];
            m.strength += t.strength(i);
            m.integrity += t.integrity(i);
        }
        let k = 1.0 / count.max(1) as f32;
        for x in [
            &mut m.bleeding,
            &mut m.clot,
            &mut m.debris,
            &mut m.bacteria,
            &mut m.biofilm,
            &mut m.depth,
            &mut m.neutrophils,
            &mut m.m1,
            &mut m.m2,
            &mut m.fibroblasts,
            &mut m.collagen,
            &mut m.maturity,
            &mut m.vessels,
            &mut m.oxygen,
            &mut m.epithelium,
            &mut m.strength,
            &mut m.integrity,
        ] {
            *x *= k;
        }
        if mus_w > 0.0 {
            m.muscle_regen /= mus_w;
            m.muscle_fibrosis /= mus_w;
        }
        if fat_w > 0.0 {
            m.fat_regen /= fat_w;
        }
        m.wound_mm2 = count as f32 * area;
        m.open_fraction = open as f32 * k;
        m.open_mm2 = open as f32 * area;
        // Долю устойчивых показываем, только когда популяция заметна.
        m.resistant_fraction = if b_all > 5.0 { br_all / b_all } else { 0.0 };
        m.measure_body(p);
        m
    }

    /// Пересчитать площади с участка на всё тело. Моделируемый участок — представительный
    /// кусок поражения: если задана доля поражённой области, доли по участку переносятся
    /// на её настоящую площадь; иначе поражение считается равным самому участку.
    fn measure_body(&mut self, p: &Params) {
        let region_mm2 = self.site.region_cm2() * 100.0;
        let lesion_mm2 = match p.lesion_percent {
            Some(pct) => region_mm2 * pct.clamp(0.0, 100.0) / 100.0,
            None => self.wound_mm2,
        };
        self.lesion_cm2 = lesion_mm2 / 100.0;
        self.lesion_percent = if region_mm2 > 0.0 { lesion_mm2 / region_mm2 * 100.0 } else { 0.0 };
        self.tbsa_percent = lesion_mm2 / (BODY_SURFACE_CM2 * 100.0) * 100.0;
        // Какая доля поражения инфицирована — меряется по участку, переносится на его площадь.
        let infected_share = if self.patch_mm2 > 0.0 { self.infected_mm2 / self.patch_mm2 } else { 0.0 };
        self.infected_tbsa_percent = self.tbsa_percent * infected_share.min(1.0);
    }

    pub fn phase(&self) -> Phase {
        // На нетронутом участке все средние по ране равны нулю (делить не на что),
        // и «нулевой коллаген» ложно выглядел бы как пролиферация.
        if self.wound_mm2 <= 0.0 {
            return Phase::Remodeling;
        }
        if self.bleeding > 0.05 {
            Phase::Hemostasis
        } else if self.neutrophils + self.m1 > self.m2 + self.fibroblasts {
            Phase::Inflammation
        } else if self.open_fraction > 0.01 || self.collagen < 0.5 {
            Phase::Proliferation
        } else {
            Phase::Remodeling
        }
    }

    /// Клиническая оценка состояния раны.
    pub fn condition(&self) -> Condition {
        if self.infected_mm2 > SEPSIS_PATCH_FRACTION * self.patch_mm2
            || self.infected_tbsa_percent > SEPSIS_TBSA_PERCENT
        {
            Condition::Sepsis
        } else if self.infected_mm2 > 1.2 * self.wound_mm2 + 10.0
            || (self.wound_mm2 > 1.5 * self.initial_wound_mm2 + 10.0 && self.infected_mm2 > 0.5 * self.wound_mm2)
        {
            Condition::Spreading
        } else if self.necrotic_mm2 > 3.0 && self.necrotic_mm2 > 0.25 * self.wound_mm2 {
            Condition::Necrosis
        } else if self.open_fraction < 0.01 {
            Condition::Healed
        } else if self.hours > 30.0 * 24.0 || self.biofilm > 0.5 {
            Condition::Chronic
        } else {
            Condition::Healing
        }
    }

    pub const CSV_HEADER: &'static str = "hours,day,phase,condition,site,patch_mm2,lesion_cm2,lesion_percent,tbsa_percent,infected_tbsa_percent,wound_mm2,open_fraction,open_mm2,depth,depth_max,necrotic_mm2,infected_mm2,bacteria,resistant_fraction,biofilm,abx_plasma,lost_muscle_mm,muscle_regen,muscle_fibrosis,lost_fat_mm,fat_regen,bleeding,clot,debris,neutrophils,m1,m2,fibroblasts,collagen,maturity,vessels,oxygen,epithelium,strength,integrity";

    pub fn csv_row(&self) -> String {
        format!(
            "{:.1},{:.3},{},{},{},{:.2},{:.3},{:.3},{:.4},{:.4},{:.2},{:.4},{:.3},{:.3},{:.3},{:.2},{:.2},{:.4},{:.4},{:.4},{:.3},{:.2},{:.4},{:.4},{:.2},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4}",
            self.hours,
            self.hours / 24.0,
            self.phase().key(),
            self.condition().key(),
            self.site.key(),
            self.patch_mm2,
            self.lesion_cm2,
            self.lesion_percent,
            self.tbsa_percent,
            self.infected_tbsa_percent,
            self.wound_mm2,
            self.open_fraction,
            self.open_mm2,
            self.depth,
            self.depth_max,
            self.necrotic_mm2,
            self.infected_mm2,
            self.bacteria,
            self.resistant_fraction,
            self.biofilm,
            self.abx_plasma,
            self.lost_muscle_mm,
            self.muscle_regen,
            self.muscle_fibrosis,
            self.lost_fat_mm,
            self.fat_regen,
            self.bleeding,
            self.clot,
            self.debris,
            self.neutrophils,
            self.m1,
            self.m2,
            self.fibroblasts,
            self.collagen,
            self.maturity,
            self.vessels,
            self.oxygen,
            self.epithelium,
            self.strength,
            self.integrity,
        )
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Phase {
    Hemostasis,
    Inflammation,
    Proliferation,
    Remodeling,
}

impl Phase {
    pub fn title(self) -> &'static str {
        match self {
            Self::Hemostasis => "гемостаз",
            Self::Inflammation => "воспаление",
            Self::Proliferation => "пролиферация",
            Self::Remodeling => "ремоделирование",
        }
    }

    pub fn key(self) -> &'static str {
        match self {
            Self::Hemostasis => "hemostasis",
            Self::Inflammation => "inflammation",
            Self::Proliferation => "proliferation",
            Self::Remodeling => "remodeling",
        }
    }
}

/// Состояние раны — от благополучного к угрожающему жизни.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub enum Condition {
    Healed,
    Healing,
    Chronic,
    Necrosis,
    Spreading,
    Sepsis,
}

impl Condition {
    pub fn title(self) -> &'static str {
        match self {
            Self::Healed => "зажила",
            Self::Healing => "заживает",
            Self::Chronic => "хроническая рана",
            Self::Necrosis => "некроз",
            Self::Spreading => "распространяющаяся инфекция",
            Self::Sepsis => "генерализация инфекции (сепсис)",
        }
    }

    pub fn key(self) -> &'static str {
        match self {
            Self::Healed => "healed",
            Self::Healing => "healing",
            Self::Chronic => "chronic",
            Self::Necrosis => "necrosis",
            Self::Spreading => "spreading",
            Self::Sepsis => "sepsis",
        }
    }
}

/// Какое поле рисовать на карте.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum View {
    Integrity,
    Strength,
    Epithelium,
    Collagen,
    Vessels,
    Oxygen,
    Clot,
    Depth,
    Necrosis,
    Bacteria,
    Resistant,
    Biofilm,
    Toxin,
    Neutrophils,
    Macrophages,
    Fibroblasts,
    Antiseptic,
    Antibiotic,
}

impl View {
    pub const NAMES: &'static str = "integrity, strength, epithelium, collagen, vessels, oxygen, clot, depth, necrosis, bacteria, resistant, biofilm, toxin, neutrophils, macrophages, fibroblasts, antiseptic, antibiotic";

    pub const ALL: [View; 18] = [
        Self::Integrity,
        Self::Strength,
        Self::Epithelium,
        Self::Collagen,
        Self::Vessels,
        Self::Oxygen,
        Self::Clot,
        Self::Depth,
        Self::Necrosis,
        Self::Bacteria,
        Self::Resistant,
        Self::Biofilm,
        Self::Toxin,
        Self::Neutrophils,
        Self::Macrophages,
        Self::Fibroblasts,
        Self::Antiseptic,
        Self::Antibiotic,
    ];

    pub fn parse(s: &str) -> Option<Self> {
        Some(match s {
            "integrity" => Self::Integrity,
            "strength" => Self::Strength,
            "epithelium" => Self::Epithelium,
            "collagen" => Self::Collagen,
            "vessels" => Self::Vessels,
            "oxygen" => Self::Oxygen,
            "clot" => Self::Clot,
            "depth" => Self::Depth,
            "necrosis" => Self::Necrosis,
            "bacteria" => Self::Bacteria,
            "resistant" => Self::Resistant,
            "biofilm" => Self::Biofilm,
            "toxin" => Self::Toxin,
            "neutrophils" => Self::Neutrophils,
            "macrophages" => Self::Macrophages,
            "fibroblasts" => Self::Fibroblasts,
            "antiseptic" => Self::Antiseptic,
            "antibiotic" => Self::Antibiotic,
            _ => return None,
        })
    }

    pub fn title(self) -> &'static str {
        match self {
            Self::Integrity => "Целостность",
            Self::Strength => "Прочность",
            Self::Epithelium => "Эпителий",
            Self::Collagen => "Коллаген",
            Self::Vessels => "Сосуды",
            Self::Oxygen => "Кислород",
            Self::Clot => "Сгусток",
            Self::Depth => "Глубина",
            Self::Necrosis => "Некроз",
            Self::Bacteria => "Бактерии",
            Self::Resistant => "Устойчивые",
            Self::Biofilm => "Биоплёнка",
            Self::Toxin => "Токсины",
            Self::Neutrophils => "Нейтрофилы",
            Self::Macrophages => "Макрофаги",
            Self::Fibroblasts => "Фибробласты",
            Self::Antiseptic => "Антисептик",
            Self::Antibiotic => "Антибиотик",
        }
    }

    /// Чем 1.0 на шкале является для этого поля.
    pub fn scale_label(self) -> &'static str {
        match self {
            Self::Depth => "5 мм",
            Self::Necrosis => "2 мм",
            Self::Antibiotic => "8 МПК",
            Self::Toxin => "3",
            _ => "1",
        }
    }

    /// Значение поля, нормированное на 0..1 для карты.
    pub fn value(self, t: &Tissue, i: usize) -> f32 {
        match self {
            Self::Integrity => t.integrity(i),
            Self::Strength => t.strength(i),
            Self::Epithelium => t.epithelium.data[i],
            Self::Collagen => t.collagen.data[i],
            Self::Vessels => t.vessels.data[i],
            Self::Oxygen => t.oxygen.data[i],
            Self::Clot => t.clot.data[i],
            Self::Depth => t.depth.data[i] / 5.0,
            Self::Necrosis => t.slough.data[i] / 2.0,
            Self::Bacteria => t.bacteria_total(i),
            Self::Resistant => t.bacteria_res.data[i],
            Self::Biofilm => t.biofilm.data[i],
            Self::Toxin => t.toxin.data[i] / 3.0,
            Self::Neutrophils => t.neutrophils.data[i],
            Self::Macrophages => t.m1.data[i] + t.m2.data[i],
            Self::Fibroblasts => t.fibroblasts.data[i],
            Self::Antiseptic => t.antiseptic.data[i],
            Self::Antibiotic => t.antibiotic.data[i] / 8.0,
        }
    }
}

const RAMP: &[u8] = b" .:-=+*#%@";

/// ASCII-карта поля; две строки сетки на строку терминала (символы вытянуты по вертикали).
pub fn render_map(t: &Tissue, view: View) -> String {
    let mut out = String::new();
    out.push('+');
    out.push_str(&"-".repeat(t.w));
    out.push_str("+\n");
    for y in (0..t.h).step_by(2) {
        out.push('|');
        for x in 0..t.w {
            let i0 = y * t.w + x;
            let v = if y + 1 < t.h { 0.5 * (view.value(t, i0) + view.value(t, i0 + t.w)) } else { view.value(t, i0) };
            let k = (v.clamp(0.0, 1.0) * (RAMP.len() - 1) as f32).round() as usize;
            out.push(RAMP[k] as char);
        }
        out.push_str("|\n");
    }
    out.push('+');
    out.push_str(&"-".repeat(t.w));
    out.push('+');
    out
}

pub fn print_snapshot(m: &Metrics, map: Option<&str>) {
    println!();
    println!(
        "=== День {:>5.1} · фаза: {} · состояние: {} ===",
        m.hours / 24.0,
        m.phase().title(),
        m.condition().title()
    );
    println!(
        "  {} · поражено {:.1} см² ({:.1}% области, {:.2}% поверхности тела), инфицировано {:.2}% тела",
        m.site.title(),
        m.lesion_cm2,
        m.lesion_percent,
        m.tbsa_percent,
        m.infected_tbsa_percent
    );
    println!(
        "  рана {:>5.1} мм², открыто {:>5.1} мм² ({:>3.0}%) | глубина {:.1} мм (макс {:.1}) | некроз {:.1} мм² | инфицировано {:.1} мм²",
        m.wound_mm2,
        m.open_mm2,
        m.open_fraction * 100.0,
        m.depth,
        m.depth_max,
        m.necrotic_mm2,
        m.infected_mm2
    );
    println!(
        "  бактерии {:.3} (устойч. {:.0}%) биоплёнка {:.2} | антибиотик в плазме {:.1} | кровотеч. {:.2} сгусток {:.2}",
        m.bacteria,
        m.resistant_fraction * 100.0,
        m.biofilm,
        m.abx_plasma,
        m.bleeding,
        m.clot
    );
    println!(
        "  нейтрофилы {:.2} | макрофаги M1 {:.2} M2 {:.2} | фибробласты {:.2} | сосуды {:.2} | O₂ {:.2}",
        m.neutrophils, m.m1, m.m2, m.fibroblasts, m.vessels, m.oxygen
    );
    println!(
        "  эпителий {:.2} | коллаген {:.2} (зрелость {:.2}) | прочность {:.2} | целостность {:.2}",
        m.epithelium, m.collagen, m.maturity, m.strength, m.integrity
    );
    if m.lost_fat_mm > 0.0 {
        let mut line =
            format!("  клетчатка: утрачено {:.1} мм, вернулось жиром {:.0}%", m.lost_fat_mm, m.fat_regen * 100.0);
        if m.lost_muscle_mm > 0.0 {
            line += &format!(
                " | мышца: утрачено {:.1} мм, новые волокна {:.0}%, фиброз {:.0}%",
                m.lost_muscle_mm,
                m.muscle_regen * 100.0,
                m.muscle_fibrosis * 100.0
            );
        }
        println!("{line}");
    }
    if let Some(map) = map {
        println!("{map}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tissue::WoundShape;

    #[test]
    fn healthy_skin_has_no_wound() {
        let p = Params::default();
        let t = Tissue::healthy(96, 48, &p);
        let m = Metrics::measure(&t, &p, 0.0);
        assert_eq!(m.wound_mm2, 0.0);
        assert_eq!(m.infected_mm2, 0.0);
        assert_eq!(m.necrotic_mm2, 0.0);
        assert_eq!(m.lost_fat_mm, 0.0);
        assert_eq!(m.condition(), Condition::Healed);
    }

    #[test]
    fn fresh_wound_metrics() {
        let p = Params::default();
        let mut t = Tissue::healthy(96, 48, &p);
        t.injure(WoundShape::Circle { radius: 16.0 }, 9.0, &p);
        let m = Metrics::measure(&t, &p, 0.0);
        let area = std::f32::consts::PI * 4.0 * 4.0;
        assert!((m.wound_mm2 - area).abs() / area < 0.05);
        assert!((m.open_fraction - 1.0).abs() < 1e-6);
        assert!((m.depth_max - 9.0).abs() < 1e-5);
        assert!((m.lost_fat_mm - 6.0).abs() < 1e-4);
        assert!((m.lost_muscle_mm - 0.6).abs() < 1e-4);
        assert_eq!(m.muscle_regen, 0.0);
        assert_eq!(m.phase(), Phase::Hemostasis);
        assert_eq!(m.condition(), Condition::Healing);
    }

    #[test]
    fn condition_classification() {
        let base = Metrics {
            wound_mm2: 50.0,
            initial_wound_mm2: 50.0,
            open_fraction: 0.8,
            hours: 24.0,
            patch_mm2: 288.0,
            ..Default::default()
        };
        let with = |f: &dyn Fn(&mut Metrics)| {
            let mut m = base.clone();
            f(&mut m);
            m.condition()
        };
        assert_eq!(with(&|_| {}), Condition::Healing);
        assert_eq!(with(&|m| m.infected_mm2 = 200.0), Condition::Sepsis);
        assert_eq!(with(&|m| m.infected_mm2 = 80.0), Condition::Spreading);
        assert_eq!(
            with(&|m| m.infected_mm2 = 50.0),
            Condition::Healing,
            "инфекция в пределах раны — не распространение"
        );
        assert_eq!(with(&|m| m.necrotic_mm2 = 20.0), Condition::Necrosis);
        assert_eq!(
            with(&|m| {
                m.wound_mm2 = 120.0;
                m.infected_mm2 = 110.0;
            }),
            Condition::Spreading,
            "рана разрастается за счёт инфекции"
        );
        assert_eq!(
            with(&|m| {
                m.wound_mm2 = 120.0;
                m.infected_mm2 = 5.0;
            }),
            Condition::Healing,
            "рану расширила хирургия, инфекции уже нет"
        );
        assert_eq!(with(&|m| m.open_fraction = 0.0), Condition::Healed);
        assert_eq!(with(&|m| m.hours = 31.0 * 24.0), Condition::Chronic);
        assert_eq!(with(&|m| m.biofilm = 0.6), Condition::Chronic);
        assert!(Condition::Sepsis > Condition::Necrosis && Condition::Necrosis > Condition::Healing);
        // Обширное поражение даёт системную реакцию даже при умеренной доле инфекции на участке.
        assert_eq!(with(&|m| m.infected_tbsa_percent = 12.0), Condition::Sepsis);
    }

    #[test]
    fn severity_threshold_does_not_depend_on_grid_size() {
        // Раньше порог был абсолютным (150 мм²), и на вдвое большем участке
        // та же доля инфицированной ткани давала другой диагноз.
        let at = |patch: f32| {
            Metrics {
                wound_mm2: 50.0,
                initial_wound_mm2: 50.0,
                open_fraction: 0.8,
                hours: 24.0,
                patch_mm2: patch,
                infected_mm2: 0.6 * patch,
                ..Default::default()
            }
            .condition()
        };
        assert_eq!(at(288.0), Condition::Sepsis);
        assert_eq!(at(2400.0), at(288.0), "диагноз не должен зависеть от размера сетки");
    }

    #[test]
    fn body_areas_scale_from_the_patch_to_the_whole_body() {
        let p = Params { site: BodySite::Back, lesion_percent: Some(50.0), ..Params::default() };
        let mut t = Tissue::healthy(96, 48, &p);
        t.injure(crate::tissue::WoundShape::Circle { radius: 16.0 }, 2.5, &p);
        let m = Metrics::measure(&t, &p, 0.0);
        assert!((m.patch_mm2 - 288.0).abs() < 1e-3, "участок 24×12 мм");
        // Половина спины — это 6.5% поверхности тела.
        assert!((m.lesion_percent - 50.0).abs() < 1e-3);
        assert!((m.tbsa_percent - 6.5).abs() < 0.01, "{}", m.tbsa_percent);
        assert!((m.lesion_cm2 - BodySite::Back.region_cm2() / 2.0).abs() < 1.0);
        // Без указанной доли поражение равно самому участку.
        let p = Params { site: BodySite::Back, ..Params::default() };
        let m = Metrics::measure(&t, &p, 0.0);
        assert!((m.lesion_cm2 - m.wound_mm2 / 100.0).abs() < 1e-4);
        assert!(m.tbsa_percent < 0.01, "рана Ø8 мм — ничтожная доля поверхности тела");
    }

    #[test]
    fn intact_skin_is_not_in_the_proliferation_phase() {
        let p = Params::default();
        let t = Tissue::healthy(96, 48, &p);
        let m = Metrics::measure(&t, &p, 0.0);
        assert_eq!(m.wound_mm2, 0.0);
        assert_eq!(m.phase(), Phase::Remodeling, "раны нет — заживать нечему");
    }

    #[test]
    fn phases() {
        let m = |f: &dyn Fn(&mut Metrics)| {
            let mut m = Metrics { wound_mm2: 50.0, open_fraction: 0.5, collagen: 0.2, ..Default::default() };
            f(&mut m);
            m.phase()
        };
        assert_eq!(m(&|m| m.bleeding = 0.5), Phase::Hemostasis);
        assert_eq!(m(&|m| m.neutrophils = 0.5), Phase::Inflammation);
        assert_eq!(m(&|m| m.fibroblasts = 0.5), Phase::Proliferation);
        assert_eq!(
            m(&|m| {
                m.open_fraction = 0.0;
                m.collagen = 0.8;
            }),
            Phase::Remodeling
        );
    }

    #[test]
    fn csv_row_matches_header() {
        let row = Metrics::default().csv_row();
        assert_eq!(row.split(',').count(), Metrics::CSV_HEADER.split(',').count());
    }

    #[test]
    fn every_view_parses_and_is_normalized() {
        let names: Vec<&str> = View::NAMES.split(", ").collect();
        assert_eq!(names.len(), View::ALL.len());
        let p = Params::default();
        let t = Tissue::healthy(8, 8, &p);
        for (name, v) in names.iter().zip(View::ALL) {
            assert_eq!(View::parse(name), Some(v), "{name}");
            assert!(!v.title().is_empty());
            let x = v.value(&t, 0);
            assert!((0.0..=1.0).contains(&x), "{name}: {x}");
        }
    }

    #[test]
    fn resistant_share_is_hidden_for_negligible_population() {
        let p = Params::default();
        let mut t = Tissue::healthy(96, 48, &p);
        t.bacteria_res.data[0] = 1e-4;
        assert_eq!(Metrics::measure(&t, &p, 0.0).resistant_fraction, 0.0);
        for i in 0..200 {
            t.bacteria.data[i] = 0.5;
            t.bacteria_res.data[i] = 0.5;
        }
        assert!((Metrics::measure(&t, &p, 0.0).resistant_fraction - 0.5).abs() < 0.01);
    }

    #[test]
    fn ascii_map_has_frame() {
        let p = Params::default();
        let t = Tissue::healthy(10, 4, &p);
        let map = render_map(&t, View::Integrity);
        assert_eq!(map.lines().count(), 4, "две строки сетки на строку терминала + рамка");
        assert!(map.lines().nth(1).unwrap().contains("@@@@"));
    }
}
