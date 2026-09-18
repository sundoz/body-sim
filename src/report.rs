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
    pub abx_plasma: f32,
}

/// Порог бактериальной нагрузки, выше которого ткань считается инфицированной.
const INFECTED: f32 = 0.2;
/// Толщина мёртвой ткани, выше которой клетка считается некротизированной (мм).
const NECROTIC: f32 = 0.2;

impl Metrics {
    pub fn measure(t: &Tissue, hours: f32, cell_mm: f32) -> Self {
        let mut m = Metrics { hours, abx_plasma: t.abx_plasma, ..Default::default() };
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
            m.depth_max = m.depth_max.max(t.depth_max.data[i]);
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
        m.wound_mm2 = count as f32 * area;
        m.open_fraction = open as f32 * k;
        m.open_mm2 = open as f32 * area;
        // Долю устойчивых показываем, только когда популяция заметна.
        m.resistant_fraction = if b_all > 5.0 { br_all / b_all } else { 0.0 };
        m
    }

    pub fn phase(&self) -> Phase {
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
        if self.infected_mm2 > 150.0 {
            Condition::Sepsis
        } else if self.infected_mm2 > 1.2 * self.wound_mm2 + 10.0 {
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

    pub const CSV_HEADER: &'static str = "hours,day,phase,condition,wound_mm2,open_fraction,open_mm2,depth,depth_max,necrotic_mm2,infected_mm2,bacteria,resistant_fraction,biofilm,abx_plasma,bleeding,clot,debris,neutrophils,m1,m2,fibroblasts,collagen,maturity,vessels,oxygen,epithelium,strength,integrity";

    pub fn csv_row(&self) -> String {
        format!(
            "{:.1},{:.3},{},{},{:.2},{:.4},{:.3},{:.3},{:.3},{:.2},{:.2},{:.4},{:.4},{:.4},{:.3},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4}",
            self.hours,
            self.hours / 24.0,
            self.phase().key(),
            self.condition().key(),
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
    Neutrophils,
    Macrophages,
    Fibroblasts,
    Antiseptic,
    Antibiotic,
}

impl View {
    pub const NAMES: &'static str = "integrity, strength, epithelium, collagen, vessels, oxygen, clot, depth, necrosis, bacteria, resistant, biofilm, neutrophils, macrophages, fibroblasts, antiseptic, antibiotic";

    pub const ALL: [View; 17] = [
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
            let v = if y + 1 < t.h {
                0.5 * (view.value(t, i0) + view.value(t, i0 + t.w))
            } else {
                view.value(t, i0)
            };
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
    if let Some(map) = map {
        println!("{map}");
    }
}
