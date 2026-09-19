//! Лечение: местные антисептики, системные антибиотики, очищение раны от некроза.
//!
//! Качественные ориентиры (in vitro и клинические обзоры):
//! - Токсичность для клеток раны: гипохлорит < полигексанид ≈ октенидин < хлоргексидин ≈ йод < перекись.
//!   Повидон-йод сильнее всех бьёт по фибробластам, хлоргексидин — по кератиноцитам и их миграции.
//! - Гипохлорит, перекись и отчасти йод инактивируются органикой (гной, кровь, некроз).
//! - Системный антибиотик приходит с кровью: в некроз и ишемизированную ткань почти не попадает.
//!   Клиндамицин проникает лучше и подавляет выработку токсинов — поэтому его добавляют при
//!   некротизирующей инфекции. Ванкомицин действует медленно и хуже проникает в ткани,
//!   но убивает штамм, устойчивый к β-лактамам.
//! - Некроз убирают: хирургически (сразу), личинками (1–2 недели), коллагеназой (дни–недели),
//!   гидрогелем (аутолиз, медленнее всего).

use crate::params::Params;
use crate::tissue::Tissue;

// ---------------------------------------------------------------- антисептики

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Antiseptic {
    Octenidine,
    Polyhexanide,
    Hypochlorous,
    Chlorhexidine,
    PovidoneIodine,
    Peroxide,
}

/// Свойства антисептика при концентрации 1.0 (сразу после нанесения).
#[derive(Clone, Copy, Debug)]
pub struct AntisepticProps {
    /// Скорость гибели бактерий, 1/ч.
    pub kill: f32,
    /// Какая доля действия проникает в биоплёнку.
    pub biofilm_pen: f32,
    /// Какая доля действия проникает в слаф и струп.
    pub necro_pen: f32,
    /// Насколько гасится органикой (гной, некроз): 0 — не гасится, 1 — полностью.
    pub organic: f32,
    /// Разрушение матрикса биоплёнки, 1/ч.
    pub disrupt: f32,
    /// Токсичность для кератиноцитов (эпителизация), 1/ч.
    pub tox_epi: f32,
    /// Токсичность для фибробластов, эндотелия и регенерирующих миобластов, 1/ч.
    pub tox_fib: f32,
    pub half_life_h: f32,
}

impl Antiseptic {
    pub const ALL: [Antiseptic; 6] = [
        Self::Octenidine,
        Self::Polyhexanide,
        Self::Hypochlorous,
        Self::Chlorhexidine,
        Self::PovidoneIodine,
        Self::Peroxide,
    ];

    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "octenidine" => Some(Self::Octenidine),
            "polyhexanide" | "phmb" => Some(Self::Polyhexanide),
            "hypochlorous" | "hocl" => Some(Self::Hypochlorous),
            "chlorhexidine" => Some(Self::Chlorhexidine),
            "povidone" | "iodine" => Some(Self::PovidoneIodine),
            "peroxide" => Some(Self::Peroxide),
            _ => None,
        }
    }

    pub fn title(self) -> &'static str {
        match self {
            Self::Octenidine => "Октенидин",
            Self::Polyhexanide => "Полигексанид",
            Self::Hypochlorous => "Гипохлорит",
            Self::Chlorhexidine => "Хлоргексидин",
            Self::PovidoneIodine => "Повидон-йод",
            Self::Peroxide => "Перекись H₂O₂",
        }
    }

    pub fn props(self) -> AntisepticProps {
        let p = |kill, biofilm_pen, necro_pen, organic, disrupt, tox_epi, tox_fib, half_life_h| AntisepticProps {
            kill,
            biofilm_pen,
            necro_pen,
            organic,
            disrupt,
            tox_epi,
            tox_fib,
            half_life_h,
        };
        // Токсичность — относительно скорости обновления клеток после калибровки (эпителий ~0.056/ч,
        // фибробласты ~0.045/ч): порядок средств — по сравнительным in vitro данным.
        match self {
            // Высокий индекс биосовместимости, остаётся на ткани несколько часов.
            Self::Octenidine => p(3.0, 0.6, 0.4, 0.2, 0.08, 0.0075, 0.01, 6.0),
            Self::Polyhexanide => p(2.2, 0.6, 0.4, 0.2, 0.06, 0.005, 0.006, 8.0),
            // Почти нетоксичен, но быстро расходуется на органике.
            Self::Hypochlorous => p(2.0, 0.5, 0.3, 0.7, 0.1, 0.0025, 0.0025, 0.5),
            // Связывается с белками и держится долго; тормозит миграцию кератиноцитов.
            Self::Chlorhexidine => p(2.0, 0.3, 0.3, 0.5, 0.02, 0.015, 0.0125, 6.0),
            // Широкий спектр, хорошо проникает, но самый токсичный для фибробластов.
            Self::PovidoneIodine => p(3.0, 0.5, 0.5, 0.5, 0.05, 0.015, 0.0225, 3.0),
            // В основном «вспенивает» и механически чистит; бактерицидность слабая, токсичность высокая.
            Self::Peroxide => p(1.0, 0.2, 0.2, 0.8, 0.1, 0.06, 0.06, 0.3),
        }
    }
}

// ---------------------------------------------------------------- антибиотики

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Antibiotic {
    Cefazolin,
    Clindamycin,
    Vancomycin,
}

/// Фармакокинетика и фармакодинамика антибиотика (концентрации — в МПК чувствительного штамма).
#[derive(Clone, Copy, Debug)]
pub struct AntibioticProps {
    pub half_life_h: f32,
    /// Прирост концентрации в плазме на одну дозу.
    pub dose: f32,
    /// Максимальная скорость гибели бактерий, 1/ч (бактериостатики — ниже).
    pub emax: f32,
    pub mic_sensitive: f32,
    /// МПК штамма, устойчивого к β-лактамам.
    pub mic_resistant: f32,
    /// Проникновение в ткани относительно плазмы.
    pub tissue_pen: f32,
    /// Какая доля доходит сквозь слаф и струп.
    pub necro_pen: f32,
    /// Насколько биоплёнка защищает от препарата.
    pub biofilm_protect: f32,
    /// Подавление выработки токсинов (0..1) при концентрации выше МПК.
    pub antitoxin: f32,
}

impl Antibiotic {
    pub const ALL: [Antibiotic; 3] = [Self::Cefazolin, Self::Clindamycin, Self::Vancomycin];

    pub fn index(self) -> usize {
        self as usize
    }

    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "cefazolin" => Some(Self::Cefazolin),
            "clindamycin" => Some(Self::Clindamycin),
            "vancomycin" => Some(Self::Vancomycin),
            _ => None,
        }
    }

    pub fn title(self) -> &'static str {
        match self {
            Self::Cefazolin => "Цефазолин",
            Self::Clindamycin => "Клиндамицин",
            Self::Vancomycin => "Ванкомицин",
        }
    }

    pub fn props(self) -> AntibioticProps {
        match self {
            // β-лактам: быстро убивает, короткий t½, на устойчивый штамм действует слабо.
            Self::Cefazolin => AntibioticProps {
                half_life_h: 2.0,
                dose: 8.0,
                emax: 0.6,
                mic_sensitive: 1.0,
                mic_resistant: 4.0,
                tissue_pen: 1.0,
                necro_pen: 0.0,
                biofilm_protect: 0.9,
                antitoxin: 0.0,
            },
            // Бактериостатик: убивает медленнее, зато хорошо проникает (в том числе в некроз)
            // и выключает синтез токсинов.
            Self::Clindamycin => AntibioticProps {
                half_life_h: 3.0,
                dose: 6.0,
                emax: 0.25,
                mic_sensitive: 1.0,
                mic_resistant: 1.5,
                tissue_pen: 1.6,
                necro_pen: 0.5,
                biofilm_protect: 0.7,
                antitoxin: 0.7,
            },
            // Гликопептид: медленный, плохо проникает в ткани, но устойчивый штамм к нему чувствителен.
            Self::Vancomycin => AntibioticProps {
                half_life_h: 6.0,
                dose: 4.0,
                emax: 0.35,
                mic_sensitive: 1.0,
                mic_resistant: 1.0,
                tissue_pen: 0.5,
                necro_pen: 0.0,
                biofilm_protect: 0.9,
                antitoxin: 0.0,
            },
        }
    }
}

// ---------------------------------------------------------------- очищение от некроза

/// Непрерывные методы очищения раны от некроза (включаются как повязка/курс).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Debrider {
    /// Аутолиз под влажной повязкой: ускоряет собственные ферменты макрофагов, размягчает струп.
    Hydrogel,
    /// Ферментная мазь: переваривает мёртвый коллаген независимо от макрофагов.
    Collagenase,
    /// Личинки Lucilia sericata: быстро и избирательно съедают некроз, выделяют антибактериальные вещества.
    Larvae,
}

impl Debrider {
    pub const ALL: [Debrider; 3] = [Self::Hydrogel, Self::Collagenase, Self::Larvae];

    pub fn index(self) -> usize {
        self as usize
    }

    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "hydrogel" => Some(Self::Hydrogel),
            "collagenase" => Some(Self::Collagenase),
            "larvae" | "maggots" => Some(Self::Larvae),
            _ => None,
        }
    }

    pub fn title(self) -> &'static str {
        match self {
            Self::Hydrogel => "Гидрогель",
            Self::Collagenase => "Коллагеназа",
            Self::Larvae => "Личинки",
        }
    }
}

// ---------------------------------------------------------------- план и журнал

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum EventKind {
    Antiseptic,
    Dose(Antibiotic),
    Debridement,
}

#[derive(Clone, Copy, Debug)]
pub struct Event {
    pub hours: f32,
    pub kind: EventKind,
}

/// План лечения и журнал процедур.
#[derive(Clone, Debug, Default)]
pub struct Therapy {
    pub antiseptic: Option<Antiseptic>,
    pub antiseptic_every_h: f32,
    /// Интервал введения каждого антибиотика (None — не назначен).
    pub antibiotic_every_h: [Option<f32>; 3],
    pub events: Vec<Event>,
    next_antiseptic_h: f32,
    next_dose_h: [f32; 3],
}

impl Therapy {
    /// Регулярные перевязки с антисептиком начиная с `now`.
    pub fn start_antiseptic(&mut self, agent: Antiseptic, every_h: f32, now: f32) {
        self.antiseptic = Some(agent);
        self.antiseptic_every_h = every_h;
        self.next_antiseptic_h = now;
    }

    pub fn stop_antiseptic(&mut self) {
        self.antiseptic = None;
    }

    /// Курс антибиотика: доза каждые `every_h` часов начиная с `now`.
    pub fn start_antibiotic(&mut self, drug: Antibiotic, every_h: f32, now: f32) {
        self.antibiotic_every_h[drug.index()] = Some(every_h);
        self.next_dose_h[drug.index()] = now;
    }

    pub fn stop_antibiotic(&mut self, drug: Antibiotic) {
        self.antibiotic_every_h[drug.index()] = None;
    }

    pub fn apply_antiseptic_now(&mut self, t: &mut Tissue, agent: Antiseptic, now: f32) {
        apply_antiseptic(t, agent);
        self.events.push(Event { hours: now, kind: EventKind::Antiseptic });
    }

    pub fn debride_now(&mut self, t: &mut Tissue, p: &Params, now: f32) {
        debride(t, p);
        self.events.push(Event { hours: now, kind: EventKind::Debridement });
    }

    /// Вызывается перед каждым шагом модели: процедуры по расписанию и фармакокинетика.
    pub fn tick(&mut self, t: &mut Tissue, p: &Params, now: f32) {
        if let Some(agent) = self.antiseptic {
            if now >= self.next_antiseptic_h {
                self.apply_antiseptic_now(t, agent, now);
                self.next_antiseptic_h = now + self.antiseptic_every_h;
            }
        }
        let ln2 = std::f32::consts::LN_2;
        for drug in Antibiotic::ALL {
            let k = drug.index();
            let props = drug.props();
            if let Some(every) = self.antibiotic_every_h[k] {
                if now >= self.next_dose_h[k] {
                    t.abx_plasma[k] += props.dose;
                    self.events.push(Event { hours: now, kind: EventKind::Dose(drug) });
                    self.next_dose_h[k] = now + every;
                }
            }
            t.abx_plasma[k] *= (-p.dt * ln2 / props.half_life_h).exp();
        }
        let as_decay = (-p.dt * ln2 / t.antiseptic_agent.props().half_life_h).exp();
        for a in t.antiseptic.data.iter_mut() {
            *a *= as_decay;
        }
    }
}

/// Нанести антисептик на открытую поверхность раны.
pub fn apply_antiseptic(t: &mut Tissue, agent: Antiseptic) {
    t.antiseptic_agent = agent;
    for i in 0..t.len() {
        let open = t.epithelium.data[i] < 0.8 || t.slough.data[i] > 0.05;
        if t.wound_mask[i] && open {
            t.antiseptic.data[i] = 1.0;
            if agent == Antiseptic::Peroxide {
                t.slough.data[i] *= 0.9;
            }
        }
    }
}

/// Бактериальная нагрузка, при которой живая ткань вокруг раны иссекается вместе с ней.
const EXCISE_INFECTED: f32 = 0.25;

/// Хирургическая обработка: иссечь некроз и инфицированную ткань, соскоблить биоплёнку.
/// Иссечение за пределами раны делает её больше, но убирает очаг.
/// Дно снова кровоточит — запускается «свежее» острое заживление.
pub fn debride(t: &mut Tissue, p: &Params) {
    let excise_depth = p.skin_bottom_mm() + 0.5;
    for i in 0..t.len() {
        let infected = t.bacteria_total(i) > EXCISE_INFECTED;
        if !t.wound_mask[i] && !infected {
            continue;
        }
        if !t.wound_mask[i] {
            t.excise(i, excise_depth, p);
        }
        let removed = t.slough.data[i] + t.biofilm.data[i];
        t.slough.data[i] = 0.0;
        t.biofilm.data[i] *= 0.1;
        t.toxin.data[i] *= 0.1;
        let keep = if infected { 0.02 } else { 0.1 };
        t.bacteria.data[i] *= keep;
        t.bacteria_res.data[i] *= keep;
        if removed > 0.05 || infected {
            t.bleeding.data[i] = t.bleeding.data[i].max(0.5);
            t.debris.data[i] = t.debris.data[i].max(0.3);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tissue::WoundShape;

    const W: usize = 48;
    const H: usize = 24;

    fn wounded() -> (Tissue, Params) {
        let p = Params::default();
        let mut t = Tissue::healthy(W, H, &p);
        t.injure(WoundShape::Circle { radius: 6.0 }, 2.5, &p);
        (t, p)
    }

    #[test]
    fn enums_parse_by_name() {
        let antiseptics = ["octenidine", "polyhexanide", "hypochlorous", "chlorhexidine", "povidone", "peroxide"];
        for (k, a) in antiseptics.iter().zip(Antiseptic::ALL) {
            assert_eq!(Antiseptic::parse(k), Some(a));
        }
        for (k, a) in ["cefazolin", "clindamycin", "vancomycin"].iter().zip(Antibiotic::ALL) {
            assert_eq!(Antibiotic::parse(k), Some(a));
        }
        for (k, d) in ["hydrogel", "collagenase", "larvae"].iter().zip(Debrider::ALL) {
            assert_eq!(Debrider::parse(k), Some(d));
        }
        assert_eq!(Antiseptic::parse("vodka"), None);
        for (k, a) in Antibiotic::ALL.iter().enumerate() {
            assert_eq!(a.index(), k);
        }
    }

    #[test]
    fn antiseptic_toxicity_ranking_matches_literature() {
        let tox = |a: Antiseptic| a.props().tox_epi + a.props().tox_fib;
        let hocl = tox(Antiseptic::Hypochlorous);
        for a in Antiseptic::ALL {
            assert!(tox(a) >= hocl, "гипохлорит — самый щадящий");
            assert!(tox(a) <= tox(Antiseptic::Peroxide), "перекись — самая токсичная");
        }
        assert!(tox(Antiseptic::Octenidine) < tox(Antiseptic::Chlorhexidine));
        assert!(tox(Antiseptic::Polyhexanide) < tox(Antiseptic::PovidoneIodine));
        // Йод сильнее всех (кроме перекиси) бьёт по фибробластам, хлоргексидин — по кератиноцитам.
        assert!(Antiseptic::PovidoneIodine.props().tox_fib > Antiseptic::Chlorhexidine.props().tox_fib);
        assert!(Antiseptic::Chlorhexidine.props().tox_epi >= Antiseptic::Octenidine.props().tox_epi);
        // Индекс биосовместимости: октенидин лучше хлоргексидина.
        let bi = |a: Antiseptic| a.props().kill / tox(a);
        assert!(bi(Antiseptic::Octenidine) > bi(Antiseptic::Chlorhexidine));
    }

    #[test]
    fn antibiotic_profiles() {
        let cef = Antibiotic::Cefazolin.props();
        let cli = Antibiotic::Clindamycin.props();
        let van = Antibiotic::Vancomycin.props();
        assert!(cef.mic_resistant > cef.mic_sensitive, "устойчивый штамм устойчив к β-лактаму");
        assert_eq!(van.mic_resistant, van.mic_sensitive, "ванкомицин его убивает");
        assert!(cli.antitoxin > 0.0 && cef.antitoxin == 0.0 && van.antitoxin == 0.0);
        assert!(cli.necro_pen > cef.necro_pen && cli.tissue_pen > cef.tissue_pen);
        assert!(van.tissue_pen < cef.tissue_pen && van.emax < cef.emax);
        assert!(cli.emax < cef.emax, "клиндамицин — бактериостатик");
    }

    #[test]
    fn antibiotic_follows_half_life() {
        let (mut t, p) = wounded();
        let mut th = Therapy::default();
        th.start_antibiotic(Antibiotic::Cefazolin, 1000.0, 0.0);
        th.tick(&mut t, &p, 0.0);
        let after_dose = t.abx_plasma[0];
        let props = Antibiotic::Cefazolin.props();
        assert!((after_dose - props.dose * (-p.dt * std::f32::consts::LN_2 / props.half_life_h).exp()).abs() < 1e-4);
        let steps = (props.half_life_h / p.dt).round() as usize;
        for k in 1..=steps {
            th.tick(&mut t, &p, k as f32 * p.dt);
        }
        assert!((t.abx_plasma[0] / after_dose - 0.5).abs() < 0.02);
        assert!(t.abx_plasma[1] == 0.0 && t.abx_plasma[2] == 0.0, "другие препараты не назначены");
        assert_eq!(th.events.iter().filter(|e| e.kind == EventKind::Dose(Antibiotic::Cefazolin)).count(), 1);
    }

    #[test]
    fn repeated_doses_and_stop() {
        let (mut t, p) = wounded();
        let mut th = Therapy::default();
        th.start_antibiotic(Antibiotic::Vancomycin, 12.0, 0.0);
        let mut now = 0.0;
        for _ in 0..(36.0 / p.dt) as usize {
            th.tick(&mut t, &p, now);
            now += p.dt;
        }
        let doses = th.events.iter().filter(|e| matches!(e.kind, EventKind::Dose(_))).count();
        assert_eq!(doses, 3, "доза в 0, 12 и 24 ч");
        th.stop_antibiotic(Antibiotic::Vancomycin);
        let level = t.abx_plasma[2];
        th.tick(&mut t, &p, 48.0);
        assert!(t.abx_plasma[2] < level);
        assert_eq!(th.events.len(), 3);
    }

    #[test]
    fn antiseptic_goes_only_on_open_wound_and_washes_out() {
        let (mut t, p) = wounded();
        let mut th = Therapy::default();
        th.apply_antiseptic_now(&mut t, Antiseptic::Octenidine, 0.0);
        let c = (H / 2) * W + W / 2;
        assert_eq!(t.antiseptic.data[c], 1.0);
        assert_eq!(t.antiseptic.data[0], 0.0, "на целую кожу не наносим");
        let hl = Antiseptic::Octenidine.props().half_life_h;
        for k in 0..(hl / p.dt).round() as usize {
            th.tick(&mut t, &p, k as f32 * p.dt);
        }
        assert!((t.antiseptic.data[c] - 0.5).abs() < 0.03);
        assert_eq!(th.events.len(), 1);
    }

    #[test]
    fn debridement_removes_dead_tissue_and_biofilm() {
        let (mut t, p) = wounded();
        let c = (H / 2) * W + W / 2;
        t.slough.data[c] = 1.0;
        t.biofilm.data[c] = 0.8;
        t.bacteria.data[c] = 0.5;
        t.bleeding.data[c] = 0.0;
        debride(&mut t, &p);
        assert_eq!(t.slough.data[c], 0.0);
        assert!(t.biofilm.data[c] <= 0.08 + 1e-6);
        assert!(t.bacteria.data[c] < 0.02);
        assert!(t.bleeding.data[c] >= 0.5, "дно снова кровоточит");
    }

    #[test]
    fn debridement_excises_infected_margin() {
        let (mut t, p) = wounded();
        let i = 2; // далеко от раны
        assert!(!t.wound_mask[i]);
        t.bacteria.data[i] = 0.5;
        debride(&mut t, &p);
        assert!(t.wound_mask[i], "инфицированная ткань иссечена вместе с раной");
        assert!(t.depth.data[i] >= p.skin_bottom_mm());
        assert!(t.bacteria_total(i) < 0.05);
    }
}