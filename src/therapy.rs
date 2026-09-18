//! Лечение: местные антисептики, системный антибиотик, хирургическая обработка.

use crate::params::Params;
use crate::tissue::Tissue;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Antiseptic {
    PovidoneIodine,
    Chlorhexidine,
    Octenidine,
    Peroxide,
}

/// Свойства антисептика при концентрации 1.0 (сразу после нанесения).
#[derive(Clone, Copy, Debug)]
pub struct AntisepticProps {
    /// Скорость гибели бактерий, 1/ч.
    pub kill: f32,
    /// Какая доля действия проникает в биоплёнку.
    pub biofilm_pen: f32,
    /// Разрушение матрикса биоплёнки, 1/ч.
    pub disrupt: f32,
    /// Токсичность для кератиноцитов и фибробластов, 1/ч.
    pub cytotox: f32,
    pub half_life_h: f32,
}

impl Antiseptic {
    pub const ALL: [Antiseptic; 4] = [Self::Octenidine, Self::Chlorhexidine, Self::PovidoneIodine, Self::Peroxide];

    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "povidone" | "iodine" => Some(Self::PovidoneIodine),
            "chlorhexidine" => Some(Self::Chlorhexidine),
            "octenidine" => Some(Self::Octenidine),
            "peroxide" => Some(Self::Peroxide),
            _ => None,
        }
    }

    pub fn title(self) -> &'static str {
        match self {
            Self::PovidoneIodine => "Повидон-йод",
            Self::Chlorhexidine => "Хлоргексидин",
            Self::Octenidine => "Октенидин",
            Self::Peroxide => "Перекись H₂O₂",
        }
    }

    pub fn props(self) -> AntisepticProps {
        match self {
            Self::PovidoneIodine => {
                AntisepticProps { kill: 3.0, biofilm_pen: 0.5, disrupt: 0.05, cytotox: 0.04, half_life_h: 3.0 }
            }
            Self::Chlorhexidine => {
                AntisepticProps { kill: 2.0, biofilm_pen: 0.3, disrupt: 0.02, cytotox: 0.04, half_life_h: 6.0 }
            }
            Self::Octenidine => {
                AntisepticProps { kill: 3.0, biofilm_pen: 0.6, disrupt: 0.08, cytotox: 0.015, half_life_h: 6.0 }
            }
            // Перекись в основном «вспенивает» и механически чистит, бактерицидность слабая.
            Self::Peroxide => {
                AntisepticProps { kill: 1.0, biofilm_pen: 0.2, disrupt: 0.1, cytotox: 0.12, half_life_h: 0.3 }
            }
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum EventKind {
    Antiseptic,
    Dose,
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
    pub antibiotic_every_h: Option<f32>,
    pub events: Vec<Event>,
    next_antiseptic_h: f32,
    next_dose_h: f32,
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
    pub fn start_antibiotic(&mut self, every_h: f32, now: f32) {
        self.antibiotic_every_h = Some(every_h);
        self.next_dose_h = now;
    }

    pub fn stop_antibiotic(&mut self) {
        self.antibiotic_every_h = None;
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
        if let Some(every) = self.antibiotic_every_h {
            if now >= self.next_dose_h {
                t.abx_plasma += p.abx_dose;
                self.events.push(Event { hours: now, kind: EventKind::Dose });
                self.next_dose_h = now + every;
            }
        }

        let ln2 = std::f32::consts::LN_2;
        t.abx_plasma *= (-p.dt * ln2 / p.abx_half_life_h).exp();
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
        let keep = if infected { 0.02 } else { 0.1 };
        t.bacteria.data[i] *= keep;
        t.bacteria_res.data[i] *= keep;
        if removed > 0.05 || infected {
            t.bleeding.data[i] = t.bleeding.data[i].max(0.5);
            t.debris.data[i] = t.debris.data[i].max(0.3);
        }
    }
}
