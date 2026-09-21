//! Место на теле: анатомия колонки кожи, местная физиология и площадь области.
//!
//! Кожа устроена по-разному в зависимости от того, где она находится, и это меняет
//! исход куда сильнее многих параметров модели: толщина дермы различается втрое,
//! подкожной клетчатки — вшестеро, плотность придатков (из которых растут островки
//! эпителия) — на порядки, а кровоснабжение голени и лица несопоставимо.
//!
//! Числа — порядковые ориентиры для взрослого человека, а не измерения конкретного
//! пациента: толщина слоёв по УЗИ- и гистологическим обзорам, доли поверхности тела —
//! по правилу девяток и таблице Ланда–Браудера. Площади считаются от опорной
//! поверхности тела 1.73 м².
//!
//! Предплечье — опорная точка: именно на этой анатомии откалибрована модель, поэтому
//! его значения совпадают с прежними значениями по умолчанию, а не измерены отдельно.

use crate::params::Params;

/// Площадь поверхности тела взрослого человека, см² (1.73 м²).
pub const BODY_SURFACE_CM2: f32 = 17_300.0;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum BodySite {
    /// Опорная анатомия, на которой калибровалась модель.
    #[default]
    Forearm,
    Face,
    Scalp,
    Back,
    Abdomen,
    /// Передняя поверхность голени — классическое место незаживающих язв.
    Shin,
    /// Подошва — безволосая кожа с очень толстым роговым слоем.
    Sole,
    /// Крестец — место пролежней.
    Sacrum,
}

/// Чем место на теле отличается от опорного предплечья.
#[derive(Clone, Copy, Debug)]
pub struct SiteProps {
    pub epidermis_mm: f32,
    pub dermis_mm: f32,
    pub fat_mm: f32,
    /// Множитель местного кровоснабжения поверх системной перфузии.
    pub perfusion: f32,
    /// Множитель плотности придатков кожи (фолликулы, железы) — источников островков эпителия.
    pub adnexal: f32,
    /// Фоновая обсеменённость кожи: нижняя граница начальной бактериальной нагрузки.
    pub flora: f32,
    /// Какую долю поверхности тела занимает эта область, % (правило девяток).
    pub tbsa_percent: f32,
}

impl BodySite {
    pub const ALL: [BodySite; 8] =
        [Self::Forearm, Self::Face, Self::Scalp, Self::Back, Self::Abdomen, Self::Shin, Self::Sole, Self::Sacrum];

    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "forearm" => Some(Self::Forearm),
            "face" => Some(Self::Face),
            "scalp" => Some(Self::Scalp),
            "back" => Some(Self::Back),
            "abdomen" => Some(Self::Abdomen),
            "shin" => Some(Self::Shin),
            "sole" | "foot" => Some(Self::Sole),
            "sacrum" => Some(Self::Sacrum),
            _ => None,
        }
    }

    pub fn key(self) -> &'static str {
        match self {
            Self::Forearm => "forearm",
            Self::Face => "face",
            Self::Scalp => "scalp",
            Self::Back => "back",
            Self::Abdomen => "abdomen",
            Self::Shin => "shin",
            Self::Sole => "sole",
            Self::Sacrum => "sacrum",
        }
    }

    pub fn title(self) -> &'static str {
        match self {
            Self::Forearm => "предплечье",
            Self::Face => "лицо",
            Self::Scalp => "волосистая часть головы",
            Self::Back => "спина",
            Self::Abdomen => "живот",
            Self::Shin => "передняя поверхность голени",
            Self::Sole => "подошва",
            Self::Sacrum => "крестец",
        }
    }

    pub fn short(self) -> &'static str {
        match self {
            Self::Forearm => "Предплечье",
            Self::Face => "Лицо",
            Self::Scalp => "Голова",
            Self::Back => "Спина",
            Self::Abdomen => "Живот",
            Self::Shin => "Голень",
            Self::Sole => "Подошва",
            Self::Sacrum => "Крестец",
        }
    }

    /// Чем это место примечательно клинически — подсказка в интерфейсе.
    pub fn note(self) -> &'static str {
        match self {
            Self::Forearm => "опорная анатомия: на ней откалибрована модель",
            Self::Face => "тонкая дерма, обильное кровоснабжение — заживает быстрее всего",
            Self::Scalp => "рекордная плотность фолликулов: эпителий растёт островками",
            Self::Back => "самая толстая дерма: глубокий дефект и грубый рубец",
            Self::Abdomen => "толстая клетчатка: до мышцы очень далеко",
            Self::Shin => "тонкая клетчатка и плохой кровоток — язвы не заживают",
            Self::Sole => "толстый роговой слой, нет волос, нагрузка весом",
            Self::Sacrum => "пролежень: давление, ишемия и богатая флора",
        }
    }

    pub fn props(self) -> SiteProps {
        let p = |epidermis_mm, dermis_mm, fat_mm, perfusion, adnexal, flora, tbsa_percent| SiteProps {
            epidermis_mm,
            dermis_mm,
            fat_mm,
            perfusion,
            adnexal,
            flora,
            tbsa_percent,
        };
        match self {
            // Значения по умолчанию модели: менять их — значит сдвинуть калибровку.
            Self::Forearm => p(0.1, 2.0, 6.0, 1.0, 1.0, 0.02, 3.0),
            // Лицо: тонкая кожа, густая сеть сосудов, много сальных желёз и пушковых волос.
            Self::Face => p(0.07, 1.2, 2.0, 1.35, 3.0, 0.03, 3.0),
            // Голова: фолликулов больше, чем где-либо; донорские места здесь заживают заметно быстрее.
            Self::Scalp => p(0.08, 1.8, 3.0, 1.3, 8.0, 0.03, 3.0),
            // Спина: дерма до 3–4 мм, отсюда склонность к гипертрофическим рубцам.
            Self::Back => p(0.1, 3.5, 8.0, 0.9, 0.8, 0.02, 13.0),
            Self::Abdomen => p(0.08, 2.0, 12.0, 0.9, 0.7, 0.02, 9.0),
            // Голень: клетчатки почти нет, кровоток зависимый — отсюда венозные язвы.
            Self::Shin => p(0.09, 1.5, 2.0, 0.55, 0.6, 0.03, 7.0),
            // Подошва: роговой слой до миллиметра, волос нет — островки дают только потовые железы.
            Self::Sole => p(1.0, 2.0, 10.0, 0.8, 0.35, 0.05, 1.5),
            // Крестец: давление на костный выступ, рядом промежностная флора.
            Self::Sacrum => p(0.1, 2.0, 4.0, 0.6, 0.6, 0.08, 1.2),
        }
    }

    /// Площадь всей этой области тела, см².
    pub fn region_cm2(self) -> f32 {
        BODY_SURFACE_CM2 * self.props().tbsa_percent / 100.0
    }

    /// Применить анатомию и местную физиологию к параметрам.
    /// Вызывается до модификаторов сценария: место задаёт ткань, сценарий — пациента.
    pub fn apply(self, p: &mut Params) {
        let s = self.props();
        p.epidermis_mm = s.epidermis_mm;
        p.dermis_mm = s.dermis_mm;
        p.fat_mm = s.fat_mm;
        p.site_perfusion = s.perfusion;
        p.adnexal_rate *= s.adnexal;
        p.initial_bacteria = p.initial_bacteria.max(s.flora);
        // Колонка ткани должна вмещать кожу, клетчатку и слой мышцы под фасцией.
        p.max_depth_mm = p.fat_bottom_mm() + p.muscle_mm;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sites_parse_and_are_described() {
        let keys = ["forearm", "face", "scalp", "back", "abdomen", "shin", "sole", "sacrum"];
        for (key, site) in keys.iter().zip(BodySite::ALL) {
            assert_eq!(BodySite::parse(key), Some(site));
            assert_eq!(site.key(), *key);
            assert!(!site.title().is_empty() && !site.short().is_empty() && !site.note().is_empty());
        }
        assert_eq!(BodySite::parse("tail"), None);
        assert_eq!(BodySite::parse("foot"), Some(BodySite::Sole), "стопа — синоним подошвы");
    }

    #[test]
    fn forearm_is_the_calibration_reference() {
        // Модель откалибрована на этой анатомии: предплечье обязано не менять ничего,
        // иначе tests/literature.rs и клиническая регрессия перестанут что-либо значить.
        let base = Params::default();
        let mut p = Params::default();
        BodySite::Forearm.apply(&mut p);
        assert_eq!(p.epidermis_mm, base.epidermis_mm);
        assert_eq!(p.dermis_mm, base.dermis_mm);
        assert_eq!(p.fat_mm, base.fat_mm);
        assert_eq!(p.site_perfusion, base.site_perfusion);
        assert_eq!(p.adnexal_rate, base.adnexal_rate);
        assert_eq!(p.initial_bacteria, base.initial_bacteria);
        assert_eq!(p.max_depth_mm, base.max_depth_mm);
    }

    #[test]
    fn anatomy_differs_between_sites() {
        let thickness = |s: BodySite| {
            let mut p = Params::default();
            s.apply(&mut p);
            p
        };
        let back = thickness(BodySite::Back);
        let face = thickness(BodySite::Face);
        let shin = thickness(BodySite::Shin);
        let sole = thickness(BodySite::Sole);
        let abdomen = thickness(BodySite::Abdomen);
        assert!(back.dermis_mm > face.dermis_mm * 2.0, "дерма спины втрое толще кожи лица");
        assert!(sole.epidermis_mm > back.epidermis_mm * 5.0, "роговой слой подошвы очень толстый");
        assert!(shin.fat_mm < abdomen.fat_mm, "на голени клетчатки почти нет");
        assert!(face.site_perfusion > shin.site_perfusion, "лицо кровоснабжается лучше голени");
        // Колонка всегда вмещает кожу с клетчаткой плюс мышцу.
        for s in BodySite::ALL {
            let p = thickness(s);
            assert!(p.max_depth_mm > p.fat_bottom_mm(), "{}: под клетчаткой должна быть мышца", s.key());
        }
    }

    #[test]
    fn appendage_density_spans_scalp_to_sole() {
        let rate = |s: BodySite| {
            let mut p = Params::default();
            s.apply(&mut p);
            p.adnexal_rate
        };
        assert!(rate(BodySite::Scalp) > rate(BodySite::Forearm));
        assert!(rate(BodySite::Sole) < rate(BodySite::Forearm), "на подошве нет волосяных фолликулов");
        assert!(rate(BodySite::Scalp) > rate(BodySite::Sole) * 20.0);
    }

    #[test]
    fn body_regions_add_up_to_a_plausible_body() {
        for s in BodySite::ALL {
            let pct = s.props().tbsa_percent;
            assert!((0.5..=20.0).contains(&pct), "{}: доля поверхности {pct}%", s.key());
            let cm2 = s.region_cm2();
            assert!((cm2 - BODY_SURFACE_CM2 * pct / 100.0).abs() < 1e-3);
        }
        // Правило девяток: спина — самая большая из перечисленных областей, подошва — одна из самых малых.
        assert!(BodySite::Back.region_cm2() > BodySite::Forearm.region_cm2());
        assert!(BodySite::Sole.region_cm2() < BodySite::Shin.region_cm2());
        // Перечисленные области не покрывают всё тело, но и не могут его превысить.
        let total: f32 = BodySite::ALL.iter().map(|s| s.props().tbsa_percent).sum();
        assert!(total < 100.0, "сумма долей {total}% не может превышать поверхность тела");
    }

    #[test]
    fn site_sets_a_floor_for_skin_flora_without_lowering_an_infected_wound() {
        let mut p = Params::default();
        BodySite::Sacrum.apply(&mut p);
        assert!(p.initial_bacteria > Params::default().initial_bacteria, "у крестца богатая флора");
        let mut dirty = Params { initial_bacteria: 0.5, ..Params::default() };
        BodySite::Face.apply(&mut dirty);
        assert_eq!(dirty.initial_bacteria, 0.5, "место не может очистить уже заражённую рану");
    }
}
