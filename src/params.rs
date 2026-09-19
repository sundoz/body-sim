//! Параметры модели. Время — часы, длина по коже — клетки сетки (cell_mm мм),
//! глубина — мм, концентрации и плотности — безразмерные, нормированы на 0..1.
//! Значения подобраны так, чтобы воспроизводить типичную хронологию
//! заживления, а не как точные физиологические константы.

#[derive(Clone, Debug)]
pub struct Params {
    pub dt: f32,
    pub cell_mm: f32,

    // Анатомия колонки ткани (мм) и кровоснабжение
    pub epidermis_mm: f32,
    pub dermis_mm: f32,
    pub fat_mm: f32,
    /// Системная перфузия конечности: 1 — норма, <0.5 — ишемия.
    pub perfusion: f32,
    pub fat_bed_perfusion: f32,
    pub muscle_bed_perfusion: f32,

    // Гемостаз
    pub k_clot: f32,
    pub k_stop: f32,
    pub k_platelet_gf: f32,
    pub k_fibrinolysis: f32,

    // Провоспалительный сигнал (DAMPs, IL-1, TNF-α)
    pub d_signal: f32,
    pub decay_signal: f32,
    pub s_debris: f32,
    pub s_bacteria: f32,
    pub s_m1: f32,
    pub s_neut: f32,
    pub s_slough: f32,
    pub s_biofilm: f32,

    // Факторы роста (PDGF, TGF-β)
    pub d_gf: f32,
    pub decay_gf: f32,
    pub gf_m2: f32,
    pub gf_m1: f32,
    pub gf_half: f32,

    // VEGF
    pub d_vegf: f32,
    pub decay_vegf: f32,
    pub vegf_hypoxia: f32,
    pub vegf_m2: f32,
    pub vegf_half: f32,

    // Кислород
    pub d_o2: f32,
    pub o2_supply: f32,
    pub o2_bed: f32,
    pub o2_consumption: f32,
    pub hypoxia_threshold: f32,
    /// Во сколько раз воспаление может увеличить приток крови.
    pub hyperemia: f32,

    // Бактерии
    pub initial_bacteria: f32,
    pub bact_growth: f32,
    pub bact_max: f32,
    pub d_bact: f32,
    /// Способность расти в живой ткани под эпителием (0 — только на поверхности раны).
    pub bact_invasion: f32,
    pub kill_neut: f32,
    pub kill_m1: f32,
    pub kill_m2: f32,
    pub kill_saturation: f32,
    /// Фоновая защита живой ткани (комплемент, резидентные макрофаги), 1/ч.
    pub resident_kill: f32,
    /// Токсины, убивающие фагоциты (лейкоцидины): иммунитет слабеет с ростом нагрузки.
    pub leukocidin: f32,
    /// Относительная скорость роста устойчивого штамма (цена устойчивости).
    pub res_fitness: f32,
    /// Доля потомства чувствительного штамма, получающая устойчивость.
    pub mutation_rate: f32,

    // Биоплёнка
    pub biofilm_rate: f32,
    pub biofilm_decay: f32,
    pub biofilm_shed: f32,
    pub biofilm_protect: f32,

    // Некроз
    pub o2_crit: f32,
    pub k_ischemia: f32,
    pub k_toxin: f32,
    pub toxin_threshold: f32,
    /// Толщина ткани (мм), погибающая за единицу «смертности».
    pub necrosis_depth_mm: f32,
    /// Экзотоксины (стрептолизин, α-токсин): выработка на единицу бактерий, 1/ч; 0 — штамм их не выделяет.
    pub toxin_production: f32,
    pub d_toxin: f32,
    pub toxin_decay: f32,
    /// Гибель ткани на единицу концентрации токсина, 1/ч.
    pub k_toxin_field: f32,
    pub slough_clear: f32,
    pub max_depth_mm: f32,

    // Клеточный дебрис
    pub clear_neut: f32,
    pub clear_mac: f32,

    // Нейтрофилы
    pub d_neut: f32,
    pub recruit_neut: f32,
    pub death_neut: f32,
    pub max_neut: f32,
    pub efferocytosis: f32,

    // Макрофаги
    pub d_mac: f32,
    pub recruit_mac: f32,
    pub max_mac: f32,
    pub switch_m1_m2: f32,
    pub death_m1: f32,
    pub death_m2: f32,

    // Фибробласты и заполнение полости грануляциями
    pub d_fib: f32,
    pub prolif_fib: f32,
    pub max_fib: f32,
    pub fib_baseline: f32,
    pub fib_return: f32,
    pub bed_fib: f32,
    /// Скорость заполнения полости (мм/ч на единицу плотности фибробластов).
    pub fill_rate: f32,

    // Коллаген и ремоделирование
    pub collagen_rate: f32,
    pub collagen_degr_m1: f32,
    pub maturation_rate: f32,
    pub scar_maturity_cap: f32,

    // Сосуды (ангиогенез)
    pub d_vessel: f32,
    pub angio_rate: f32,

    // Эпителий (реэпителизация)
    pub d_epi: f32,
    pub epi_rate: f32,
    pub epi_bact_half: f32,
    /// Эпителий из придатков кожи (фолликулы, потовые железы) в неглубоких ранах.
    pub adnexal_rate: f32,
    /// Как круто падает число придатков (фолликулов, желёз) с глубиной раны:
    /// островков эпителия ∝ (доля уцелевшей дермы)^показатель.
    pub adnexal_exponent: f32,
    /// До какой остаточной глубины полости эпителий может наползать.
    pub epi_level_mm: f32,

    // Регенерация глубоких слоёв
    /// Скорость образования новых мышечных волокон из клеток-сателлитов, 1/ч.
    pub muscle_regen: f32,
    /// Характерная глубина потери мышцы (мм), при которой способность к регенерации падает в e раз:
    /// большие дефекты (volumetric muscle loss) не восстанавливаются, а рубцуются.
    pub muscle_loss_scale: f32,
    /// Скорость фиброза мышцы при большом дефекте, воспалении, инфекции, ишемии, 1/ч.
    pub muscle_fibrosis: f32,
    /// Скорость появления новых адипоцитов (зависит от уцелевших фолликулов), 1/ч.
    pub fat_regen: f32,
    /// Какая доля невосстановленной клетчатки остаётся вдавлением после заживления.
    pub atrophy: f32,

    // Очищение от некроза (мм мёртвой ткани в час)
    pub hydrogel_boost: f32,
    pub collagenase_rate: f32,
    pub larvae_rate: f32,
    pub larvae_kill: f32,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            dt: 0.1,
            cell_mm: 0.25,

            epidermis_mm: 0.1,
            dermis_mm: 2.0,
            fat_mm: 6.0,
            perfusion: 1.0,
            fat_bed_perfusion: 0.6,
            muscle_bed_perfusion: 0.9,

            k_clot: 2.0,
            k_stop: 3.0,
            k_platelet_gf: 0.5,
            k_fibrinolysis: 0.02,

            d_signal: 0.3,
            decay_signal: 0.2,
            s_debris: 0.2,
            s_bacteria: 0.5,
            // Слабый вклад M1 и нулевой вклад нейтрофилов: иначе петля
            // «нейтрофилы → сигнал → новые нейтрофилы» не даёт воспалению угаснуть.
            s_m1: 0.02,
            s_neut: 0.0,
            s_slough: 0.1,
            s_biofilm: 0.1,

            d_gf: 0.3,
            decay_gf: 0.1,
            gf_m2: 0.1,
            gf_m1: 0.02,
            gf_half: 0.3,

            d_vegf: 0.3,
            decay_vegf: 0.1,
            vegf_hypoxia: 0.1,
            vegf_m2: 0.03,
            vegf_half: 0.3,

            d_o2: 1.0,
            o2_supply: 1.0,
            o2_bed: 0.25,
            o2_consumption: 0.5,
            hypoxia_threshold: 0.6,
            hyperemia: 1.5,

            initial_bacteria: 0.02,
            bact_growth: 0.15,
            bact_max: 1.0,
            d_bact: 0.02,
            bact_invasion: 0.0,
            kill_neut: 1.0,
            kill_m1: 0.6,
            kill_m2: 0.1,
            kill_saturation: 0.1,
            resident_kill: 0.3,
            leukocidin: 0.0,
            res_fitness: 0.9,
            mutation_rate: 1e-6,

            biofilm_rate: 0.01,
            biofilm_decay: 0.01,
            biofilm_shed: 0.02,
            biofilm_protect: 0.8,

            o2_crit: 0.2,
            k_ischemia: 0.05,
            k_toxin: 0.03,
            toxin_threshold: 0.8,
            necrosis_depth_mm: 0.5,
            toxin_production: 0.0,
            d_toxin: 0.3,
            toxin_decay: 0.3,
            k_toxin_field: 0.15,
            slough_clear: 0.01,
            max_depth_mm: 12.0,

            clear_neut: 0.02,
            clear_mac: 0.08,

            d_neut: 0.1,
            recruit_neut: 0.5,
            death_neut: 0.064,
            max_neut: 1.0,
            efferocytosis: 0.05,

            d_mac: 0.05,
            recruit_mac: 0.0772,
            max_mac: 1.0,
            switch_m1_m2: 0.0432,
            death_m1: 0.01,
            death_m2: 0.00742,

            d_fib: 0.1,
            prolif_fib: 0.0448,
            max_fib: 1.0,
            fib_baseline: 0.1,
            fib_return: 0.0072,
            bed_fib: 0.01,
            fill_rate: 0.1,

            collagen_rate: 0.00636,
            collagen_degr_m1: 0.02,
            // τ ≈ 60 дней: ~20% прочности к 3-й неделе, ~50% к 3-му месяцу.
            maturation_rate: 0.00152,
            scar_maturity_cap: 0.8,

            d_vessel: 0.05,
            angio_rate: 0.08,

            d_epi: 0.0404,
            epi_rate: 0.0564,
            epi_bact_half: 0.1,
            adnexal_rate: 0.00002,
            adnexal_exponent: 5.9,
            epi_level_mm: 0.4,

            // Небольшая травма: 50% новых волокон к ~10-му дню, ~90% к 3–4 неделям.
            muscle_regen: 0.01,
            muscle_loss_scale: 2.5,
            muscle_fibrosis: 0.004,
            // Жир во взрослой ране почти не возвращается: несколько процентов за месяц.
            fat_regen: 0.00015,
            atrophy: 0.35,

            hydrogel_boost: 2.5,
            collagenase_rate: 0.004,
            larvae_rate: 0.012,
            larvae_kill: 0.2,
        }
    }
}

fn lerp_points(x: f32, pts: &[(f32, f32)]) -> f32 {
    if x <= pts[0].0 {
        return pts[0].1;
    }
    for w in pts.windows(2) {
        let ((x0, y0), (x1, y1)) = (w[0], w[1]);
        if x <= x1 {
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
        }
    }
    pts[pts.len() - 1].1
}

impl Params {
    pub fn skin_bottom_mm(&self) -> f32 {
        self.epidermis_mm + self.dermis_mm
    }

    pub fn fat_bottom_mm(&self) -> f32 {
        self.skin_bottom_mm() + self.fat_mm
    }

    /// Кровоснабжение дна раны в зависимости от того, до какого слоя она дошла:
    /// дерма богата сосудами, клетчатка бедна, мышца снова хорошо снабжается.
    pub fn bed_perfusion(&self, depth_mm: f32) -> f32 {
        let (s, f) = (self.skin_bottom_mm(), self.fat_bottom_mm());
        lerp_points(
            depth_mm,
            &[(s, 1.0), (s + 0.5, self.fat_bed_perfusion), (f, self.fat_bed_perfusion), (f + 0.5, self.muscle_bed_perfusion)],
        )
    }

    /// Доля уцелевшей дермы (с фолликулами и сосудами) под раной глубиной `depth_mm`.
    pub fn residual_dermis(&self, depth_mm: f32) -> f32 {
        (1.0 - (depth_mm - self.epidermis_mm) / self.dermis_mm).clamp(0.0, 1.0)
    }

    /// Сколько миллиметров жировой клетчатки утрачено в колонке с поражением до dmax.
    pub fn lost_fat_mm(&self, dmax: f32) -> f32 {
        dmax.clamp(self.skin_bottom_mm(), self.fat_bottom_mm()) - self.skin_bottom_mm()
    }

    /// Сколько миллиметров мышцы утрачено (ниже клетчатки и фасции).
    pub fn lost_muscle_mm(&self, dmax: f32) -> f32 {
        (dmax - self.fat_bottom_mm() - 0.3).max(0.0)
    }

    /// Название слоя на глубине `depth_mm`.
    pub fn layer_title(&self, depth_mm: f32) -> &'static str {
        if depth_mm <= self.epidermis_mm + 1e-3 {
            "эпидермис"
        } else if depth_mm < self.skin_bottom_mm() - 1e-3 {
            "дерма (пограничная)"
        } else if depth_mm <= self.skin_bottom_mm() + 0.5 {
            "полнослойная"
        } else if depth_mm <= self.fat_bottom_mm() {
            "жировая клетчатка"
        } else {
            "мышца"
        }
    }
}

/// Готовые клинические сценарии — модификаторы поверх базовых параметров.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Scenario {
    Healthy,
    Diabetic,
    Infected,
    Elderly,
    Ischemic,
    DiabeticFoot,
    Necrotizing,
}

impl Scenario {
    pub const ALL: [Scenario; 7] = [
        Self::Healthy,
        Self::Diabetic,
        Self::Infected,
        Self::Elderly,
        Self::Ischemic,
        Self::DiabeticFoot,
        Self::Necrotizing,
    ];

    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "healthy" => Some(Self::Healthy),
            "diabetic" => Some(Self::Diabetic),
            "infected" => Some(Self::Infected),
            "elderly" => Some(Self::Elderly),
            "ischemic" => Some(Self::Ischemic),
            "diabetic-foot" => Some(Self::DiabeticFoot),
            "necrotizing" => Some(Self::Necrotizing),
            _ => None,
        }
    }

    pub fn title(self) -> &'static str {
        match self {
            Self::Healthy => "здоровый взрослый",
            Self::Diabetic => "сахарный диабет",
            Self::Infected => "инфицированная рана",
            Self::Elderly => "пожилой возраст",
            Self::Ischemic => "ишемия конечности",
            Self::DiabeticFoot => "диабетическая стопа",
            Self::Necrotizing => "некротизирующая инфекция",
        }
    }

    pub fn short(self) -> &'static str {
        match self {
            Self::Healthy => "Здоровый",
            Self::Diabetic => "Диабет",
            Self::Infected => "Инфекция",
            Self::Elderly => "Пожилой",
            Self::Ischemic => "Ишемия",
            Self::DiabeticFoot => "Диаб. стопа",
            Self::Necrotizing => "Некротиз. инф.",
        }
    }

    pub fn params(self) -> Params {
        let mut p = Params::default();
        match self {
            Self::Healthy => {}
            Self::Diabetic => diabetic(&mut p),
            Self::Infected => {
                p.initial_bacteria = 0.5;
                p.bact_growth = 0.35;
            }
            Self::Elderly => {
                p.prolif_fib *= 0.6;
                p.epi_rate *= 0.7;
                p.angio_rate *= 0.7;
                p.recruit_mac *= 0.7;
                p.maturation_rate *= 0.7;
            }
            Self::Ischemic => {
                // Атеросклероз артерий: ткань живёт на грани, любая рана «съедает» резерв.
                p.perfusion = 0.2;
                p.angio_rate *= 0.5;
            }
            Self::DiabeticFoot => {
                diabetic(&mut p);
                p.perfusion = 0.45;
                p.initial_bacteria = 0.3;
                p.bact_growth = 0.3;
            }
            Self::Necrotizing => {
                // Высоковирулентный токсинпродуцирующий штамм, растёт в живой ткани.
                p.initial_bacteria = 0.3;
                p.bact_growth = 0.4;
                p.bact_invasion = 0.35;
                p.d_bact = 0.1;
                p.leukocidin = 4.0;
                // Токсины расходятся впереди бактерий и убивают ткань, которую те потом заселяют.
                p.toxin_production = 1.0;
                p.k_toxin = 0.12;
                p.toxin_threshold = 0.25;
            }
        }
        p
    }
}

/// Микроангиопатия, дисфункция нейтрофилов, «застревание» макрофагов в M1,
/// нарушение миграции кератиноцитов.
fn diabetic(p: &mut Params) {
    p.perfusion = 0.7;
    p.epi_rate *= 0.7;
    p.angio_rate *= 0.4;
    p.kill_neut *= 0.6;
    p.switch_m1_m2 *= 0.35;
    p.gf_m2 *= 0.6;
    p.prolif_fib *= 0.7;
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f32, b: f32) -> bool {
        (a - b).abs() < 1e-4
    }

    #[test]
    fn layer_boundaries() {
        let p = Params::default();
        assert!(close(p.skin_bottom_mm(), 2.1));
        assert!(close(p.fat_bottom_mm(), 8.1));
    }

    #[test]
    fn bed_perfusion_follows_layers() {
        let p = Params::default();
        assert!(close(p.bed_perfusion(0.0), 1.0));
        assert!(close(p.bed_perfusion(1.5), 1.0), "дерма хорошо кровоснабжается");
        assert!(close(p.bed_perfusion(5.0), p.fat_bed_perfusion), "клетчатка — плохо");
        assert!(close(p.bed_perfusion(11.0), p.muscle_bed_perfusion), "мышца — снова хорошо");
        assert!(p.bed_perfusion(5.0) < p.bed_perfusion(1.0));
        assert!(p.bed_perfusion(11.0) > p.bed_perfusion(5.0));
    }

    #[test]
    fn residual_dermis_shrinks_with_depth() {
        let p = Params::default();
        assert!(close(p.residual_dermis(0.0), 1.0));
        assert!(close(p.residual_dermis(0.1), 1.0));
        assert!(close(p.residual_dermis(1.1), 0.5));
        assert!(close(p.residual_dermis(2.1), 0.0));
        assert!(close(p.residual_dermis(9.0), 0.0));
    }

    #[test]
    fn lost_layers() {
        let p = Params::default();
        assert!(close(p.lost_fat_mm(1.0), 0.0));
        assert!(close(p.lost_fat_mm(5.0), 2.9));
        assert!(close(p.lost_fat_mm(12.0), 6.0));
        assert!(close(p.lost_muscle_mm(8.0), 0.0));
        assert!(close(p.lost_muscle_mm(9.0), 0.6));
    }

    #[test]
    fn layer_titles() {
        let p = Params::default();
        assert_eq!(p.layer_title(0.1), "эпидермис");
        assert_eq!(p.layer_title(1.0), "дерма (пограничная)");
        assert_eq!(p.layer_title(2.5), "полнослойная");
        assert_eq!(p.layer_title(5.0), "жировая клетчатка");
        assert_eq!(p.layer_title(9.0), "мышца");
    }

    #[test]
    fn scenarios_parse_and_change_the_patient() {
        let keys = ["healthy", "diabetic", "infected", "elderly", "ischemic", "diabetic-foot", "necrotizing"];
        for (key, s) in keys.iter().zip(Scenario::ALL) {
            assert_eq!(Scenario::parse(key), Some(s));
            assert!(!s.title().is_empty() && !s.short().is_empty());
        }
        assert_eq!(Scenario::parse("unknown"), None);
        let healthy = Scenario::Healthy.params();
        assert!(Scenario::Ischemic.params().perfusion < 0.5 * healthy.perfusion);
        assert!(Scenario::Diabetic.params().perfusion < healthy.perfusion);
        assert!(Scenario::Infected.params().initial_bacteria > healthy.initial_bacteria);
        let nf = Scenario::Necrotizing.params();
        assert!(nf.leukocidin > 0.0 && nf.bact_invasion > 0.0 && nf.k_toxin > healthy.k_toxin);
    }

    #[test]
    fn explicit_scheme_is_stable() {
        // Явная схема устойчива при D·dt ≤ 0.25 (5-точечный лапласиан).
        let p = Params::default();
        for d in [p.d_signal, p.d_gf, p.d_vegf, p.d_o2, p.d_bact, p.d_neut, p.d_mac, p.d_fib, p.d_vessel, p.d_epi] {
            assert!(d * p.dt <= 0.25, "D = {d}");
        }
        assert!(Scenario::Necrotizing.params().d_bact * p.dt <= 0.25);
    }
}