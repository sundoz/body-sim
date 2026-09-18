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
    /// До какой остаточной глубины полости эпителий может наползать.
    pub epi_level_mm: f32,

    // Антибиотик (системный, условный β-лактам)
    pub abx_half_life_h: f32,
    pub abx_dose: f32,
    pub abx_emax: f32,
    pub mic_sensitive: f32,
    pub mic_resistant: f32,
    pub abx_biofilm_protect: f32,
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
            slough_clear: 0.01,
            max_depth_mm: 12.0,

            clear_neut: 0.02,
            clear_mac: 0.08,

            d_neut: 0.1,
            recruit_neut: 0.5,
            death_neut: 0.03,
            max_neut: 1.0,
            efferocytosis: 0.05,

            d_mac: 0.05,
            recruit_mac: 0.05,
            max_mac: 1.0,
            switch_m1_m2: 0.03,
            death_m1: 0.01,
            death_m2: 0.008,

            d_fib: 0.1,
            prolif_fib: 0.08,
            max_fib: 1.0,
            fib_baseline: 0.1,
            fib_return: 0.02,
            bed_fib: 0.01,
            fill_rate: 0.03,

            collagen_rate: 0.015,
            collagen_degr_m1: 0.02,
            // τ ≈ 60 дней: ~20% прочности к 3-й неделе, ~50% к 3-му месяцу.
            maturation_rate: 0.0007,
            scar_maturity_cap: 0.8,

            d_vessel: 0.05,
            angio_rate: 0.08,

            d_epi: 0.03,
            epi_rate: 0.12,
            epi_bact_half: 0.1,
            adnexal_rate: 0.002,
            epi_level_mm: 0.4,

            abx_half_life_h: 2.0,
            abx_dose: 8.0,
            abx_emax: 0.6,
            mic_sensitive: 1.0,
            mic_resistant: 4.0,
            abx_biofilm_protect: 0.9,
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
