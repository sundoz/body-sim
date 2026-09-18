//! Один шаг модели: реакция–диффузия, явная схема Эйлера.
//! Все локальные члены читают значения клетки до её обновления,
//! лапласианы считаются заранее по старому состоянию.

use crate::params::Params;
use crate::tissue::Tissue;

/// Буферы лапласианов, переиспользуемые между шагами.
pub struct Scratch {
    signal: Vec<f32>,
    gf: Vec<f32>,
    vegf: Vec<f32>,
    oxygen: Vec<f32>,
    bacteria: Vec<f32>,
    bacteria_res: Vec<f32>,
    neutrophils: Vec<f32>,
    m1: Vec<f32>,
    m2: Vec<f32>,
    fibroblasts: Vec<f32>,
    vessels: Vec<f32>,
    epithelium: Vec<f32>,
}

impl Scratch {
    pub fn new(n: usize) -> Self {
        let z = || vec![0.0; n];
        Self {
            signal: z(),
            gf: z(),
            vegf: z(),
            oxygen: z(),
            bacteria: z(),
            bacteria_res: z(),
            neutrophils: z(),
            m1: z(),
            m2: z(),
            fibroblasts: z(),
            vessels: z(),
            epithelium: z(),
        }
    }
}

#[inline]
fn sat(x: f32, half: f32) -> f32 {
    x / (x + half)
}

/// Эффект антибиотика по Хиллу: доля максимального убийства при концентрации `c`.
#[inline]
fn hill(c: f32, mic: f32) -> f32 {
    let (c2, m2) = (c * c, mic * mic);
    c2 / (c2 + m2)
}

pub fn step(t: &mut Tissue, p: &Params, s: &mut Scratch) {
    t.signal.laplacian_into(&mut s.signal);
    t.growth_factor.laplacian_into(&mut s.gf);
    t.vegf.laplacian_into(&mut s.vegf);
    t.oxygen.laplacian_into(&mut s.oxygen);
    t.bacteria.laplacian_into(&mut s.bacteria);
    t.bacteria_res.laplacian_into(&mut s.bacteria_res);
    t.neutrophils.laplacian_into(&mut s.neutrophils);
    t.m1.laplacian_into(&mut s.m1);
    t.m2.laplacian_into(&mut s.m2);
    t.fibroblasts.laplacian_into(&mut s.fibroblasts);
    t.vessels.laplacian_into(&mut s.vessels);
    t.epithelium.laplacian_into(&mut s.epithelium);

    let dt = p.dt;
    let asp = t.antiseptic_agent.props();
    for i in 0..t.len() {
        let bl = t.bleeding.data[i];
        let cf = t.clot.data[i];
        let debris = t.debris.data[i];
        let b = t.bacteria.data[i];
        let br = t.bacteria_res.data[i];
        let bf = t.biofilm.data[i];
        let depth = t.depth.data[i];
        let dmax = t.depth_max.data[i];
        let slough = t.slough.data[i];
        let n = t.neutrophils.data[i];
        let m1 = t.m1.data[i];
        let m2 = t.m2.data[i];
        let f = t.fibroblasts.data[i];
        let c = t.collagen.data[i];
        let q = t.maturity.data[i];
        let v = t.vessels.data[i];
        let e = t.epithelium.data[i];
        let o = t.oxygen.data[i];
        let sig = t.signal.data[i];
        let g = t.growth_factor.data[i];
        let a = t.vegf.data[i];
        let asc = t.antiseptic.data[i];

        let bt = b + br;
        let wounded = dmax > 0.05;
        let bq = p.bed_perfusion(dmax);
        // Мёртвая ткань отрезает дно раны: нет доступа ни клеткам, ни крови, ни лекарствам.
        let block = sat(slough, 0.3);

        // --- Гемостаз: кровотечение → фибриновый сгусток → тромбоцитарные факторы роста.
        let clot_form = p.k_clot * bl * (1.0 - cf);
        // Сгусток растворяется плазмином и замещается коллагеном.
        let d_clot = clot_form - p.k_fibrinolysis * cf * c - 0.005 * cf;
        let d_bleed = -p.k_stop * cf * bl;

        // --- Кислород: подаётся сосудами и дном раны, тратится клетками и бактериями.
        let load = 0.1 + n + m1 + m2 + f + 0.5 * bt;
        // Воспаление расширяет сосуды (гиперемия) — но только если артерии на это способны.
        let hyperemia = 1.0 + p.hyperemia * sig.min(1.0) * p.perfusion;
        let supply = p.perfusion * hyperemia * (p.o2_supply * v + p.o2_bed * bq * (1.0 - block));
        let d_o2 = p.d_o2 * s.oxygen[i] + supply * (1.0 - o) - p.o2_consumption * o * load;
        let hypoxia = ((p.hypoxia_threshold - o) / p.hypoxia_threshold).max(0.0);
        let oxy = sat(o, 0.2);

        // --- Лекарства в ткани. Антибиотик приходит с кровью: в некроз и в дно
        // с плохой перфузией он почти не попадает.
        let bed_blood = if wounded { 0.5 * bq } else { 0.0 };
        let perf_local = (v + bed_blood).min(1.0) * p.perfusion * (1.0 - block);
        let abx = t.abx_plasma * perf_local;
        let abx_shield = 1.0 - p.abx_biofilm_protect * bf;
        let abx_kill_s = p.abx_emax * hill(abx, p.mic_sensitive) * abx_shield;
        let abx_kill_r = p.abx_emax * hill(abx, p.mic_resistant) * abx_shield;
        let as_kill = asp.kill * asc * (1.0 - (1.0 - asp.biofilm_pen) * bf);

        // --- Бактерии: растут на открытой поверхности и в мёртвой ткани;
        // инвазивные штаммы — и в живой ткани.
        let habitat = (1.0 - e).max(block).max(p.bact_invasion);
        let crowd = (1.0 - bt / p.bact_max).max(0.0);
        let grow_s = p.bact_growth * b * crowd * habitat;
        let grow_r = p.bact_growth * p.res_fitness * br * crowd * habitat;
        let mutate = p.mutation_rate * grow_s;
        // Фагоцитоз насыщается при высокой нагрузке; биоплёнка прячет бактерии.
        // Живое дно раны тоже защищено резидентными клетками — если не закрыто некрозом.
        let resident = p.resident_kill * perf_local;
        let immune = (p.kill_neut * n + p.kill_m1 * m1 + p.kill_m2 * m2 + resident) / (1.0 + bt / p.kill_saturation)
            * (1.0 - p.biofilm_protect * bf)
            / (1.0 + p.leukocidin * bt);
        // Биоплёнка постоянно «выпускает» планктонные бактерии.
        let shed = p.biofilm_shed * bf * crowd;
        let frac_r = if bt > 1e-6 { br / bt } else { 0.0 };
        let d_bact = p.d_bact * s.bacteria[i] + grow_s - mutate - (immune + as_kill + abx_kill_s) * b
            + shed * (1.0 - frac_r);
        let d_bact_res = p.d_bact * s.bacteria_res[i] + grow_r + mutate - (immune + as_kill + abx_kill_r) * br
            + shed * frac_r;

        // --- Биоплёнка: матрикс на открытой поверхности, где бактерии держатся долго.
        let surface = (1.0 - e).max(block);
        let d_biofilm = p.biofilm_rate * bt * (1.0 - bf) * surface
            - p.biofilm_decay * bf * (1.0 - bt / 0.05).max(0.0)
            - 0.05 * bf * e
            - asp.disrupt * asc * bf;

        // --- Некроз: гипоксия и бактериальные токсины убивают ткань колонки сверху вниз.
        // Слой мёртвой ткани сам не потребляет кислород и отделяет живое дно от бактерий,
        // поэтому под толстым струпом некроз замедляется.
        let death = p.k_ischemia * ((p.o2_crit - o) / p.o2_crit).max(0.0) * (1.0 - block)
            + p.k_toxin * (bt - p.toxin_threshold).max(0.0) * (1.0 - 0.5 * block);
        let room = (p.max_depth_mm - depth).max(0.0);
        let dead_mm = (death * p.necrosis_depth_mm).min(room / dt);
        // Аутолиз: макрофаги постепенно растворяют мёртвую ткань.
        let autolysis = p.slough_clear * (m1 + m2) * sat(slough, 0.5);

        // --- Дебрис: погибшие клетки и обломки матрикса, убираются фагоцитами.
        let d_debris = -(p.clear_neut * n + p.clear_mac * (m1 + m2)) * debris + death;

        // --- Сигналы.
        let d_sig = p.d_signal * s.signal[i]
            + p.s_debris * debris
            + p.s_bacteria * bt
            + p.s_m1 * m1
            + p.s_neut * n
            + p.s_slough * block
            + p.s_biofilm * bf
            - p.decay_signal * sig;
        let gfh = sat(g, p.gf_half);
        let d_gf = p.d_gf * s.gf[i] + p.k_platelet_gf * clot_form + p.gf_m2 * m2 + p.gf_m1 * m1
            - p.decay_gf * g;
        let d_vegf = p.d_vegf * s.vegf[i]
            + p.vegf_hypoxia * hypoxia * (0.1 + m1 + m2 + f)
            + p.vegf_m2 * m2
            - p.decay_vegf * a;

        // --- Воспаление. Клетки выходят из сосудов края и дна раны.
        let access = (0.3 * bq + v) * p.perfusion * (1.0 - 0.7 * block);
        let recruit_n = p.recruit_neut * sig * access * (1.0 - n / p.max_neut).max(0.0);
        let effero = p.efferocytosis * (m1 + m2) * n;
        let d_neut = p.d_neut * s.neutrophils[i] + recruit_n - p.death_neut * n - effero;

        let recruit_m = p.recruit_mac * sig * access * (1.0 - (m1 + m2) / p.max_mac).max(0.0);
        // Переключение M1→M2 запускается поеданием апоптотических нейтрофилов
        // и подавляется бактериальной нагрузкой.
        let switch = p.switch_m1_m2 * (0.3 + sat(effero, 0.005)) / (1.0 + 10.0 * bt) * m1;
        let d_m1 = p.d_mac * s.m1[i] + recruit_m - switch - p.death_m1 * m1;
        let d_m2 = p.d_mac * s.m2[i] + switch - p.death_m2 * m2;

        // --- Пролиферация. Клетки мигрируют только по матриксу (сгусток или коллаген),
        // в глубоких ранах грануляции растут и со дна.
        let scaffold = (cf + c).clamp(0.05, 1.0) * (1.0 - block);
        let bed_source = if depth > 0.05 { p.bed_fib * bq * gfh * (1.0 - f).max(0.0) * (1.0 - block) } else { 0.0 };
        let d_fib = p.d_fib * s.fibroblasts[i] * scaffold
            + p.prolif_fib * f * gfh * oxy * (1.0 - f / p.max_fib)
            + bed_source
            - p.fib_return * (f - p.fib_baseline) * c
            - asp.cytotox * asc * f * 0.5;
        // Грануляции заполняют полость снизу; через мёртвую ткань не прорастают.
        let fill = if depth > slough { p.fill_rate * f * oxy * (1.0 - block) } else { 0.0 };

        let deposit = p.collagen_rate * f * oxy * (1.0 - c).max(0.0) * (0.5 + 0.5 * gfh);
        // MMP макрофагов перестраивают только повреждённую дерму.
        let degrade = p.collagen_degr_m1 * m1 * c * (1.0 - p.residual_dermis(dmax)) + 0.5 * death * c;

        let d_vessel = (p.d_vessel * s.vessels[i] + p.angio_rate * sat(a, p.vegf_half) * v * (1.0 - v)) * scaffold;

        // Эпителий наползает только на выровненное грануляциями дно;
        // в неглубоких ранах он прорастает ещё и островками из придатков кожи.
        let level = 1.0 / (1.0 + ((depth - p.epi_level_mm).max(0.0) / 0.4).powi(2));
        let epi_scaffold = (cf + c).min(1.0) * level * (1.0 - block);
        let islands = if wounded { p.residual_dermis(dmax) } else { 0.0 };
        let d_epi = (p.d_epi * s.epithelium[i]
            + p.epi_rate * e * (1.0 - e) * (0.5 + 0.5 * gfh) * oxy / (1.0 + bt / p.epi_bact_half))
            * epi_scaffold
            + p.adnexal_rate * islands * (1.0 - e) * level * (1.0 - block) * oxy / (1.0 + bt / p.epi_bact_half)
            - asp.cytotox * asc * e;

        // --- Запись.
        t.bleeding.data[i] = (bl + dt * d_bleed).clamp(0.0, 1.0);
        t.clot.data[i] = (cf + dt * d_clot).clamp(0.0, 1.0);
        t.debris.data[i] = (debris + dt * d_debris).max(0.0);
        let b_new = (b + dt * d_bact).max(0.0);
        t.bacteria.data[i] = if b_new < 1e-5 { 0.0 } else { b_new };
        let br_new = (br + dt * d_bact_res).max(0.0);
        t.bacteria_res.data[i] = if br_new < 1e-9 { 0.0 } else { br_new };
        t.biofilm.data[i] = (bf + dt * d_biofilm).clamp(0.0, 1.0);
        t.oxygen.data[i] = (o + dt * d_o2).clamp(0.0, 1.0);
        t.signal.data[i] = (sig + dt * d_sig).max(0.0);
        t.growth_factor.data[i] = (g + dt * d_gf).max(0.0);
        t.vegf.data[i] = (a + dt * d_vegf).max(0.0);
        t.neutrophils.data[i] = (n + dt * d_neut).max(0.0);
        t.m1.data[i] = (m1 + dt * d_m1).max(0.0);
        t.m2.data[i] = (m2 + dt * d_m2).max(0.0);
        t.fibroblasts.data[i] = ((f + dt * d_fib) * (1.0 - dt * death)).max(0.0);
        t.vessels.data[i] = ((v + dt * d_vessel) * (1.0 - dt * death)).clamp(0.0, 1.0);
        t.epithelium.data[i] = ((e + dt * d_epi) * (1.0 - dt * death)).clamp(0.0, 1.0);
        t.antibiotic.data[i] = abx;

        let slough_new = (slough + dt * (dead_mm - autolysis)).max(0.0);
        let depth_new = (depth + dt * (dead_mm - fill)).clamp(0.0, p.max_depth_mm).max(slough_new);
        t.slough.data[i] = slough_new;
        t.depth.data[i] = depth_new;
        t.depth_max.data[i] = dmax.max(depth_new);
        if slough_new > 0.3 || depth_new > 0.3 {
            t.wound_mask[i] = true;
        }

        // --- Коллаген и ремоделирование: новый коллаген незрелый (тип III),
        // со временем созревает. В полнослойном рубце — не выше scar_maturity_cap,
        // при уцелевшей дерме заживает почти без рубца.
        let c_new = (c + dt * (deposit - degrade)).clamp(0.0, 1.0);
        let mature = (q * c * (1.0 - dt * (degrade / c.max(1e-6)))).max(0.0);
        let mut q_new = if c_new > 1e-6 { (mature / c_new).min(1.0) } else { 0.0 };
        let cap = p.scar_maturity_cap + (1.0 - p.scar_maturity_cap) * p.residual_dermis(dmax);
        q_new += dt * p.maturation_rate * (cap - q_new).max(0.0) / (1.0 + 5.0 * m1);
        t.collagen.data[i] = c_new;
        t.maturity.data[i] = q_new;
    }
}
