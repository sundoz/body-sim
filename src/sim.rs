//! Один шаг модели: реакция–диффузия, явная схема Эйлера.
//! Все локальные члены читают значения клетки до её обновления,
//! лапласианы считаются заранее по старому состоянию.

use rayon::prelude::*;

use crate::grid::Field;
use crate::params::Params;
use crate::therapy::{Antibiotic, AntibioticProps, AntisepticProps, Debrider};
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
    toxin: Vec<f32>,
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
            toxin: z(),
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

/// Сколько клеток сетки обрабатывает один поток за раз (4 строки по 96).
const CHUNK: usize = 384;

/// Изменяемые срезы всех полей ткани для одного блока клеток: блоки не пересекаются,
/// поэтому их можно считать параллельно без блокировок.
macro_rules! chunked_fields {
    ($($name:ident),* $(,)?) => {
        struct Chunk<'a> {
            base: usize,
            len: usize,
            wound_mask: &'a mut [bool],
            $($name: &'a mut [f32],)*
        }

        fn split_chunks(t: &mut Tissue, size: usize) -> Vec<Chunk<'_>> {
            let mut wound_mask = t.wound_mask.chunks_mut(size);
            $(let mut $name = t.$name.data.chunks_mut(size);)*
            let mut out = Vec::new();
            let mut base = 0;
            while let (Some(wound_mask), $(Some($name),)*) = (wound_mask.next(), $($name.next(),)*) {
                let len = wound_mask.len();
                out.push(Chunk { base, len, wound_mask, $($name,)* });
                base += len;
            }
            out
        }
    };
}

chunked_fields!(
    antibiotic, antiseptic, bacteria, bacteria_res, biofilm, bleeding, clot, collagen, debris, depth,
    depth_max, epithelium, fibroblasts, growth_factor, m1, m2, maturity, neutrophils, oxygen, signal,
    slough, vegf, vessels, myo, muscle_scar, fat_new, toxin,
);

/// Всё, что шаг берёт из состояния ткани целиком (лечение).
struct Globals {
    asp: AntisepticProps,
    abx: [AntibioticProps; 3],
    abx_plasma: [f32; 3],
    hydrogel: bool,
    collagenase: bool,
    larvae: bool,
}

pub fn step(t: &mut Tissue, p: &Params, s: &mut Scratch) {
    let lap: [(&Field, &mut Vec<f32>); 13] = [
        (&t.signal, &mut s.signal),
        (&t.growth_factor, &mut s.gf),
        (&t.vegf, &mut s.vegf),
        (&t.oxygen, &mut s.oxygen),
        (&t.bacteria, &mut s.bacteria),
        (&t.bacteria_res, &mut s.bacteria_res),
        (&t.neutrophils, &mut s.neutrophils),
        (&t.m1, &mut s.m1),
        (&t.m2, &mut s.m2),
        (&t.fibroblasts, &mut s.fibroblasts),
        (&t.vessels, &mut s.vessels),
        (&t.epithelium, &mut s.epithelium),
        (&t.toxin, &mut s.toxin),
    ];
    lap.into_par_iter().for_each(|(field, out)| field.laplacian_into(out));

    let gl = Globals {
        asp: t.antiseptic_agent.props(),
        abx: Antibiotic::ALL.map(|d| d.props()),
        abx_plasma: t.abx_plasma,
        hydrogel: t.debriders[Debrider::Hydrogel.index()],
        collagenase: t.debriders[Debrider::Collagenase.index()],
        larvae: t.debriders[Debrider::Larvae.index()],
    };
    let s = &*s;
    split_chunks(t, CHUNK).into_par_iter().for_each(|mut ch| step_chunk(&mut ch, s, p, &gl));
}

/// Доля восстановленной ткани при углублении дефекта: новое «пустое» место ещё не восстановлено.
fn rescale(frac: f32, lost_old: f32, lost_new: f32) -> f32 {
    if lost_new > lost_old + 1e-6 {
        frac * lost_old / lost_new
    } else {
        frac
    }
}

fn step_chunk(ch: &mut Chunk, s: &Scratch, p: &Params, gl: &Globals) {
    let dt = p.dt;
    let asp = &gl.asp;
    for j in 0..ch.len {
        let i = ch.base + j;
        let bl = ch.bleeding[j];
        let cf = ch.clot[j];
        let debris = ch.debris[j];
        let b = ch.bacteria[j];
        let br = ch.bacteria_res[j];
        let bf = ch.biofilm[j];
        let depth = ch.depth[j];
        let dmax = ch.depth_max[j];
        let slough = ch.slough[j];
        let n = ch.neutrophils[j];
        let m1 = ch.m1[j];
        let m2 = ch.m2[j];
        let f = ch.fibroblasts[j];
        let c = ch.collagen[j];
        let q = ch.maturity[j];
        let v = ch.vessels[j];
        let e = ch.epithelium[j];
        let o = ch.oxygen[j];
        let sig = ch.signal[j];
        let g = ch.growth_factor[j];
        let a = ch.vegf[j];
        let asc = ch.antiseptic[j];
        let tox = ch.toxin[j];
        let myo = ch.myo[j];
        let mscar = ch.muscle_scar[j];
        let fat_new = ch.fat_new[j];

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
        let blood = (v + bed_blood).min(1.0) * p.perfusion;
        let perf_local = blood * (1.0 - block);
        let (mut abx_kill_s, mut abx_kill_r, mut antitoxin, mut abx_level) = (0.0f32, 0.0f32, 0.0f32, 0.0f32);
        for (k, pr) in gl.abx.iter().enumerate() {
            let conc = gl.abx_plasma[k] * blood * pr.tissue_pen * (1.0 - block * (1.0 - pr.necro_pen));
            if conc < 1e-4 {
                continue;
            }
            let shield = 1.0 - pr.biofilm_protect * bf;
            abx_kill_s += pr.emax * hill(conc, pr.mic_sensitive) * shield;
            abx_kill_r += pr.emax * hill(conc, pr.mic_resistant) * shield;
            antitoxin = antitoxin.max(pr.antitoxin * hill(conc, pr.mic_sensitive));
            abx_level += conc / pr.mic_sensitive;
        }
        // Антисептик гасится органикой (гной, некроз) и лишь частично проникает в слаф и биоплёнку.
        let organic = block.max((2.0 * bt).min(1.0) * (1.0 - e));
        let asc_eff = asc * (1.0 - asp.organic * organic);
        let as_kill = asp.kill * asc_eff * (1.0 - (1.0 - asp.biofilm_pen) * bf) * (1.0 - (1.0 - asp.necro_pen) * block);
        // Личинки поедают бактерии вместе с мёртвой тканью.
        let larvae_kill = if gl.larvae { p.larvae_kill * block } else { 0.0 };

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
            / (1.0 + p.leukocidin * bt * (1.0 - antitoxin));
        // Биоплёнка постоянно «выпускает» планктонные бактерии.
        let shed = p.biofilm_shed * bf * crowd;
        let frac_r = if bt > 1e-6 { br / bt } else { 0.0 };
        let d_bact = p.d_bact * s.bacteria[i] + grow_s - mutate - (immune + as_kill + abx_kill_s + larvae_kill) * b
            + shed * (1.0 - frac_r);
        let d_bact_res = p.d_bact * s.bacteria_res[i] + grow_r + mutate - (immune + as_kill + abx_kill_r + larvae_kill) * br
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
            + p.k_toxin * (bt - p.toxin_threshold).max(0.0) * (1.0 - 0.5 * block) * (1.0 - antitoxin)
            + p.k_toxin_field * tox;
        // Экзотоксины выделяют бактерии некротизирующих штаммов; клиндамицин выключает их синтез.
        let d_tox = p.d_toxin * s.toxin[i] + p.toxin_production * bt * (1.0 - antitoxin) - p.toxin_decay * tox;
        let room = (p.max_depth_mm - depth).max(0.0);
        let dead_mm = (death * p.necrosis_depth_mm).min(room / dt);
        // Очищение от некроза. Аутолиз: макрофаги растворяют мёртвую ткань (под влажной повязкой — быстрее).
        // Коллагеназа и личинки переваривают её сами, но хуже справляются с сухим струпом.
        let dryness = (slough / 1.2).min(1.0) * (1.0 - 0.7 * bt.min(1.0)) * if gl.hydrogel { 0.3 } else { 1.0 };
        let mut autolysis = p.slough_clear * (m1 + m2) * sat(slough, 0.5);
        if gl.hydrogel {
            autolysis *= p.hydrogel_boost;
        }
        if gl.collagenase {
            autolysis += p.collagenase_rate * sat(slough, 0.3) * (1.0 - 0.5 * dryness);
        }
        if gl.larvae {
            autolysis += p.larvae_rate * sat(slough, 0.3) * (1.0 - 0.7 * dryness);
        }

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
            - asp.tox_fib * asc * f;
        // --- Глубокие слои. Мышца: клетки-сателлиты строят новые волокна, если дефект небольшой,
        // а воспаление утихает; большой дефект, инфекция, некроз и ишемия ведут к фиброзу.
        let lost_mus = p.lost_muscle_mm(dmax);
        let bed_in_muscle = depth > p.fat_bottom_mm() + 0.3;
        let (mut d_myo, mut d_mscar, mut mus_gate) = (0.0, 0.0, 0.0);
        if lost_mus > 0.0 {
            let capacity = (-lost_mus / p.muscle_loss_scale).exp();
            let free = (1.0 - myo - mscar).max(0.0);
            let gate = (1.0 - block) * oxy * (1.0 - sat(bt, 0.1)) * (0.3 + 0.7 * sat(m2, 0.05));
            let ischemia = ((0.3 - o) / 0.3).max(0.0);
            let drive = (1.0 - capacity) + 0.5 * (bt.min(1.0) + block + ischemia) + 0.1 * sat(m1, 0.2);
            mus_gate = capacity * gate;
            let exposed = if bed_in_muscle { 1.0 } else { 0.0 };
            d_myo = p.muscle_regen * mus_gate * free - asp.tox_fib * asc * myo * exposed;
            d_mscar = p.muscle_fibrosis * drive * free;
        }
        // Жир во взрослой ране почти не возвращается (адипоциты появляются лишь около новых фолликулов).
        let lost_fat = p.lost_fat_mm(dmax);
        let d_fat = if lost_fat > 0.0 { p.fat_regen * e * oxy * (1.0 - block) * (1.0 - fat_new) } else { 0.0 };

        // Грануляции заполняют полость снизу; через мёртвую ткань не прорастают.
        // В мышечной части полость быстрее заполняют и регенерирующие волокна.
        let muscle_fill = if bed_in_muscle { 1.0 + 1.5 * mus_gate } else { 1.0 };
        let fill = if depth > slough { p.fill_rate * f * oxy * (1.0 - block) * muscle_fill } else { 0.0 };

        let deposit = p.collagen_rate * f * oxy * (1.0 - c).max(0.0) * (0.5 + 0.5 * gfh);
        // MMP макрофагов перестраивают только повреждённую дерму.
        let degrade = p.collagen_degr_m1 * m1 * c * (1.0 - p.residual_dermis(dmax)) + 0.5 * death * c;

        let d_vessel = (p.d_vessel * s.vessels[i] + p.angio_rate * sat(a, p.vegf_half) * v * (1.0 - v)) * scaffold
            - 0.5 * asp.tox_fib * asc * v * (1.0 - e);

        // Эпителий наползает только на выровненное грануляциями дно;
        // в неглубоких ранах он прорастает ещё и островками из придатков кожи.
        let level = 1.0 / (1.0 + ((depth - p.epi_level_mm).max(0.0) / 0.4).powi(2));
        let epi_scaffold = (cf + c).min(1.0) * level * (1.0 - block);
        let islands = if wounded { p.residual_dermis(dmax) } else { 0.0 };
        // Во влажной среде (гидрогель) кератиноциты мигрируют быстрее.
        let moist = if gl.hydrogel { 1.3 } else { 1.0 };
        let d_epi = moist * (p.d_epi * s.epithelium[i]
            + p.epi_rate * e * (1.0 - e) * (0.5 + 0.5 * gfh) * oxy / (1.0 + bt / p.epi_bact_half))
            * epi_scaffold
            + p.adnexal_rate * islands * (1.0 - e) * level * (1.0 - block) * oxy / (1.0 + bt / p.epi_bact_half)
            - asp.tox_epi * asc * e;

        // --- Запись.
        ch.bleeding[j] = (bl + dt * d_bleed).clamp(0.0, 1.0);
        ch.clot[j] = (cf + dt * d_clot).clamp(0.0, 1.0);
        ch.debris[j] = (debris + dt * d_debris).max(0.0);
        let b_new = (b + dt * d_bact).max(0.0);
        ch.bacteria[j] = if b_new < 1e-5 { 0.0 } else { b_new };
        let br_new = (br + dt * d_bact_res).max(0.0);
        ch.bacteria_res[j] = if br_new < 1e-9 { 0.0 } else { br_new };
        ch.biofilm[j] = (bf + dt * d_biofilm).clamp(0.0, 1.0);
        ch.oxygen[j] = (o + dt * d_o2).clamp(0.0, 1.0);
        ch.signal[j] = (sig + dt * d_sig).max(0.0);
        ch.growth_factor[j] = (g + dt * d_gf).max(0.0);
        ch.vegf[j] = (a + dt * d_vegf).max(0.0);
        ch.neutrophils[j] = (n + dt * d_neut).max(0.0);
        ch.m1[j] = (m1 + dt * d_m1).max(0.0);
        ch.m2[j] = (m2 + dt * d_m2).max(0.0);
        ch.fibroblasts[j] = ((f + dt * d_fib) * (1.0 - dt * death)).max(0.0);
        ch.vessels[j] = ((v + dt * d_vessel) * (1.0 - dt * death)).clamp(0.0, 1.0);
        ch.epithelium[j] = ((e + dt * d_epi) * (1.0 - dt * death)).clamp(0.0, 1.0);
        ch.antibiotic[j] = abx_level;
        ch.toxin[j] = (tox + dt * d_tox).max(0.0);

        let slough_new = (slough + dt * (dead_mm - autolysis)).max(0.0);
        let depth_new = (depth + dt * (dead_mm - fill)).clamp(0.0, p.max_depth_mm).max(slough_new);
        ch.slough[j] = slough_new;
        ch.depth[j] = depth_new;
        let dmax_new = dmax.max(depth_new);
        ch.depth_max[j] = dmax_new;
        ch.myo[j] = rescale((myo + dt * d_myo).clamp(0.0, 1.0), lost_mus, p.lost_muscle_mm(dmax_new));
        ch.muscle_scar[j] = rescale((mscar + dt * d_mscar).clamp(0.0, 1.0), lost_mus, p.lost_muscle_mm(dmax_new));
        ch.fat_new[j] = rescale((fat_new + dt * d_fat).clamp(0.0, 1.0), lost_fat, p.lost_fat_mm(dmax_new));
        if slough_new > 0.3 || depth_new > 0.3 {
            ch.wound_mask[j] = true;
        }

        // --- Коллаген и ремоделирование: новый коллаген незрелый (тип III),
        // со временем созревает. В полнослойном рубце — не выше scar_maturity_cap,
        // при уцелевшей дерме заживает почти без рубца.
        let c_new = (c + dt * (deposit - degrade)).clamp(0.0, 1.0);
        let mature = (q * c * (1.0 - dt * (degrade / c.max(1e-6)))).max(0.0);
        let mut q_new = if c_new > 1e-6 { (mature / c_new).min(1.0) } else { 0.0 };
        let cap = p.scar_maturity_cap + (1.0 - p.scar_maturity_cap) * p.residual_dermis(dmax);
        q_new += dt * p.maturation_rate * (cap - q_new).max(0.0) / (1.0 + 5.0 * m1);
        ch.collagen[j] = c_new;
        ch.maturity[j] = q_new;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::Scenario;
    use crate::tissue::WoundShape;

    const W: usize = 96;
    const H: usize = 48;

    fn run(t: &mut Tissue, p: &Params, steps: usize) {
        let mut s = Scratch::new(t.len());
        for _ in 0..steps {
            step(t, p, &mut s);
        }
    }

    #[test]
    fn hill_curve() {
        assert_eq!(hill(0.0, 1.0), 0.0);
        assert!((hill(1.0, 1.0) - 0.5).abs() < 1e-6);
        assert!(hill(4.0, 1.0) > hill(2.0, 1.0));
        assert!(hill(100.0, 1.0) > 0.99);
    }

    #[test]
    fn healthy_tissue_is_stationary() {
        let p = Params::default();
        let mut t = Tissue::healthy(W, H, &p);
        run(&mut t, &p, 500);
        for i in [0, (H / 2) * W + W / 2, t.len() - 1] {
            assert!((t.epithelium.data[i] - 1.0).abs() < 1e-4);
            assert!((t.collagen.data[i] - 1.0).abs() < 1e-4);
            assert!((t.vessels.data[i] - 1.0).abs() < 1e-4);
            assert!((t.fibroblasts.data[i] - p.fib_baseline).abs() < 1e-4);
            assert_eq!(t.bacteria_total(i), 0.0);
            assert!(t.neutrophils.data[i] < 1e-6);
            assert_eq!(t.depth.data[i], 0.0);
            assert_eq!(t.slough.data[i], 0.0);
            assert!(t.oxygen.data[i] > 0.85, "здоровая ткань не должна задыхаться");
        }
    }

    #[test]
    fn fields_stay_physical_during_infected_healing() {
        let p = Scenario::Necrotizing.params();
        let mut t = Tissue::healthy(W, H, &p);
        t.injure(WoundShape::Circle { radius: 16.0 }, 5.0, &p);
        run(&mut t, &p, 24 * 10 * 4);
        let unit = [&t.epithelium, &t.collagen, &t.vessels, &t.bleeding, &t.clot, &t.oxygen, &t.biofilm, &t.myo, &t.fat_new];
        for i in 0..t.len() {
            for f in unit {
                let v = f.data[i];
                assert!(v.is_finite() && (-1e-6..=1.0 + 1e-6).contains(&v), "значение вне 0..1: {v}");
            }
            for f in [&t.bacteria, &t.bacteria_res, &t.neutrophils, &t.m1, &t.m2, &t.fibroblasts, &t.signal, &t.debris] {
                assert!(f.data[i].is_finite() && f.data[i] >= 0.0);
            }
            assert!(t.slough.data[i] <= t.depth.data[i] + 1e-5, "мёртвая ткань не может быть толще полости");
            assert!(t.depth.data[i] <= t.depth_max.data[i] + 1e-5);
            assert!(t.depth_max.data[i] <= p.max_depth_mm + 1e-5);
            assert!(t.myo.data[i] + t.muscle_scar.data[i] <= 1.0 + 1e-4);
        }
    }

    #[test]
    fn parallel_step_is_deterministic() {
        let p = Scenario::Infected.params();
        let mut a = Tissue::healthy(W, H, &p);
        a.injure(WoundShape::Circle { radius: 12.0 }, 2.5, &p);
        let mut b = Tissue::healthy(W, H, &p);
        b.injure(WoundShape::Circle { radius: 12.0 }, 2.5, &p);
        run(&mut a, &p, 300);
        run(&mut b, &p, 300);
        assert_eq!(a.bacteria.data, b.bacteria.data);
        assert_eq!(a.epithelium.data, b.epithelium.data);
        assert_eq!(a.depth.data, b.depth.data);
    }

    #[test]
    fn clindamycin_suppresses_toxic_necrosis() {
        // Одинаковая нагрузка токсинпродуцирующих бактерий; клиндамицин в плазме снижает гибель ткани.
        let p = Scenario::Necrotizing.params();
        let necrosis = |clinda: bool| {
            let mut t = Tissue::healthy(W, H, &p);
            t.injure(WoundShape::Circle { radius: 12.0 }, 2.5, &p);
            for i in 0..t.len() {
                if t.wound_mask[i] {
                    t.bacteria.data[i] = 0.9;
                }
            }
            if clinda {
                t.abx_plasma[Antibiotic::Clindamycin.index()] = 6.0;
            }
            run(&mut t, &p, 5);
            t.slough.data.iter().sum::<f32>()
        };
        assert!(necrosis(true) < 0.8 * necrosis(false));
    }

    #[test]
    fn hydrogel_and_larvae_speed_up_clearing_of_dead_tissue() {
        let p = Params::default();
        let remaining = |debriders: [bool; 3]| {
            let mut t = Tissue::healthy(W, H, &p);
            t.injure(WoundShape::Circle { radius: 12.0 }, 2.5, &p);
            for i in 0..t.len() {
                if t.wound_mask[i] {
                    t.slough.data[i] = 1.0;
                    t.depth.data[i] = t.depth.data[i].max(1.0);
                }
            }
            t.debriders = debriders;
            run(&mut t, &p, 240);
            t.slough.data.iter().sum::<f32>()
        };
        let none = remaining([false; 3]);
        let hydrogel = remaining([true, false, false]);
        let larvae = remaining([false, false, true]);
        assert!(hydrogel < none, "гидрогель ускоряет аутолиз");
        assert!(larvae < hydrogel, "личинки быстрее гидрогеля");
    }

    #[test]
    fn antiseptic_is_toxic_to_new_epithelium() {
        use crate::therapy::Antiseptic;
        let p = Params::default();
        let epi = |agent: Option<Antiseptic>| {
            let mut t = Tissue::healthy(W, H, &p);
            t.injure(WoundShape::Circle { radius: 12.0 }, 0.1, &p);
            if let Some(a) = agent {
                t.antiseptic_agent = a;
                t.antiseptic.data.iter_mut().for_each(|v| *v = 1.0);
            }
            run(&mut t, &p, 60);
            t.epithelium.data.iter().sum::<f32>()
        };
        let clean = epi(None);
        let peroxide = epi(Some(Antiseptic::Peroxide));
        let hocl = epi(Some(Antiseptic::Hypochlorous));
        assert!(peroxide < hocl && hocl <= clean + 1e-3, "{peroxide} < {hocl} <= {clean}");
    }
}