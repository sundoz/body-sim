use std::fs::File;
use std::io::{BufWriter, Write};
use std::process::ExitCode;

use body_sim::params::Scenario;
use body_sim::report::{self, Condition, Metrics, Phase, View};
use body_sim::simulation::Simulation;
use body_sim::therapy::{Antibiotic, Antiseptic, Debrider};
use body_sim::tissue::WoundShape;

const HELP: &str = "\
body-sim — симулятор заживления кожной раны (консольный режим)

Использование: body-sim [опции]

Пациент и рана:
  --scenario <s>   healthy | diabetic | infected | elderly |
                   ischemic | diabetic-foot | necrotizing      (healthy)
  --wound <w>      circle | cut                              (circle)
  --size <мм>      радиус круглой раны / четверть длины разреза (4)
  --depth <мм>     глубина: 0.1 эпидермис, 1 дерма, 2.5 полнослойная,
                   5 жировая клетчатка, 9 мышца                (2.5)

Лечение:
  --antiseptic <a>        octenidine | polyhexanide | hypochlorous |
                          chlorhexidine | povidone | peroxide
  --antiseptic-every <ч>  интервал перевязок                  (24)
  --antibiotic-every <ч>  курс цефазолина, доза каждые N часов
  --antibiotics <список>  несколько препаратов: cefazolin:8,clindamycin:8,vancomycin:12
  --debriders <список>    очищение от некроза: hydrogel,collagenase,larvae
  --treat-from <день>     когда начать лечение                (0)
  --debride <дни>         хирургическая обработка, напр. 3,10

Вывод:
  --days <n>       длительность моделирования в днях         (28)
  --every <n>      печатать срез каждые n дней                (2)
  --view <поле>    что рисовать на карте                      (integrity)
                   integrity, strength, epithelium, collagen, vessels, oxygen,
                   clot, depth, necrosis, bacteria, resistant, biofilm, toxin,
                   neutrophils, macrophages, fibroblasts, antiseptic, antibiotic
  --no-map         не рисовать карты, только цифры
  --csv <файл>     записать почасовую динамику в CSV
  --grid <WxH>     размер сетки в клетках по 0.25 мм         (96x48)
  -h, --help       эта справка
";

struct Cli {
    scenario: Scenario,
    wound: String,
    size_mm: f32,
    depth_mm: f32,
    antiseptic: Option<Antiseptic>,
    antiseptic_every: f32,
    antibiotics: Vec<(Antibiotic, f32)>,
    debriders: Vec<Debrider>,
    treat_from: f32,
    debride: Vec<f32>,
    days: f32,
    every: f32,
    view: View,
    map: bool,
    csv: Option<String>,
    w: usize,
    h: usize,
}

fn parse_args(argv: impl IntoIterator<Item = String>) -> Result<Option<Cli>, String> {
    let mut cli = Cli {
        scenario: Scenario::Healthy,
        wound: "circle".into(),
        size_mm: 4.0,
        depth_mm: 2.5,
        antiseptic: None,
        antiseptic_every: 24.0,
        antibiotics: Vec::new(),
        debriders: Vec::new(),
        treat_from: 0.0,
        debride: Vec::new(),
        days: 28.0,
        every: 2.0,
        view: View::Integrity,
        map: true,
        csv: None,
        w: 96,
        h: 48,
    };
    let mut args = argv.into_iter();
    while let Some(a) = args.next() {
        let mut val = || args.next().ok_or_else(|| format!("для {a} нужно значение"));
        match a.as_str() {
            "-h" | "--help" => return Ok(None),
            "--scenario" => {
                let v = val()?;
                cli.scenario = Scenario::parse(&v).ok_or_else(|| format!("неизвестный сценарий: {v}"))?;
            }
            "--wound" => {
                let v = val()?;
                if v != "circle" && v != "cut" {
                    return Err(format!("неизвестная форма раны: {v}"));
                }
                cli.wound = v;
            }
            "--size" => cli.size_mm = num(&val()?)?,
            "--depth" => cli.depth_mm = num(&val()?)?,
            "--antiseptic" => {
                let v = val()?;
                cli.antiseptic = Some(Antiseptic::parse(&v).ok_or_else(|| format!("неизвестный антисептик: {v}"))?);
            }
            "--antiseptic-every" => cli.antiseptic_every = num(&val()?)?,
            "--antibiotic-every" => cli.antibiotics.push((Antibiotic::Cefazolin, num(&val()?)?)),
            "--antibiotics" => {
                for item in val()?.split(',') {
                    let (name, every) = item.split_once(':').ok_or("формат --antibiotics: cefazolin:8")?;
                    let drug =
                        Antibiotic::parse(name.trim()).ok_or_else(|| format!("неизвестный антибиотик: {name}"))?;
                    cli.antibiotics.push((drug, num(every.trim())?));
                }
            }
            "--debriders" => {
                for name in val()?.split(',') {
                    cli.debriders.push(
                        Debrider::parse(name.trim()).ok_or_else(|| format!("неизвестный метод очищения: {name}"))?,
                    );
                }
            }
            "--treat-from" => cli.treat_from = non_negative(&val()?, "--treat-from")?,
            "--debride" => {
                let mut days: Vec<f32> =
                    val()?.split(',').map(|s| non_negative(s.trim(), "--debride")).collect::<Result<_, _>>()?;
                // Обработки выполняются строго по порядку списка, поэтому неотсортированный
                // список раньше молча терял ранние операции: `3,1` вело себя как `3,3`.
                days.sort_by(f32::total_cmp);
                days.dedup();
                cli.debride = days;
            }
            "--days" => cli.days = num(&val()?)?,
            "--every" => cli.every = num(&val()?)?,
            "--view" => {
                let v = val()?;
                cli.view = View::parse(&v).ok_or_else(|| format!("неизвестное поле: {v} (есть: {})", View::NAMES))?;
            }
            "--no-map" => cli.map = false,
            "--csv" => cli.csv = Some(val()?),
            "--grid" => {
                let v = val()?;
                let (w, h) = v.split_once('x').ok_or("формат --grid: 96x48")?;
                cli.w = w.parse().map_err(|_| "формат --grid: 96x48")?;
                cli.h = h.parse().map_err(|_| "формат --grid: 96x48")?;
                if cli.w == 0 || cli.h == 0 {
                    return Err("--grid: размеры сетки должны быть больше нуля".into());
                }
            }
            _ => return Err(format!("неизвестная опция: {a}")),
        }
    }
    Ok(Some(cli))
}

fn num(s: &str) -> Result<f32, String> {
    s.parse::<f32>()
        .ok()
        .filter(|v| *v > 0.0 && v.is_finite())
        .ok_or_else(|| format!("ожидалось положительное число, получено: {s}"))
}

/// Число ≥ 0 (день начала лечения, день обработки: ноль — «сразу»).
fn non_negative(s: &str, opt: &str) -> Result<f32, String> {
    s.parse::<f32>()
        .ok()
        .filter(|v| *v >= 0.0 && v.is_finite())
        .ok_or_else(|| format!("{opt}: ожидалось неотрицательное число, получено: {s}"))
}

/// Момент и величина пика какой-либо популяции.
#[derive(Default)]
struct Peak {
    value: f32,
    hours: f32,
}

impl Peak {
    fn track(&mut self, value: f32, hours: f32) {
        if value > self.value {
            self.value = value;
            self.hours = hours;
        }
    }
}

fn main() -> ExitCode {
    let cli = match parse_args(std::env::args().skip(1)) {
        Ok(Some(c)) => c,
        Ok(None) => {
            print!("{HELP}");
            return ExitCode::SUCCESS;
        }
        Err(e) => {
            eprintln!("ошибка: {e}\n\n{HELP}");
            return ExitCode::FAILURE;
        }
    };

    let p = cli.scenario.params();
    if cli.depth_mm > p.max_depth_mm {
        eprintln!(
            "ошибка: --depth {:.1} мм глубже моделируемой колонки ткани ({:.1} мм)",
            cli.depth_mm, p.max_depth_mm
        );
        return ExitCode::FAILURE;
    }
    let size_cells = cli.size_mm / p.cell_mm;
    let shape = if cli.wound == "cut" {
        WoundShape::Cut { half_length: 2.0 * size_cells, half_width: 1.0 / p.cell_mm }
    } else {
        WoundShape::Circle { radius: size_cells }
    };
    // Рана наносится в центр участка; вылезая за край, она молча обрезается,
    // а условие непротекания превращает обрезанный край в «бесконечную» рану.
    let (need_w, need_h) = match shape {
        WoundShape::Cut { half_length, half_width } => (2.0 * half_length, 2.0 * half_width),
        WoundShape::Circle { radius } => (2.0 * radius, 2.0 * radius),
    };
    if need_w > cli.w as f32 || need_h > cli.h as f32 {
        eprintln!(
            "внимание: рана ({:.0}×{:.0} клеток) не помещается в сетку {}×{} и будет обрезана —              увеличьте --grid или уменьшите --size",
            need_w, need_h, cli.w, cli.h
        );
    }
    let layer = p.layer_title(cli.depth_mm);
    let mut sim = Simulation::new(p, cli.w, cli.h, shape, cli.depth_mm);

    // Ошибки записи CSV не должны теряться: иначе при полном диске файл молча окажется обрезанным.
    let mut csv_error: Option<std::io::Error> = None;
    let mut csv = match &cli.csv {
        Some(path) => match File::create(path) {
            Ok(f) => {
                let mut w = BufWriter::new(f);
                if let Err(e) = writeln!(w, "{}", Metrics::CSV_HEADER) {
                    csv_error = Some(e);
                }
                Some(w)
            }
            Err(e) => {
                eprintln!("не удалось создать {path}: {e}");
                return ExitCode::FAILURE;
            }
        },
        None => None,
    };

    let per_hour = sim.steps_per_hour();
    let total_hours = (cli.days * 24.0).round() as usize;
    let snap_every = ((cli.every * 24.0).round() as usize).max(1);

    let initial = sim.latest().clone();
    println!(
        "Сценарий: {} | рана: {} {:.1} мм, глубина {:.1} мм ({}), площадь {:.1} мм² | {} дней",
        cli.scenario.title(),
        if cli.wound == "cut" { "разрез, полудлина" } else { "круглая, радиус" },
        if cli.wound == "cut" { 2.0 * cli.size_mm } else { cli.size_mm },
        cli.depth_mm,
        layer,
        initial.wound_mm2,
        cli.days
    );
    let mut plan = Vec::new();
    if let Some(a) = cli.antiseptic {
        plan.push(format!("{} каждые {:.0} ч", a.title(), cli.antiseptic_every));
    }
    for (drug, h) in &cli.antibiotics {
        plan.push(format!("{} каждые {h:.0} ч", drug.title()));
    }
    for d in &cli.debriders {
        plan.push(d.title().to_lowercase());
    }
    if !plan.is_empty() {
        println!("Лечение с {:.1} дня: {}", cli.treat_from, plan.join(", "));
    }
    if !cli.debride.is_empty() {
        let days: Vec<String> = cli.debride.iter().map(|d| format!("{d}")).collect();
        println!("Хирургическая обработка на дни: {}", days.join(", "));
    }
    if cli.map {
        println!("Карта: {:?}, шкала ' .:-=+*#%@' = 0 … {}", cli.view, cli.view.scale_label());
    }

    let mut phase_starts: Vec<(Phase, f32)> = Vec::new();
    let mut conditions: Vec<(Condition, f32)> = Vec::new();
    let mut worst = Condition::Healed;
    let [mut pk_neut, mut pk_bact, mut pk_necr, mut pk_inf, mut pk_depth] =
        std::array::from_fn::<Peak, 5, _>(|_| Peak::default());
    let mut treatment_started = false;
    let mut debride_idx = 0;

    for hour in 0..=total_hours {
        let hours = hour as f32;
        if !treatment_started && hours >= cli.treat_from * 24.0 {
            treatment_started = true;
            if let Some(a) = cli.antiseptic {
                sim.therapy.start_antiseptic(a, cli.antiseptic_every, hours);
            }
            for (drug, every) in &cli.antibiotics {
                sim.therapy.start_antibiotic(*drug, *every, hours);
            }
            for d in &cli.debriders {
                sim.tissue.debriders[d.index()] = true;
            }
        }
        while debride_idx < cli.debride.len() && hours >= cli.debride[debride_idx] * 24.0 {
            sim.therapy.debride_now(&mut sim.tissue, &sim.p, hours);
            debride_idx += 1;
        }

        let m = sim.latest().clone();
        if let Some(w) = csv.as_mut() {
            if let (Err(e), None) = (writeln!(w, "{}", m.csv_row()), &csv_error) {
                csv_error = Some(e);
            }
        }
        let phase = m.phase();
        if phase_starts.last().map(|(ph, _)| *ph) != Some(phase) {
            phase_starts.push((phase, hours));
        }
        let cond = m.condition();
        if conditions.last().map(|(c, _)| *c) != Some(cond) {
            conditions.push((cond, hours));
        }
        if cond > worst {
            worst = cond;
        }
        pk_neut.track(m.neutrophils, hours);
        pk_bact.track(m.bacteria, hours);
        pk_necr.track(m.necrotic_mm2, hours);
        pk_inf.track(m.infected_mm2, hours);
        pk_depth.track(m.depth_max, hours);

        if hour % snap_every == 0 || hour == total_hours {
            let map = cli.map.then(|| report::render_map(&sim.tissue, cli.view));
            report::print_snapshot(&m, map.as_deref());
        }
        if hour < total_hours {
            sim.run_steps(per_hour);
        }
    }

    if let Some(w) = csv.as_mut() {
        if let (Err(e), None) = (w.flush(), &csv_error) {
            csv_error = Some(e);
        }
    }

    let last = sim.latest().clone();
    println!("\n=== Итог ===");
    match sim.closed_at {
        Some(h) => println!("  Рана эпителизировалась на {:.1} день", h / 24.0),
        None => println!("  Рана НЕ закрылась за {} дней (открыто {:.1} мм²)", cli.days, last.open_mm2),
    }
    let fmt = |v: &[(String, f32)]| {
        v.iter().map(|(s, h)| format!("{s} с {:.1} дн.", h / 24.0)).collect::<Vec<_>>().join(" → ")
    };
    let ph: Vec<(String, f32)> = phase_starts.iter().map(|(p, h)| (p.title().to_string(), *h)).collect();
    let co: Vec<(String, f32)> = conditions.iter().map(|(c, h)| (c.title().to_string(), *h)).collect();
    println!("  Фазы: {}", fmt(&ph));
    println!("  Состояние: {}", fmt(&co));
    println!("  Худшее состояние: {}; итоговое: {}", worst.title(), last.condition().title());
    println!("  Максимальная глубина поражения: {:.1} мм", pk_depth.value);
    println!("  Пик некроза: {:.1} мм² на {:.1} день", pk_necr.value, pk_necr.hours / 24.0);
    println!("  Пик инфицированной площади: {:.1} мм² на {:.1} день", pk_inf.value, pk_inf.hours / 24.0);
    println!("  Пик бактерий в ране: {:.2} на {:.1} день", pk_bact.value, pk_bact.hours / 24.0);
    println!("  Пик нейтрофилов: {:.2} на {:.1} день", pk_neut.value, pk_neut.hours / 24.0);
    println!(
        "  Бактерии в конце: {:.3}, из них устойчивых {:.0}%, биоплёнка {:.2}",
        last.bacteria,
        last.resistant_fraction * 100.0,
        last.biofilm
    );
    if last.lost_fat_mm > 0.0 {
        println!(
            "  Клетчатка: утрачено {:.1} мм, вернулось жиром {:.0}% — остальное фиброзный рубец",
            last.lost_fat_mm,
            last.fat_regen * 100.0
        );
    }
    if last.lost_muscle_mm > 0.0 {
        println!(
            "  Мышца: утрачено {:.1} мм, новые волокна {:.0}%, фиброз {:.0}%, ещё не заполнено {:.0}%",
            last.lost_muscle_mm,
            last.muscle_regen * 100.0,
            last.muscle_fibrosis * 100.0,
            (1.0 - last.muscle_regen - last.muscle_fibrosis).max(0.0) * 100.0
        );
    }
    println!(
        "  Прочность рубца: {:.0}% от здоровой кожи (коллаген {:.2}, зрелость {:.2})",
        last.strength * 100.0,
        last.collagen,
        last.maturity
    );
    if let Some(path) = &cli.csv {
        if let Some(e) = &csv_error {
            eprintln!("ошибка записи CSV {path}: {e}");
            return ExitCode::FAILURE;
        }
        println!("  CSV: {path}");
    }
    ExitCode::SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(args: &[&str]) -> Result<Option<Cli>, String> {
        parse_args(args.iter().map(|s| s.to_string()))
    }

    #[test]
    fn defaults() {
        let cli = parse(&[]).unwrap().unwrap();
        assert_eq!(cli.scenario, Scenario::Healthy);
        assert_eq!(cli.depth_mm, 2.5);
        assert_eq!(cli.days, 28.0);
        assert!(cli.antibiotics.is_empty() && cli.debriders.is_empty() && cli.antiseptic.is_none());
        assert!(cli.map);
    }

    #[test]
    fn help() {
        assert!(parse(&["--help"]).unwrap().is_none());
    }

    #[test]
    fn full_treatment_plan() {
        let cli = parse(&[
            "--scenario",
            "necrotizing",
            "--depth",
            "5",
            "--antiseptic",
            "phmb",
            "--antibiotic-every",
            "8",
            "--antibiotics",
            "clindamycin:8, vancomycin:12",
            "--debriders",
            "larvae,hydrogel",
            "--debride",
            "0.5, 3",
            "--treat-from",
            "0.5",
            "--no-map",
        ])
        .unwrap()
        .unwrap();
        assert_eq!(cli.scenario, Scenario::Necrotizing);
        assert_eq!(cli.antiseptic, Some(Antiseptic::Polyhexanide));
        assert_eq!(
            cli.antibiotics,
            vec![(Antibiotic::Cefazolin, 8.0), (Antibiotic::Clindamycin, 8.0), (Antibiotic::Vancomycin, 12.0)]
        );
        assert_eq!(cli.debriders, vec![Debrider::Larvae, Debrider::Hydrogel]);
        assert_eq!(cli.debride, vec![0.5, 3.0]);
        assert_eq!(cli.treat_from, 0.5);
        assert!(!cli.map);
    }

    #[test]
    fn grid_and_view() {
        let cli = parse(&["--grid", "40x20", "--view", "toxin"]).unwrap().unwrap();
        assert_eq!((cli.w, cli.h), (40, 20));
        assert_eq!(cli.view, View::Toxin);
    }

    #[test]
    fn errors_are_reported() {
        for bad in [
            vec!["--scenario", "zombie"],
            vec!["--antibiotics", "aspirin:8"],
            vec!["--antibiotics", "cefazolin"],
            vec!["--debriders", "leeches"],
            vec!["--depth", "-1"],
            vec!["--days"],
            vec!["--grid", "40"],
            vec!["--grid", "0x48"],
            vec!["--grid", "96x0"],
            vec!["--treat-from", "-1"],
            vec!["--debride", "1,-2"],
            vec!["--wat"],
        ] {
            assert!(parse(&bad).is_err(), "{bad:?} должно быть ошибкой");
        }
    }

    #[test]
    fn debridement_days_are_sorted_and_deduplicated() {
        // Обработки выполняются по порядку списка, поэтому «3,1» раньше молча вело себя как «3,3».
        let cli = parse(&["--debride", "3, 1, 10, 3"]).unwrap().unwrap();
        assert_eq!(cli.debride, vec![1.0, 3.0, 10.0]);
        assert_eq!(parse(&["--debride", "0"]).unwrap().unwrap().debride, vec![0.0], "ноль — «сразу»");
    }

    #[test]
    fn treatment_may_start_on_day_zero() {
        assert_eq!(parse(&["--treat-from", "0"]).unwrap().unwrap().treat_from, 0.0);
    }
}
