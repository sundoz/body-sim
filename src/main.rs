use std::fs::File;
use std::io::{BufWriter, Write};
use std::process::ExitCode;

use body_sim::params::Scenario;
use body_sim::report::{self, Condition, Metrics, Phase, View};
use body_sim::simulation::Simulation;
use body_sim::therapy::Antiseptic;
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
  --antiseptic <a>        octenidine | chlorhexidine | povidone | peroxide
  --antiseptic-every <ч>  интервал перевязок                  (24)
  --antibiotic-every <ч>  курс антибиотика, доза каждые N часов
  --treat-from <день>     когда начать лечение                (0)
  --debride <дни>         хирургическая обработка, напр. 3,10

Вывод:
  --days <n>       длительность моделирования в днях         (28)
  --every <n>      печатать срез каждые n дней                (2)
  --view <поле>    что рисовать на карте                      (integrity)
                   integrity, strength, epithelium, collagen, vessels, oxygen,
                   clot, depth, necrosis, bacteria, resistant, biofilm,
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
    antibiotic_every: Option<f32>,
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

fn parse_args() -> Result<Option<Cli>, String> {
    let mut cli = Cli {
        scenario: Scenario::Healthy,
        wound: "circle".into(),
        size_mm: 4.0,
        depth_mm: 2.5,
        antiseptic: None,
        antiseptic_every: 24.0,
        antibiotic_every: None,
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
    let mut args = std::env::args().skip(1);
    while let Some(a) = args.next() {
        let mut val = || args.next().ok_or_else(|| format!("для {a} нужно значение"));
        match a.as_str() {
            "-h" | "--help" => return Ok(None),
            "--scenario" => {
                let v = val()?;
                cli.scenario =
                    Scenario::parse(&v).ok_or_else(|| format!("неизвестный сценарий: {v}"))?;
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
                cli.antiseptic =
                    Some(Antiseptic::parse(&v).ok_or_else(|| format!("неизвестный антисептик: {v}"))?);
            }
            "--antiseptic-every" => cli.antiseptic_every = num(&val()?)?,
            "--antibiotic-every" => cli.antibiotic_every = Some(num(&val()?)?),
            "--treat-from" => cli.treat_from = val()?.parse().map_err(|_| "--treat-from: ожидалось число")?,
            "--debride" => {
                cli.debride = val()?
                    .split(',')
                    .map(|s| s.trim().parse::<f32>().map_err(|_| format!("--debride: плохой день {s}")))
                    .collect::<Result<_, _>>()?;
            }
            "--days" => cli.days = num(&val()?)?,
            "--every" => cli.every = num(&val()?)?,
            "--view" => {
                let v = val()?;
                cli.view = View::parse(&v)
                    .ok_or_else(|| format!("неизвестное поле: {v} (есть: {})", View::NAMES))?;
            }
            "--no-map" => cli.map = false,
            "--csv" => cli.csv = Some(val()?),
            "--grid" => {
                let v = val()?;
                let (w, h) = v.split_once('x').ok_or("формат --grid: 96x48")?;
                cli.w = w.parse().map_err(|_| "формат --grid: 96x48")?;
                cli.h = h.parse().map_err(|_| "формат --grid: 96x48")?;
            }
            _ => return Err(format!("неизвестная опция: {a}")),
        }
    }
    Ok(Some(cli))
}

fn num(s: &str) -> Result<f32, String> {
    s.parse::<f32>()
        .ok()
        .filter(|v| *v > 0.0)
        .ok_or_else(|| format!("ожидалось положительное число, получено: {s}"))
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
    let cli = match parse_args() {
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
    let size_cells = cli.size_mm / p.cell_mm;
    let shape = if cli.wound == "cut" {
        WoundShape::Cut { half_length: 2.0 * size_cells, half_width: 1.0 / p.cell_mm }
    } else {
        WoundShape::Circle { radius: size_cells }
    };
    let layer = p.layer_title(cli.depth_mm);
    let mut sim = Simulation::new(p, cli.w, cli.h, shape, cli.depth_mm);

    let mut csv = match &cli.csv {
        Some(path) => match File::create(path) {
            Ok(f) => {
                let mut w = BufWriter::new(f);
                let _ = writeln!(w, "{}", Metrics::CSV_HEADER);
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
    if let Some(h) = cli.antibiotic_every {
        plan.push(format!("антибиотик каждые {h:.0} ч"));
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
            if let Some(every) = cli.antibiotic_every {
                sim.therapy.start_antibiotic(every, hours);
            }
        }
        while debride_idx < cli.debride.len() && hours >= cli.debride[debride_idx] * 24.0 {
            sim.therapy.debride_now(&mut sim.tissue, &sim.p, hours);
            debride_idx += 1;
        }

        let m = sim.latest().clone();
        if let Some(w) = csv.as_mut() {
            let _ = writeln!(w, "{}", m.csv_row());
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
        let _ = w.flush();
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
    println!(
        "  Прочность рубца: {:.0}% от здоровой кожи (коллаген {:.2}, зрелость {:.2})",
        last.strength * 100.0,
        last.collagen,
        last.maturity
    );
    if let Some(path) = &cli.csv {
        println!("  CSV: {path}");
    }
    ExitCode::SUCCESS
}
