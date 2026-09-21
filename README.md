# body-sim — a human tissue regeneration simulator

[![CI](https://github.com/sundoz/body-sim/actions/workflows/ci.yml/badge.svg)](https://github.com/sundoz/body-sim/actions/workflows/ci.yml)

Skin wound healing on a 24×12 mm patch: depth down to muscle, necrosis, infection,
treatment with antiseptics, antibiotics and surgical debridement. Rust, graphics on `macroquad`.

> The user interface, code comments and CLI output are in Russian; this README is in English.

![Necrotizing infection on day 3](docs/screenshot-necrotizing.png)

## Graphical mode

```bash
cargo run --release
```

- **A 3D tissue block styled like an anatomical atlas.** A relief surface: the wound is a real
  depression filled with blood, slough or granulation tissue; skin with pores and skin markings
  (a scar has neither), the wet sheen of a wound, soft subsurface scattering, tone mapping.
  The faces of the block are sections through the layers with callouts: epidermis, dermis with
  hair follicles, sweat glands and capillaries, fat lobules, fascia, muscle fibre bundles;
  inside the wound — clot, dead tissue, biofilm, individual bacteria and neutrophils.
  Plus heat maps of 18 fields drawn straight onto the surface.
- **Mouse**: left button — inflict a wound of the selected depth, right button — orbit the camera,
  wheel — zoom, Shift+wheel — brush size, ↑/↓ — move the cutting plane, Tab — section on/off.
  Hovering shows every value at that point.
- **"Body" tab**: 8 body sites (each with its own layer thicknesses, local blood supply,
  appendage density and skin flora), the size of the simulated skin patch, and 7 scenarios.
- **"Wound" tab**: wound shape, size and depth.
- **"Treatment" tab**: 6 antiseptics with their profile (potency, penetration into biofilm and
  necrosis, toxicity to epithelium and fibroblasts) — a single application or dressing changes
  every 12/24 h; 3 systemic antibiotics that can be combined (a dose every 8/12/24 h, plasma
  concentrations); debridement — surgery, hydrogel, collagenase, maggots.
- **Header**: healing phase and clinical condition — healing / chronic / necrosis /
  spreading infection / sepsis.
- **Chart**: 14 series (toggled by clicking the legend), phase and condition bands, procedure marks.
- Space — pause, R — restart, → — step one hour while paused.
- Debug options: `--scenario necrotizing --depth 5 --days 3 --cam 0.5,0.7,20 --top --no-cut
  --treat --antibiotic-every 8 --debride 0.5,1.5 --screenshot out.png`.
- Performance: `--bench 300` — 300 frames at maximum speed with a per-section timing table.
  The model step and section rendering run in parallel (`rayon`), and the scene is rebuilt only
  when the tissue changes; on a Ryzen 5 5600H — ~13 ms per frame at 10 days/s (60 FPS),
  90 simulated days in the console — 3.7 s.

## Console mode

```bash
cargo run --release --bin body-sim -- --scenario diabetic-foot --depth 5 --days 42 --no-map
cargo run --release --bin body-sim -- --scenario necrotizing --antibiotics cefazolin:8,clindamycin:8 --treat-from 0.5 --debride 0.5,1.5,3
cargo run --release --bin body-sim -- --scenario infected --antiseptic hypochlorous --antiseptic-every 12
cargo run --release --bin body-sim -- --scenario diabetic-foot --antibiotics vancomycin:12 --debriders larvae,hydrogel
cargo run --release --bin body-sim -- --site shin --lesion 12 --days 60 --no-map
cargo run --release --bin body-sim -- --site sacrum --scenario diabetic --depth 4 --days 60
cargo run --release --bin body-sim -- --area 60x40 --size 12 --days 40 --no-map   # a larger patch
cargo run --release --bin body-sim -- --depth 9 --days 60 --no-map     # muscle and fat regeneration
cargo run --release --bin body-sim -- --help
```

## Tests

```bash
cargo test
```

106 tests, ~2.5 min; CI (GitHub Actions) runs them on Ubuntu and Windows on every push,
plus `cargo fmt --check` and `cargo clippy -D warnings`.
Unit test coverage: `cargo llvm-cov --release --lib --bins` (~20 s) — the model core is 83–100%,
58% overall (the window, the buttons and the 3D rendering need a GPU and are not covered).

What is checked: unit tests (the Laplacian and conservation of mass, the anatomy of a tissue
column, inflicting a wound, the stationarity of healthy tissue, field bounds, determinism of the
parallel step, pharmacokinetics, antiseptics, surgery, metrics and condition classification,
argument parsing, section rendering) and the clinical regression in `tests/clinical.rs` — the model
is required to reproduce known outcomes: healing times by depth, muscle regeneration and fibrosis
after a large defect, the near-absence of fat regeneration, dry necrosis under ischemia, sepsis in
untreated necrotizing infection, cure by early surgery, selection for resistance, the effect of
clindamycin and vancomycin, the harm of frequent cytotoxic antiseptics; `tests/literature.rs` —
the default parameters reproduce the calibration targets from the literature.

## Approach

Regeneration of a whole body cannot be modelled in one go, so we start with a single tissue — skin —
and build the core so that further tissues and scales can be added later.

**A continuum "2.5D" model.** A grid of 96×48 cells of 0.25 mm (seen from above). Each cell is a
column of tissue: epidermis 0.1 mm → dermis 2 mm → subcutaneous fat 6 mm → fascia → muscle.
A cell stores cell densities and signal concentrations, plus the depth of the defect, the deepest
level the damage has reached, and the thickness of dead tissue. Full 3D (voxels through the depth)
will be needed for undermining and tracking under the wound edges — not yet.

### What is modelled

| Process | How |
|---|---|
| Hemostasis | bleeding (heavier in deep wounds) → clot → platelet-derived growth factors |
| Inflammation | debris, bacteria, necrosis, biofilm → signal → neutrophils and M1 macrophages; efferocytosis switches M1→M2 |
| Depth | in shallow wounds the surviving dermis provides epithelial islands from the follicles and leaves almost no scar; deep ones first fill with granulation tissue and only then epithelialize; fat is poorly perfused |
| Oxygen | from the vessels and the wound bed; inflammation dilates the vessels (hyperemia) — if the arteries are able to |
| Necrosis | tissue dies from hypoxia and bacterial toxins; dead tissue blocks healing, feeds bacteria and is out of reach of the immune system and antibiotics; macrophages dissolve it slowly |
| Infection | growth on the open surface and in necrosis (invasive strains — in living tissue too); biofilm hides bacteria; living tissue has a background defence; necrotizing strains release leukocidins and diffusing exotoxins that kill tissue ahead of the bacteria |
| Antiseptics | octenidine, polyhexanide, hypochlorous acid, chlorhexidine, povidone-iodine, hydrogen peroxide: potency, penetration into biofilm and necrosis, inactivation by organic matter (pus, necrosis), residence time on the tissue, separate toxicity to keratinocytes and fibroblasts (the ordering follows in vitro data: hypochlorous acid < PHMB ≈ octenidine < chlorhexidine ≈ iodine < peroxide) |
| Antibiotics | cefazolin (a fast β-lactam), clindamycin (bacteriostatic, penetrates well — into necrosis too — and suppresses toxin synthesis), vancomycin (slow, penetrates poorly, but kills the β-lactam-resistant strain); each has its own pharmacokinetics and a Hill effect on MIC; they reach the tissue with the blood; they can be combined |
| Debridement | surgery — immediate (excision of necrosis and infected tissue, the wound gets larger); maggots — 1–2 weeks (they eat necrosis and the bacteria in it selectively); collagenase — days to weeks; hydrogel — speeds up autolysis, softens dry eschar and speeds up epithelialization in a moist environment |
| Muscle | satellite cells build new fibres: after a small defect ~50% by day 7–8 and ~75–80% by weeks 3–4; the capacity falls exponentially with the volume lost (volumetric muscle loss) — a large defect scars instead; infection, necrosis and ischemia shift the outcome towards fibrosis |
| Fat | barely returns in an adult wound (adipocytes appear only near new follicles): a few per cent over 1–2 months, the rest is fibrous scar and a contour depression where the wound healed |
| Remodeling | collagen maturation; a scar is never stronger than ~84% of healthy skin |
| Body site | the skin of a forearm, a face, a shin or a sole differs in layer thickness, local blood supply, appendage density and resident flora; the site also fixes how large the body region is, which turns the wound into a share of the body surface |

## Body site

`--site` picks where on the body the wound is. The site sets the anatomy of the tissue column
and the local physiology, and the scenario (the patient's disease) is applied on top of it —
the site describes the tissue, the scenario describes the person.

| Site | Epidermis | Dermis | Fat | Blood flow | Appendages | % of body surface |
|---|---|---|---|---|---|---|
| forearm | 0.10 mm | 2.0 mm | 6 mm | ×1.00 | ×1.0 | 3% |
| face | 0.07 mm | 1.2 mm | 2 mm | ×1.35 | ×3.0 | 3% |
| scalp | 0.08 mm | 1.8 mm | 3 mm | ×1.30 | ×8.0 | 3% |
| back | 0.10 mm | 3.5 mm | 8 mm | ×0.90 | ×0.8 | 13% |
| abdomen | 0.08 mm | 2.0 mm | 12 mm | ×0.90 | ×0.7 | 9% |
| shin | 0.09 mm | 1.5 mm | 2 mm | ×0.55 | ×0.6 | 7% |
| sole | 1.00 mm | 2.0 mm | 10 mm | ×0.80 | ×0.35 | 1.5% |
| sacrum | 0.10 mm | 2.0 mm | 4 mm | ×0.60 | ×0.6 | 1.2% |

The forearm is the reference: the model was calibrated on this anatomy, so its numbers are the
previous defaults rather than a separate measurement. The rest are order-of-magnitude figures for
an adult from ultrasound and histology reviews; the body-surface shares follow the rule of nines
against a reference body surface of 1.73 m².

A full-thickness wound (the whole dermis lost, so neither site keeps any appendages) heals at a
speed set by the blood supply and the wound bed:

| Site | Skin thickness | Full-thickness depth | Epithelialization |
|---|---|---|---|
| face | 1.27 mm | 1.77 mm | 27.0 d |
| scalp | 1.88 mm | 2.38 mm | 27.9 d |
| forearm | 2.10 mm | 2.60 mm | 29.5 d |
| abdomen | 2.08 mm | 2.58 mm | 30.0 d |
| back | 3.60 mm | 4.10 mm | 31.2 d |
| sole | 3.00 mm | 3.50 mm | 31.2 d |
| shin | 1.59 mm | 2.09 mm | 33.4 d |
| sacrum | 2.10 mm | 2.60 mm | 34.5 d |

Note that wounds of the same depth in millimetres are **not** comparable across sites: 1 mm takes
three quarters of the thin dermis of a face but only three fifths of a shin's, so on the face the
better perfusion is cancelled out by the wound simply being deeper in relative terms. This falls
out of the model rather than being put in by hand.

## Size of the affected area

The simulated patch is a representative window into the lesion, not the whole of it.
`--area` sets the physical size of that window in millimetres and `--cell` its resolution;
`--lesion` states what share of the body region is affected, which turns into a share of the
whole body surface and drives the systemic severity:

```bash
cargo run --release --bin body-sim -- --site shin --area 60x40 --lesion 12
```

```
Место: передняя поверхность голени | участок 60×40 мм (240×160 клеток по 0.25 мм), область тела 1211 см² (7.0% поверхности)
Поражено 12% этой области = 145 см², это 0.84% поверхности тела; участок — её представительный кусок
```

The diffusion coefficients are calibrated in cells²/h on the 0.25 mm grid, so changing the
resolution converts them to mm²/h; otherwise a coarser grid would silently change the physics.
The conversion is applied once, to the Laplacians, and is exactly 1.0 at the reference resolution,
so previous runs do not change by a single bit. A clinical regression test pins the property: the
same wound closes within 20% of the same day at 0.25 mm and 0.5 mm.

The time step follows the grid: the coefficients grow as 1/cell_mm², so `fit_dt` shrinks `dt` to
stay below the stability limit of the explicit scheme with a margin. The margin matters — exactly
at the limit the checkerboard mode has an amplification factor of −1, so it neither grows nor
decays, the oxygen field rings from cell to cell, dips below the critical level locally and
produces sterile necrosis that is a numerical artefact rather than physiology.

## Calibration

The parameters are calibrated against published data: `cargo run --release --bin calibrate [iterations] [from stage]`
(Nelder–Mead in a logarithmic parameter space with bounds; ~15 min for all stages).
The error is the sum of squared deviations in units of the tolerance. The test `tests/literature.rs`
makes sure the default parameters do not drift away from the targets.

| Observation | Literature | Before calibration | After |
|---|---|---|---|
| neutrophil peak | 1 ± 0.5 days | 1.17 | 0.96 |
| macrophage peak | 2.5 ± 0.75 days | 2.67 | 2.50 |
| neutrophils on day 7 / peak | 0.1 ± 0.1 | 0.22 | 0.10 |
| fibroblast peak | 10 ± 3 days | 8.2 | 10.0 |
| epithelialization of a 0.1 mm abrasion | 8 ± 3 days | 2.4 | 9.0 |
| epithelialization of a 0.3 mm graft donor site | 12 ± 4 days | 2.4 | 9.9 |
| epithelialization of a 1 mm deep dermal wound | 21 ± 6 days | 4.7 | 15.1 |
| epithelialization of a 2.5 mm full-thickness wound | 24 ± 6 days | 20.5 | 28.5 |
| scar strength on day 7 | ~0.03 | 0.10 | 0.07 |
| strength on day 21 | 0.20 ± 0.07 | 0.27 | 0.23 |
| strength on day 42 | 0.40 ± 0.10 | 0.41 | 0.41 |
| strength on day 90 | 0.62 ± 0.08 | 0.60 | 0.58 |

What changed: neutrophils live half as long, macrophages arrive sooner, fibroblasts divide more
slowly; the epithelium grows half as fast, and epithelial islands from the skin appendages thin out
sharply with depth; collagen is deposited more slowly but matures faster. Antiseptic toxicity was
rescaled against the new cell turnover rates.

Sources: healing phases ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC2933384/),
[Healogics](https://www.healogics.com/wound-care-patient-information/understanding-the-stages-of-wound-healing-from-inflammation-to-remodeling/)),
epithelialization of donor sites and abrasions ([systematic review](https://www.sciencedirect.com/science/article/pii/S0305417921000504),
[ScienceDirect Topics](https://www.sciencedirect.com/topics/pharmacology-toxicology-and-pharmaceutical-science/skin-abrasion)),
scar strength (Levenson et al., 1965; [PRS Global Open](https://journals.lww.com/prsgo/fulltext/2013/04000/the_role_of_wound_healing_and_its_everyday.4.aspx)).

### Results (healthy patient, Ø8 mm wound)

| Depth | Epithelialization | Strength on day 90 |
|---|---|---|
| 0.1 mm — epidermis | 9 d | 100% |
| 1 mm — dermis | 15 d | 77% |
| 2.5 mm — full thickness | 28.5 d | 57% |
| 5 mm — subcutaneous fat | 31 d | 57% |
| 9 mm — muscle | 31 d | 59% |

Deep layers: with a 9 mm wound (0.6 mm of muscle lost) — 82% new fibres and 18% fibrosis; with a
12 mm wound (3.6 mm of muscle) fibrosis dominates. Fat: ~18% comes back by day 90, and in place of
the rest there is scar and a contour depression.

### Scenarios without treatment (depth 2.5 mm, 90 days)

| Scenario | Outcome |
|---|---|
| healthy | healed in 28.5 days |
| elderly | 37.7 days |
| diabetic | chronic, healed in 88.8 days |
| infected | biofilm, chronic, healed in 61 days |
| ischemic | dry necrosis of the bed from day 2.5, later a secondary infection, never heals |
| diabetic-foot | necrosis + infection, never heals |
| necrotizing | spreading from day 2.2, sepsis from day 6.1 |

### Antiseptics (infected wound, 90 days; 61 days without treatment)

| Antiseptic | every 24 h | every 12 h |
|---|---|---|
| hypochlorous acid | 44 d (consumed quickly) | **33.5 d** |
| polyhexanide | 35.9 d | 42.5 d |
| octenidine | 39.3 d | 60.4 d |
| povidone-iodine | 41.7 d | never closed |
| chlorhexidine | 62 d — no better than no treatment | never closed |
| hydrogen peroxide | 56.7 d | 56.2 d |

Frequent application of cytotoxic antiseptics holds healing back more than it helps.

### Antibiotics and debridement (90 days)

| Case | Outcome |
|---|---|
| necrotizing + cefazolin alone | sepsis from day 8.5, 79% of the bacteria resistant |
| necrotizing + cefazolin + clindamycin | spreading delayed to day 7.6, but without surgery sepsis by day 25 |
| necrotizing + early surgery (days 0.5, 1.5, 3) + cefazolin + clindamycin | healed in 27.2 days |
| necrotizing + late surgery (days 3, 4) | necrosis recurs, sepsis by day 17.5 |
| diabetic-foot + cefazolin alone | the necrosis stays, 58% of the bacteria resistant |
| diabetic-foot + vancomycin alone | the necrosis stays: without debridement an antibiotic does not help |
| diabetic-foot + surgery + octenidine + an antibiotic | infection and necrosis cleared, the wound stays chronic because of the ischemia |
| ischemic + hydrogel | dry gangrene remains: blood flow was never restored |

## Architecture

```
src/
  lib.rs         the body_sim library — the model core
  params.rs      model constants, the anatomy of a column, 7 clinical scenarios
  grid.rs        Field — a scalar field on the grid, Laplacian (Neumann boundary)
  tissue.rs      Tissue — the state of a tissue patch, inflicting a wound of a given depth
  sim.rs         one model step: explicit Euler, dt = 0.1 h
  therapy.rs     antiseptics, antibiotics (PK/PD), debridement, a log of procedures
  body.rs        body sites: skin anatomy, local physiology, region areas and body-surface shares
  calibration.rs calibration: literature targets, stages, the Nelder–Mead optimizer
  bin/calibrate.rs  runs the calibration
  simulation.rs  the run over time: tissue + treatment + an hourly history of metrics
  report.rs      metrics, phase, clinical condition, ASCII map, CSV
  main.rs        the console binary
  bin/gui/       the graphical binary (macroquad)
    main.rs      application state, layout, tabs, input, camera
    render3d.rs  the 3D scene: surface mesh, section faces, the GLSL skin shader, orbit camera, ray picking
    paint.rs     how the tissue looks: surface colour/wetness/relief, histological section textures
    chart.rs     the timeline chart with phase/condition bands and procedure marks
    ui.rs        buttons, tabs, panels, a Cyrillic font taken from the system
    prof.rs      the per-frame profiler for `--bench`
```

## Known limitations

- The calibration uses literature targets for healthy skin; the parameters of the scenarios
  (diabetes, ischemia, infection) and of the treatments were chosen qualitatively. Two
  epithelialization parameters ran into their bounds (cavity filling at the maximum, appendage
  islands at almost zero): for a Ø8 mm wound the epithelium advances mostly from the edges, so the
  model probably heals large superficial abrasions more slowly than it should.
- Necrotizing infection is cured by a single early operation even without an antibiotic — in
  reality both are needed.
- One resistant strain (to β-lactams); no other resistance mechanisms and no antibiotic side effects.
- Muscle and fat regeneration are described as fractions of the volume restored in a column,
  without individual fibres and lobules.
- No wound contraction, mechanics, edema, pressure (pressure ulcers) or revascularization
  in ischemia.
- No systemic level: sepsis here is an estimate of the infected area, not a model of the organism.
  It fires when the infection takes over the simulated patch, or when the affected body surface
  exceeds 10% — both are area thresholds, not physiology.
- The site figures are order-of-magnitude values for an adult, not a specific patient, and only
  the forearm is tied to the calibration. The effect of appendage density is weaker than in
  reality: `adnexal_rate` sits at its lower calibration bound, so even the scalp's follicles
  cannot reproduce how much faster donor sites heal there.
- The patch is treated as a representative window into the lesion, so `--lesion` scales the
  measured fractions to the real area rather than simulating it. A lesion large enough for its
  own gradients (a burn with a healing edge and a dying centre) is not captured.
- Cells are not agents: no chemotaxis and no directed migration.

## Roadmap

1. ~~**Calibration** against literature curves~~ — done (cells, epithelialization, scar strength);
   next — calibrating the scenarios (diabetes, ischemia) and the treatments against clinical studies.
2. **More treatment**: NPWT (negative pressure), revascularization in ischemia, skin grafting,
   drug side effects.
3. ~~**CI**: run `cargo test` on every push~~ — done (GitHub Actions + coverage).
4. **Systemic level**: temperature, white cell count, lactate, SOFA — sepsis as a state of the patient.
5. **Mechanics**: myofibroblasts, contraction, scar quality.
6. **A hybrid model**: agent-based cells on top of the fields.
7. **3D** and performance: `rayon`, then `wgpu` compute.
