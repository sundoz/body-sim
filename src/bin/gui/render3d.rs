//! 3D-блок ткани как в анатомическом атласе: рельефная поверхность кожи с раной
//! и грани-срезы через слои. Мир в миллиметрах: x — вдоль ширины, z — вдоль высоты
//! участка, y — вверх (поверхность кожи y = 0, мышца ниже y = -8).

use body_sim::params::Params;
use body_sim::report::View;
use body_sim::tissue::Tissue;
use macroquad::prelude::*;
use rayon::prelude::*;

use crate::paint::{self, CellLook, SectionGeom};

/// Во сколько раз сетка поверхности детальнее сетки модели.
const SUB: usize = 3;
/// Глубина блока, мм.
pub const BLOCK_DEPTH: f32 = 12.0;
const TEX_PX_PER_MM: f32 = 36.0;
const TEX_TOP_MM: f32 = 1.0;
const TEX_H: usize = 480;
const TEX_PX_PER_CELL: usize = 10;

const VERTEX: &str = r#"#version 100
attribute vec3 position;
attribute vec2 texcoord;
attribute vec4 color0;
attribute vec4 normal;

varying lowp vec4 v_color;
varying vec2 v_uv;
varying vec4 v_normal;
varying vec3 v_world;

uniform mat4 Model;
uniform mat4 Projection;

void main() {
    vec4 wp = Model * vec4(position, 1.0);
    v_world = wp.xyz;
    v_color = color0 / 255.0;
    v_uv = texcoord;
    v_normal = normal;
    gl_Position = Projection * wp;
}
"#;

const FRAGMENT: &str = r#"#version 100
#ifdef GL_FRAGMENT_PRECISION_HIGH
precision highp float;
#else
precision mediump float;
#endif

varying lowp vec4 v_color;
varying vec2 v_uv;
varying vec4 v_normal;
varying vec3 v_world;

uniform sampler2D Texture;
uniform vec3 CamPos;
uniform vec3 LightDir;

float hash(vec2 p) {
    p = fract(p * vec2(123.34, 456.21));
    p += dot(p, p + 45.32);
    return fract(p.x * p.y);
}

float vnoise(vec2 p) {
    vec2 i = floor(p);
    vec2 f = fract(p);
    vec2 u = f * f * (3.0 - 2.0 * f);
    float a = hash(i);
    float b = hash(i + vec2(1.0, 0.0));
    float c = hash(i + vec2(0.0, 1.0));
    float d = hash(i + vec2(1.0, 1.0));
    return mix(mix(a, b, u.x), mix(c, d, u.x), u.y);
}

float fbm(vec2 p) {
    float s = 0.0;
    float a = 0.5;
    for (int k = 0; k < 4; k++) {
        s += a * vnoise(p);
        p *= 2.03;
        a *= 0.5;
    }
    return s;
}

// Кожный рельеф: мелкие борозды (кожный рисунок) + поры.
float skin_height(vec2 p) {
    float lines = abs(vnoise(p * vec2(1.3, 3.8)) - 0.5);
    float pores = smoothstep(0.82, 0.95, vnoise(p * 10.0));
    return 0.5 * fbm(p * 3.0) - 0.6 * smoothstep(0.0, 0.06, 0.06 - lines) - 0.9 * pores;
}

vec3 tonemap(vec3 x) {
    return clamp((x * (2.51 * x + 0.03)) / (x * (2.43 * x + 0.59) + 0.14), 0.0, 1.0);
}

void main() {
    vec4 tex = texture2D(Texture, v_uv);
    vec3 albedo = pow(v_color.rgb * tex.rgb, vec3(2.2));
    float skin = v_color.a;
    float wet = v_normal.w;
    vec3 n = normalize(v_normal.xyz);
    vec2 p = v_world.xz;

    if (skin > 0.01) {
        float e = 0.015;
        float h0 = skin_height(p);
        float hx = skin_height(p + vec2(e, 0.0));
        float hz = skin_height(p + vec2(0.0, e));
        vec3 grad = vec3(hx - h0, 0.0, hz - h0) / e;
        n = normalize(n - grad * 0.045 * skin);
        float pore = smoothstep(0.82, 0.95, vnoise(p * 10.0));
        albedo *= 1.0 - skin * (0.4 * pore - 0.14 * (fbm(p * 0.7) - 0.5));
    } else {
        // Мелкая неровность раневой поверхности, чтобы блики не были «пластиковыми».
        float e = 0.02;
        float h0 = fbm(p * 5.0);
        vec3 grad = vec3(fbm((p + vec2(e, 0.0)) * 5.0) - h0, 0.0, fbm((p + vec2(0.0, e)) * 5.0) - h0) / e;
        n = normalize(n - grad * 0.03 * wet * step(0.5, n.y));
    }

    vec3 l = normalize(LightDir);
    vec3 v = normalize(CamPos - v_world);
    vec3 h = normalize(l + v);
    float ndl = dot(n, l);

    // Кожа пропускает свет: мягкий «обёрнутый» диффуз и красноватое рассеяние на терминаторе.
    float wrap = 0.4 * skin;
    float diff = max((ndl + wrap) / (1.0 + wrap), 0.0);
    vec3 sss = albedo * vec3(0.9, 0.25, 0.18) * skin * 0.35 * clamp(1.0 - abs(ndl), 0.0, 1.0);

    float hemi = 0.5 + 0.5 * n.y;
    vec3 ambient = albedo * mix(vec3(0.06, 0.05, 0.06), vec3(0.26, 0.27, 0.30), hemi);
    vec3 l2 = normalize(vec3(-l.x, 0.5, -l.z));
    float fill = max(dot(n, l2), 0.0) * 0.22;

    float shin = mix(24.0, 260.0, wet);
    float fres = 0.04 + 0.96 * pow(1.0 - max(dot(n, v), 0.0), 5.0);
    float spec = pow(max(dot(n, h), 0.0), shin) * mix(0.06, 1.6, wet) * (0.35 + fres) * step(0.0, ndl);

    vec3 col = ambient
        + albedo * diff * vec3(1.0, 0.96, 0.9) * 1.25
        + albedo * fill * vec3(0.8, 0.85, 1.0)
        + sss
        + vec3(spec);
    gl_FragColor = vec4(pow(tonemap(col * 0.82), vec3(1.0 / 2.2)), 1.0);
}
"#;

/// Орбитальная камера вокруг центра блока.
pub struct OrbitCam {
    pub yaw: f32,
    pub pitch: f32,
    pub dist: f32,
}

impl OrbitCam {
    pub fn atlas() -> Self {
        Self { yaw: 0.45, pitch: 0.6, dist: 31.0 }
    }

    pub fn top() -> Self {
        Self { yaw: 0.0, pitch: 1.5, dist: 30.0 }
    }

    fn target() -> Vec3 {
        vec3(0.0, -3.5, 0.0)
    }

    fn position(&self) -> Vec3 {
        let (sy, cy) = self.yaw.sin_cos();
        let (sp, cp) = self.pitch.sin_cos();
        Self::target() + self.dist * vec3(cp * sy, sp, cp * cy)
    }
}

/// Грань блока с текстурой среза.
struct Face {
    tex: Texture2D,
    rgba: Vec<u8>,
    /// Нетронутые слои кожи — копируются в «спокойные» столбцы.
    base: Vec<u8>,
    geom: SectionGeom,
    mesh: Mesh,
    seed: u32,
    /// Для какого состояния сцены текстура нарисована.
    painted: Option<(u64, usize)>,
}

impl Face {
    fn new(p: &Params, cells_along: usize, seed: u32) -> Self {
        let geom =
            SectionGeom { w: cells_along * TEX_PX_PER_CELL, h: TEX_H, top_mm: TEX_TOP_MM, px_per_mm: TEX_PX_PER_MM };
        let mut base = vec![0u8; geom.w * geom.h * 4];
        paint::paint_base(p, cells_along, &geom, &mut base, seed);
        let tex = Texture2D::from_rgba8(geom.w as u16, geom.h as u16, &base);
        tex.set_filter(FilterMode::Linear);
        let rgba = base.clone();
        Self {
            tex,
            rgba,
            base,
            geom,
            mesh: Mesh { vertices: Vec::new(), indices: Vec::new(), texture: None },
            seed,
            painted: None,
        }
    }

    /// Перерисовать срез по линии клеток, если сцена изменилась с прошлого раза.
    fn repaint(&mut self, t: &Tissue, p: &Params, cells: &[usize], key: (u64, usize)) {
        if self.painted == Some(key) {
            return;
        }
        paint::paint_section(t, p, cells, &self.geom, &self.base, &mut self.rgba, self.seed);
        self.tex.update_from_bytes(self.geom.w as u32, self.geom.h as u32, &self.rgba);
        self.painted = Some(key);
    }
}
pub struct Scene3D {
    material: Material,
    w: usize,
    h: usize,
    cell_mm: f32,
    looks: Vec<CellLook>,
    zs: Vec<f32>,
    surface: Mesh,
    shadow: Mesh,
    front: Face,
    back: Face,
    left: Face,
    right: Face,
    pub cam: OrbitCam,
    frame: u64,
    /// Состояние, по которому сцена собрана последний раз: ревизия ткани, слой, срез.
    last_key: Option<(u64, Option<View>, Option<usize>)>,
    /// Время последнего update() по этапам, мс: вид клеток, меш поверхности, текстуры срезов.
    pub timings: [f64; 3],
}

/// Точка на поверхности под курсором.
pub struct Pick {
    /// Непрерывные координаты сетки (центр клетки i — ровно i).
    pub gx: f32,
    pub gy: f32,
}

fn vertex(pos: Vec3, uv: Vec2, col: [u8; 4], normal: Vec4) -> Vertex {
    Vertex { position: pos, uv, color: col, normal }
}

fn to_u8(c: f32) -> u8 {
    (c.clamp(0.0, 1.0) * 255.0) as u8
}

/// Полупрозрачное тёмное пятно под блоком, чтобы он «стоял», а не висел в воздухе.
fn shadow_mesh(hw: f32, hh: f32) -> Mesh {
    let n = 24;
    let (ex, ez) = (hw + 6.0, hh + 6.0);
    let mut vertices = Vec::with_capacity((n + 1) * (n + 1));
    for j in 0..=n {
        for k in 0..=n {
            let x = -ex + 2.0 * ex * k as f32 / n as f32;
            let z = -ez + 2.0 * ez * j as f32 / n as f32;
            // Расстояние за пределы прямоугольника блока.
            let dx = (x.abs() - hw).max(0.0);
            let dz = (z.abs() - hh).max(0.0);
            let d = (dx * dx + dz * dz).sqrt();
            let a = (0.55 * (1.0 - d / 6.0).clamp(0.0, 1.0).powi(2) * 255.0) as u8;
            vertices.push(vertex(vec3(x, -BLOCK_DEPTH - 0.05, z), vec2(0.0, 0.0), [0, 0, 0, a], Vec4::ZERO));
        }
    }
    let mut indices = Vec::with_capacity(n * n * 6);
    for j in 0..n {
        for k in 0..n {
            let a = (j * (n + 1) + k) as u16;
            let c = a + (n + 1) as u16;
            indices.extend_from_slice(&[a, c, a + 1, a + 1, c, c + 1]);
        }
    }
    Mesh { vertices, indices, texture: None }
}

/// Сплайн Кэтмулла–Рома по четырём точкам.
fn catmull(p0: f32, p1: f32, p2: f32, p3: f32, t: f32) -> f32 {
    let t2 = t * t;
    let t3 = t2 * t;
    0.5 * ((2.0 * p1)
        + (-p0 + p2) * t
        + (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t2
        + (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t3)
}

impl Scene3D {
    pub fn new(p: &Params, w: usize, h: usize) -> Self {
        let cell_mm = p.cell_mm;
        let material = load_material(
            ShaderSource::Glsl { vertex: VERTEX, fragment: FRAGMENT },
            MaterialParams {
                pipeline_params: PipelineParams {
                    depth_test: Comparison::LessOrEqual,
                    depth_write: true,
                    ..Default::default()
                },
                uniforms: vec![
                    UniformDesc::new("CamPos", UniformType::Float3),
                    UniformDesc::new("LightDir", UniformType::Float3),
                ],
                textures: vec![],
            },
        )
        .expect("шейдер ткани не скомпилировался");
        Self {
            material,
            w,
            h,
            cell_mm,
            looks: vec![CellLook::default(); w * h],
            zs: vec![0.0; w * h],
            surface: Mesh { vertices: Vec::new(), indices: Vec::new(), texture: None },
            shadow: shadow_mesh(w as f32 * cell_mm / 2.0, h as f32 * cell_mm / 2.0),
            front: Face::new(p, w, 1),
            back: Face::new(p, w, 2),
            left: Face::new(p, h, 3),
            right: Face::new(p, h, 4),
            cam: OrbitCam::atlas(),
            frame: 0,
            last_key: None,
            timings: [0.0; 3],
        }
    }

    fn half_w(&self) -> f32 {
        self.w as f32 * self.cell_mm / 2.0
    }

    fn half_h(&self) -> f32 {
        self.h as f32 * self.cell_mm / 2.0
    }

    /// Глубина видимой поверхности в непрерывных координатах сетки (центр клетки i — i).
    fn z_at(&self, gx: f32, gy: f32) -> f32 {
        let gx = gx.clamp(0.0, (self.w - 1) as f32);
        let gy = gy.clamp(0.0, (self.h - 1) as f32);
        let (ix, iy) = (gx.floor() as i32, gy.floor() as i32);
        let (fx, fy) = (gx - ix as f32, gy - iy as f32);
        let get = |x: i32, y: i32| {
            let x = x.clamp(0, self.w as i32 - 1) as usize;
            let y = y.clamp(0, self.h as i32 - 1) as usize;
            self.zs[y * self.w + x]
        };
        let mut rows = [0.0f32; 4];
        for (k, dy) in (-1..=2).enumerate() {
            rows[k] = catmull(get(ix - 1, iy + dy), get(ix, iy + dy), get(ix + 1, iy + dy), get(ix + 2, iy + dy), fx);
        }
        catmull(rows[0], rows[1], rows[2], rows[3], fy)
    }

    fn look_at(&self, gx: f32, gy: f32) -> CellLook {
        let gx = gx.clamp(0.0, (self.w - 1) as f32);
        let gy = gy.clamp(0.0, (self.h - 1) as f32);
        let x0 = gx.floor() as usize;
        let y0 = gy.floor() as usize;
        let x1 = (x0 + 1).min(self.w - 1);
        let y1 = (y0 + 1).min(self.h - 1);
        let (fx, fy) = (gx - x0 as f32, gy - y0 as f32);
        let l = |x: usize, y: usize| self.looks[y * self.w + x];
        let (a, b, c, d) = (l(x0, y0), l(x1, y0), l(x0, y1), l(x1, y1));
        let lerp = |p: f32, q: f32, r: f32, s: f32| {
            let top = p + (q - p) * fx;
            let bot = r + (s - r) * fx;
            top + (bot - top) * fy
        };
        let albedo = std::array::from_fn(|k| lerp(a.albedo[k], b.albedo[k], c.albedo[k], d.albedo[k]));
        CellLook {
            albedo,
            wet: lerp(a.wet, b.wet, c.wet, d.wet),
            skin: lerp(a.skin, b.skin, c.skin, d.skin),
            ao: lerp(a.ao, b.ao, c.ao, d.ao),
        }
    }

    /// Мировые координаты ↔ непрерывные координаты сетки.
    fn world_x(&self, gx: f32) -> f32 {
        (gx + 0.5) * self.cell_mm - self.half_w()
    }

    fn world_z(&self, gy: f32) -> f32 {
        (gy + 0.5) * self.cell_mm - self.half_h()
    }

    /// Пересобрать сцену по состоянию ткани. `revision` меняется при любом изменении ткани,
    /// `cut` — строка сетки, через которую проходит передний срез. Если ничего не изменилось —
    /// ничего не делаем.
    pub fn update(&mut self, t: &Tissue, p: &Params, revision: u64, view: Option<View>, cut: Option<usize>) {
        self.timings = [0.0; 3];
        self.frame += 1;
        let key = (revision, view, cut);
        let ms = |t: std::time::Instant| t.elapsed().as_secs_f64() * 1000.0;
        let cut_row = cut.unwrap_or(self.h - 1);
        // Координата плоскости переднего среза (в клетках): через центр строки или по краю блока.
        let gy_cut = if cut.is_some() { cut_row as f32 } else { self.h as f32 - 0.5 };

        if self.last_key != Some(key) {
            let t0 = std::time::Instant::now();
            for i in 0..t.len() {
                self.looks[i] = paint::cell_look(t, p, i, view);
                self.zs[i] = paint::surface_z(t, p, i);
            }
            self.timings[0] = ms(t0);
            let t1 = std::time::Instant::now();
            self.build_surface(gy_cut);
            self.build_faces(gy_cut);
            self.timings[1] = ms(t1);
        }

        // Передняя грань — сразу; остальные (там обычно здоровая ткань) — по очереди в разных кадрах,
        // чтобы не рисовать всё в одном кадре.
        let t2 = std::time::Instant::now();
        let face_key = (revision, cut_row);
        // Передняя грань — раз в 3 кадра (ткань меняется плавно), при сдвиге среза — сразу.
        let cut_moved = self.front.painted.map(|k| k.1) != Some(cut_row);
        if cut_moved || self.frame.is_multiple_of(3) {
            let cells: Vec<usize> = (0..self.w).map(|x| cut_row * self.w + x).collect();
            self.front.repaint(t, p, &cells, face_key);
        }
        match self.frame % 30 {
            0 => {
                let cells: Vec<usize> = (0..self.w).collect();
                self.back.repaint(t, p, &cells, (revision, 0));
            }
            10 => {
                let cells: Vec<usize> = (0..self.h).map(|y| y * self.w).collect();
                self.left.repaint(t, p, &cells, (revision, 0));
            }
            20 => {
                let cells: Vec<usize> = (0..self.h).map(|y| y * self.w + self.w - 1).collect();
                self.right.repaint(t, p, &cells, (revision, 0));
            }
            _ => {}
        }
        self.timings[2] = ms(t2);
        self.last_key = Some(key);
    }
    /// Непрерывные координаты сетки вдоль x: от левого края блока (-0.5) до правого (w - 0.5).
    fn xs(&self) -> Vec<f32> {
        let n = self.w * SUB;
        (0..=n).map(|k| k as f32 / SUB as f32 - 0.5).collect()
    }

    fn zs_until(&self, gy_cut: f32) -> Vec<f32> {
        let mut v: Vec<f32> = Vec::new();
        let mut k = 0usize;
        loop {
            let g = k as f32 / SUB as f32 - 0.5;
            if g >= gy_cut - 1e-3 {
                break;
            }
            v.push(g);
            k += 1;
        }
        v.push(gy_cut);
        v
    }

    fn build_surface(&mut self, gy_cut: f32) {
        let xs = self.xs();
        let zs = self.zs_until(gy_cut);
        let (nx, nz) = (xs.len(), zs.len());
        // Высоты и вершины считаются параллельно по строкам меша.
        let this = &*self;
        let mut heights = vec![0.0f32; nx * nz];
        heights.par_chunks_mut(nx).zip(zs.par_iter()).for_each(|(row, &gy)| {
            for (k, &gx) in xs.iter().enumerate() {
                row[k] = -this.z_at(gx, gy);
            }
        });
        // Буфер берём из меша и возвращаем в конце, чтобы не выделять память каждый кадр.
        let mut verts = std::mem::take(&mut self.surface.vertices);
        verts.clear();
        verts.resize(nx * nz, vertex(Vec3::ZERO, Vec2::ZERO, [0; 4], Vec4::ZERO));
        let this = &*self;
        verts.par_chunks_mut(nx).enumerate().for_each(|(j, out)| {
            let gy = zs[j];
            for (k, &gx) in xs.iter().enumerate() {
                let hl = heights[j * nx + k.saturating_sub(1)];
                let hr = heights[j * nx + (k + 1).min(nx - 1)];
                let hd = heights[j.saturating_sub(1) * nx + k];
                let hu = heights[(j + 1).min(nz - 1) * nx + k];
                let dx = (xs[(k + 1).min(nx - 1)] - xs[k.saturating_sub(1)]) * this.cell_mm;
                let dz = (zs[(j + 1).min(nz - 1)] - zs[j.saturating_sub(1)]) * this.cell_mm;
                let n = vec3(-(hr - hl) / dx.max(1e-4), 1.0, -(hu - hd) / dz.max(1e-4)).normalize();
                let look = this.look_at(gx, gy);
                let col = [
                    to_u8(look.albedo[0] * look.ao),
                    to_u8(look.albedo[1] * look.ao),
                    to_u8(look.albedo[2] * look.ao),
                    to_u8(look.skin),
                ];
                let pos = vec3(this.world_x(gx), heights[j * nx + k], this.world_z(gy));
                out[k] = vertex(pos, vec2(0.0, 0.0), col, vec4(n.x, n.y, n.z, look.wet));
            }
        });
        self.surface.vertices = verts;
        let idx = &mut self.surface.indices;
        idx.clear();
        for j in 0..nz - 1 {
            for k in 0..nx - 1 {
                let a = (j * nx + k) as u16;
                let b = a + 1;
                let c = a + nx as u16;
                let d = c + 1;
                idx.extend_from_slice(&[a, c, b, b, c, d]);
            }
        }
    }

    fn build_faces(&mut self, gy_cut: f32) {
        let xs = self.xs();
        let zs = self.zs_until(gy_cut);
        let bottom = -BLOCK_DEPTH;
        let white = [255, 255, 255, 0];

        // Передняя грань (срез) и задняя: вдоль x.
        let mut front = Vec::with_capacity(xs.len() * 2);
        let mut back = Vec::with_capacity(xs.len() * 2);
        let zf = self.world_z(gy_cut);
        let zb = self.world_z(-0.5);
        for &gx in &xs {
            let u = (gx + 0.5) / self.w as f32;
            let top_f = -self.z_at(gx, gy_cut);
            let top_b = -self.z_at(gx, 0.0);
            let x = self.world_x(gx);
            let g = &self.front.geom;
            front.push(vertex(vec3(x, top_f, zf), vec2(u, g.v_of_z(-top_f)), white, vec4(0.0, 0.0, 1.0, 0.0)));
            front.push(vertex(vec3(x, bottom, zf), vec2(u, g.v_of_z(-bottom)), white, vec4(0.0, 0.0, 1.0, 0.0)));
            back.push(vertex(vec3(x, top_b, zb), vec2(u, g.v_of_z(-top_b)), white, vec4(0.0, 0.0, -1.0, 0.0)));
            back.push(vertex(vec3(x, bottom, zb), vec2(u, g.v_of_z(-bottom)), white, vec4(0.0, 0.0, -1.0, 0.0)));
        }
        // Боковые грани: вдоль z до плоскости среза.
        let mut left = Vec::with_capacity(zs.len() * 2);
        let mut right = Vec::with_capacity(zs.len() * 2);
        let xl = self.world_x(-0.5);
        let xr = self.world_x(self.w as f32 - 0.5);
        for &gy in &zs {
            let u = (gy + 0.5) / self.h as f32;
            let z = self.world_z(gy);
            let g = &self.left.geom;
            let top_l = -self.z_at(0.0, gy);
            let top_r = -self.z_at((self.w - 1) as f32, gy);
            left.push(vertex(vec3(xl, top_l, z), vec2(u, g.v_of_z(-top_l)), white, vec4(-1.0, 0.0, 0.0, 0.0)));
            left.push(vertex(vec3(xl, bottom, z), vec2(u, g.v_of_z(-bottom)), white, vec4(-1.0, 0.0, 0.0, 0.0)));
            right.push(vertex(vec3(xr, top_r, z), vec2(u, g.v_of_z(-top_r)), white, vec4(1.0, 0.0, 0.0, 0.0)));
            right.push(vertex(vec3(xr, bottom, z), vec2(u, g.v_of_z(-bottom)), white, vec4(1.0, 0.0, 0.0, 0.0)));
        }
        for (face, verts) in
            [(&mut self.front, front), (&mut self.back, back), (&mut self.left, left), (&mut self.right, right)]
        {
            let n = verts.len() / 2;
            face.mesh.indices.clear();
            for k in 0..n - 1 {
                let a = (2 * k) as u16;
                face.mesh.indices.extend_from_slice(&[a, a + 1, a + 2, a + 2, a + 1, a + 3]);
            }
            face.mesh.vertices = verts;
            face.mesh.texture = Some(face.tex.clone());
        }
    }

    pub fn camera(&self, r: Rect) -> Camera3D {
        let sh = screen_height();
        Camera3D {
            position: self.cam.position(),
            target: OrbitCam::target(),
            up: vec3(0.0, 1.0, 0.0),
            fovy: 32f32.to_radians(),
            aspect: Some(r.w / r.h),
            viewport: Some((r.x as i32, (sh - r.y - r.h) as i32, r.w as i32, r.h as i32)),
            ..Default::default()
        }
    }

    pub fn draw(&self, r: Rect, pick: Option<&Pick>, brush_mm: f32) {
        let cam = self.camera(r);
        set_camera(&cam);
        let pos = self.cam.position();
        self.material.set_uniform("CamPos", [pos.x, pos.y, pos.z]);
        let light = vec3(-0.45, 0.85, 0.55).normalize();
        self.material.set_uniform("LightDir", [light.x, light.y, light.z]);
        // Мягкая контактная тень под блоком.
        draw_mesh(&self.shadow);
        gl_use_material(&self.material);
        draw_mesh(&self.surface);
        for face in [&self.front, &self.back, &self.left, &self.right] {
            draw_mesh(&face.mesh);
        }
        gl_use_default_material();

        // Контур кисти, повторяющий рельеф.
        if let Some(p) = pick {
            let r_cells = brush_mm / self.cell_mm;
            let mut prev: Option<Vec3> = None;
            for k in 0..=48 {
                let a = k as f32 / 48.0 * std::f32::consts::TAU;
                let gx = p.gx + r_cells * a.cos();
                let gy = p.gy + r_cells * a.sin();
                let pt = vec3(self.world_x(gx), -self.z_at(gx, gy) + 0.05, self.world_z(gy));
                if let Some(q) = prev {
                    draw_line_3d(q, pt, Color::new(1.0, 1.0, 1.0, 0.9));
                }
                prev = Some(pt);
            }
        }
        set_default_camera();
    }

    /// Экранные координаты мировой точки (или None, если за камерой).
    pub fn project(&self, r: Rect, p: Vec3) -> Option<Vec2> {
        let m = self.camera(r).matrix();
        let c = m * p.extend(1.0);
        if c.w <= 0.0 {
            return None;
        }
        let ndc = c.truncate() / c.w;
        Some(vec2(r.x + (ndc.x + 1.0) / 2.0 * r.w, r.y + (1.0 - ndc.y) / 2.0 * r.h))
    }

    /// Луч из-под курсора до пересечения с поверхностью (трассировка по карте высот).
    pub fn pick(&self, r: Rect, mouse: Vec2, gy_limit: f32) -> Option<Pick> {
        if !r.contains(mouse) {
            return None;
        }
        let inv = self.camera(r).matrix().inverse();
        let nx = (mouse.x - r.x) / r.w * 2.0 - 1.0;
        let ny = 1.0 - (mouse.y - r.y) / r.h * 2.0;
        let unproject = |z: f32| {
            let v = inv * vec4(nx, ny, z, 1.0);
            v.truncate() / v.w
        };
        let a = unproject(-1.0);
        let b = unproject(1.0);
        let dir = (b - a).normalize();
        let mut pt = a;
        let step = 0.05;
        for _ in 0..4000 {
            pt += dir * step;
            if pt.y < -BLOCK_DEPTH - 1.0 {
                break;
            }
            let gx = (pt.x + self.half_w()) / self.cell_mm - 0.5;
            let gy = (pt.z + self.half_h()) / self.cell_mm - 0.5;
            let inside = gx >= -0.5 && gx <= self.w as f32 - 0.5 && gy >= -0.5 && gy <= gy_limit;
            if inside && pt.y <= -self.z_at(gx, gy) {
                return Some(Pick { gx, gy });
            }
        }
        None
    }

    /// Точка на левом или правом ребре передней грани на глубине `depth_mm` — для подписей слоёв.
    pub fn front_edge_point(&self, right: bool, gy_cut: f32, depth_mm: f32) -> Vec3 {
        let gx = if right { self.w as f32 - 0.5 } else { -0.5 };
        vec3(self.world_x(gx), -depth_mm, self.world_z(gy_cut))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn catmull_passes_through_control_points() {
        assert!((catmull(0.0, 1.0, 3.0, 4.0, 0.0) - 1.0).abs() < 1e-6);
        assert!((catmull(0.0, 1.0, 3.0, 4.0, 1.0) - 3.0).abs() < 1e-6);
        assert!((catmull(2.0, 2.0, 2.0, 2.0, 0.37) - 2.0).abs() < 1e-6, "константа остаётся константой");
        let mid = catmull(0.0, 1.0, 2.0, 3.0, 0.5);
        assert!((mid - 1.5).abs() < 1e-6, "прямая остаётся прямой");
    }

    #[test]
    fn orbit_camera_looks_at_the_block() {
        for cam in [OrbitCam::atlas(), OrbitCam::top()] {
            let d = (cam.position() - OrbitCam::target()).length();
            assert!((d - cam.dist).abs() < 1e-3);
            assert!(cam.position().y > OrbitCam::target().y, "камера над блоком");
        }
    }

    #[test]
    fn shadow_is_darkest_under_the_block() {
        let m = shadow_mesh(12.0, 6.0);
        let center = m.vertices.iter().min_by(|a, b| a.position.length().total_cmp(&b.position.length())).unwrap();
        let corner = m.vertices.iter().max_by(|a, b| a.position.length().total_cmp(&b.position.length())).unwrap();
        assert!(center.color[3] > 100 && corner.color[3] == 0);
        assert!(m.vertices.iter().all(|v| (v.position.y + BLOCK_DEPTH).abs() < 0.1));
    }
}
