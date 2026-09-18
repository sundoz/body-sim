/// Скалярное поле на регулярной 2D-сетке (вид сверху на участок кожи,
/// значения усреднены по толщине дермы).
#[derive(Clone, Debug)]
pub struct Field {
    pub w: usize,
    pub h: usize,
    pub data: Vec<f32>,
}

impl Field {
    pub fn new(w: usize, h: usize, value: f32) -> Self {
        Self { w, h, data: vec![value; w * h] }
    }

    /// 5-точечный лапласиан (в единицах клеток сетки) с условием
    /// непротекания (Неймана) на краях: за краем — такая же ткань.
    pub fn laplacian_into(&self, out: &mut [f32]) {
        let (w, h) = (self.w, self.h);
        let d = &self.data;
        for y in 0..h {
            let up = if y == 0 { 0 } else { y - 1 };
            let dn = if y + 1 == h { y } else { y + 1 };
            for x in 0..w {
                let lf = if x == 0 { 0 } else { x - 1 };
                let rt = if x + 1 == w { x } else { x + 1 };
                let c = d[y * w + x];
                out[y * w + x] =
                    d[up * w + x] + d[dn * w + x] + d[y * w + lf] + d[y * w + rt] - 4.0 * c;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn laplacian_of_constant_field_is_zero() {
        let f = Field::new(8, 5, 3.0);
        let mut out = vec![1.0; 40];
        f.laplacian_into(&mut out);
        assert!(out.iter().all(|v| v.abs() < 1e-6));
    }

    #[test]
    fn laplacian_conserves_mass_with_neumann_boundary() {
        let (w, h) = (7, 6);
        let mut f = Field::new(w, h, 0.0);
        f.data[0] = 1.0; // угол
        f.data[3 * w + 3] = 2.0; // внутренняя клетка
        let mut out = vec![0.0; w * h];
        f.laplacian_into(&mut out);
        let sum: f32 = out.iter().sum();
        assert!(sum.abs() < 1e-5, "диффузия не должна создавать или терять вещество: {sum}");
        assert!((out[3 * w + 3] + 8.0).abs() < 1e-6);
        assert!((out[3 * w + 4] - 2.0).abs() < 1e-6);
        // В углу два соседа «за краем» равны самой клетке.
        assert!((out[0] + 2.0).abs() < 1e-6);
    }
}