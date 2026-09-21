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
    ///
    /// Считается по строкам: срезы строк дают компилятору известную длину, поэтому из
    /// внутреннего цикла уходят и проверки границ, и ветвления на краевые столбцы —
    /// края обрабатываются отдельно. На профиле это был самый дорогой цикл в программе.
    pub fn laplacian_into(&self, out: &mut [f32]) {
        let (w, h) = (self.w, self.h);
        assert_eq!(out.len(), w * h, "буфер лапласиана не по размеру поля");
        if w == 0 || h == 0 {
            return;
        }
        for y in 0..h {
            let up = if y == 0 { 0 } else { y - 1 };
            let dn = if y + 1 == h { y } else { y + 1 };
            let row = &self.data[y * w..y * w + w];
            let row_up = &self.data[up * w..up * w + w];
            let row_dn = &self.data[dn * w..dn * w + w];
            let out_row = &mut out[y * w..y * w + w];
            if w == 1 {
                // Единственный столбец: оба горизонтальных соседа — сама клетка.
                out_row[0] = row_up[0] + row_dn[0] - 2.0 * row[0];
                continue;
            }
            // Края: сосед за границей участка равен самой клетке.
            out_row[0] = row_up[0] + row_dn[0] + row[0] + row[1] - 4.0 * row[0];
            out_row[w - 1] = row_up[w - 1] + row_dn[w - 1] + row[w - 2] + row[w - 1] - 4.0 * row[w - 1];
            // Середина строки — через срезы одинаковой длины: так компилятор видит,
            // что они не пересекаются и что итераций ровно w-2, и векторизует цикл.
            // При индексации со смещениями он этого не доказывал и считал по одному числу.
            let mid = w - 2;
            let left = &row[..mid];
            let centre = &row[1..1 + mid];
            let right = &row[2..2 + mid];
            let up = &row_up[1..1 + mid];
            let dn = &row_dn[1..1 + mid];
            let out_mid = &mut out_row[1..1 + mid];
            for k in 0..mid {
                out_mid[k] = up[k] + dn[k] + left[k] + right[k] - 4.0 * centre[k];
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
    fn laplacian_handles_degenerate_grids() {
        // Края считаются отдельно от середины, поэтому узкие сетки — особый случай.
        for (w, h) in [(1, 1), (1, 5), (5, 1), (2, 2)] {
            let mut f = Field::new(w, h, 2.0);
            let mut out = vec![9.0; w * h];
            f.laplacian_into(&mut out);
            assert!(out.iter().all(|v| v.abs() < 1e-6), "{w}×{h}: константа даёт нулевой лапласиан");
            f.data[0] = 3.0;
            f.laplacian_into(&mut out);
            let sum: f32 = out.iter().sum();
            assert!(sum.abs() < 1e-5, "{w}×{h}: вещество не сохраняется ({sum})");
        }
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
