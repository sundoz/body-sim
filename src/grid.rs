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
