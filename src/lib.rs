//! Ядро симулятора: модель ткани, параметры, шаг интегрирования, лечение, метрики.
//! Используется консольным (`body-sim`) и графическим (`body-sim-gui`) бинарниками.

pub mod grid;
pub mod params;
pub mod report;
pub mod sim;
pub mod simulation;
pub mod therapy;
pub mod tissue;
