use nalgebra::{DMatrix, DVector};


pub trait DampingModel {
    fn compute(&self, mass: &DVector<f64>, k: &DMatrix<f64>) -> DMatrix<f64>;
}


pub struct RayleighDamping {
    pub alpha: f64,
    pub beta: f64,
}

impl DampingModel for  RayleighDamping {
    fn compute(&self, mass: &DVector<f64>, k: &DMatrix<f64>) -> DMatrix<f64> {
        let n = k.nrows();
        let mut c = DMatrix::zeros(n, n);
        for i in 0..n {
            c[(i, i)] += self.alpha * mass[i];
        }
        c += self.beta * k;
        c
    }
}