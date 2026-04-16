use nalgebra::Matrix6;


/// Generic definition of a material law
/// This trait define the fonction stifness matrix that take into account the deformation, even if the material is linear 
/// It use a reference of the deformation to limit the copying
pub trait MaterialLaw {
    fn stiffness_matrix(&self, deformation: Option<&Matrix6<f64>>) -> Matrix6<f64>;
    fn density(&self)->f64;
}


/// The linear elastic law
/// It define a youngs_modulus et poisson_ratio to calcul the stiffness
pub struct LinearElastic {
    pub youngs_modulus : f64,
    pub poisson_ratio : f64,
    pub density : f64
}

impl MaterialLaw for LinearElastic {

    fn stiffness_matrix(&self, _deformation: Option<&Matrix6<f64>>) -> Matrix6<f64>{
        let nu = self.poisson_ratio;
        let e = self.youngs_modulus;
        e / ((1.0+nu) * (1.0-2.0*nu)) * Matrix6::new( 
            1.0-nu, nu, nu, 0.0, 0.0, 0.0,
            nu, 1.0-nu, nu, 0.0, 0.0, 0.0,
            nu, nu, 1.0-nu, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, (1.0-2.0*nu)/2.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, (1.0-2.0*nu)/2.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, (1.0-2.0*nu)/2.0,
        )
    }

    fn density(&self)->f64{
        self.density
    }


}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_linear_elastic_stifness_matrix_calculation() {
        let material: LinearElastic = LinearElastic { youngs_modulus: 1000.0, poisson_ratio: 0.3, density : 1000. };
        let stiffness = material.stiffness_matrix(None);
        assert!((stiffness[(0,0)] - 1346.1538).abs() < 1e-3);
    }

}