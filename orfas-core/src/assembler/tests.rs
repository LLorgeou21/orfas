// UTF-8
// assembler/tests.rs — assembler test suite.

use super::*;
use nalgebra::{DMatrix, DVector, Vector3};
use crate::element::Tet4;
use crate::material::{SaintVenantKirchhoff, SimulationContext};
use crate::mesh::Tet4Mesh;

fn unit_mesh() -> Tet4Mesh {
    Tet4Mesh::generate(2, 2, 2, 1.0, 1.0, 1.0)
}

fn svk() -> SaintVenantKirchhoff {
    SaintVenantKirchhoff { youngs_modulus: 1e6, poisson_ratio: 0.3, density: 1000.0 }
}

#[test]
fn test_tangent_symmetric_at_zero() {
    let mesh      = unit_mesh();
    let mat       = svk();
    let assembler = Assembler::<Tet4>::new(&mesh);
    let u         = DVector::zeros(3 * mesh.nodes.len());
    let k         = assembler.assemble_tangent(&mesh, &mat, &u, &SimulationContext::isotropic_static(mesh.elements.len()));
    let diff      = (&k - k.transpose()).abs().max();
    assert!(diff < 1e-10, "K must be symmetric, diff = {:.2e}", diff);
}

#[test]
fn test_internal_forces_zero_at_rest() {
    let mesh      = unit_mesh();
    let mat       = svk();
    let assembler = Assembler::<Tet4>::new(&mesh);
    let u         = DVector::zeros(3 * mesh.nodes.len());
    let f_int     = assembler.assemble_internal_forces(&mesh, &mat, &u, &SimulationContext::isotropic_static(mesh.elements.len()));
    assert!(f_int.norm() < 1e-10, "f_int must be zero at rest, norm = {:.2e}", f_int.norm());
}

#[test]
fn test_internal_forces_consistent_with_tangent() {
    let mesh      = unit_mesh();
    let mat       = svk();
    let assembler = Assembler::<Tet4>::new(&mesh);
    let n         = 3 * mesh.nodes.len();
    let u_zero    = DVector::zeros(n);
    let k         = assembler.assemble_tangent(&mesh, &mat, &u_zero, &SimulationContext::isotropic_static(mesh.elements.len()));

    let mut u_small = DVector::zeros(n);
    u_small[3] = 1e-6;
    u_small[4] = 1e-6;

    let f_int  = assembler.assemble_internal_forces(&mesh, &mat, &u_small, &SimulationContext::isotropic_static(mesh.elements.len()));
    let f_lin  = &k * &u_small;
    let error  = (&f_int - &f_lin).norm() / f_lin.norm();
    assert!(error < 1e-4, "f_int must linearize to K*u for small u, error = {:.2e}", error);
}

#[test]
fn test_mass_assembly() {
    let mesh      = unit_mesh();
    let mat       = svk();
    let assembler = Assembler::<Tet4>::new(&mesh);
    let mass      = assembler.assemble_mass(&mesh, &mat);
    let total     = mass.sum() / 3.0;
    assert!((total - 1000.0).abs() / 1000.0 < 1e-10);
}

#[test]
fn test_assemble_method_comparaison() {
    let mesh      = unit_mesh();
    let mat       = svk();
    let assembler = Assembler::<Tet4>::new(&mesh);
    let u_zero    = DVector::zeros(3 * mesh.nodes.len());
    let k_dense   = assembler.assemble_tangent(&mesh, &mat, &u_zero, &SimulationContext::isotropic_static(mesh.elements.len()));
    let k_sparse  = assembler.assemble_tangent_sparse(&mesh, &mat, &u_zero, &SimulationContext::isotropic_static(mesh.elements.len()));
    let diff      = (k_dense - DMatrix::from(&k_sparse)).abs().max();
    assert!(diff < 1e-9, "diff = {:.2e}", diff);
}

#[test]
fn test_assemble_parrallel_method_comparaison() {
    let mesh      = unit_mesh();
    let mat       = svk();
    let assembler = Assembler::<Tet4>::new(&mesh);
    let u_zero    = DVector::zeros(3 * mesh.nodes.len());
    let k_par     = assembler.assemble_tangent_sparse_parallel(&mesh, &mat, &u_zero, &SimulationContext::isotropic_static(mesh.elements.len()));
    let k_seq     = assembler.assemble_tangent_sparse(&mesh, &mat, &u_zero, &SimulationContext::isotropic_static(mesh.elements.len()));
    let diff      = (DMatrix::from(&k_seq) - DMatrix::from(&k_par)).abs().max();
    assert!(diff < 1e-9, "diff = {:.2e}", diff);
}

#[test]
fn test_assemble_parallel_speedup() {
    let mesh      = Tet4Mesh::generate(30, 30, 30, 1.0, 1.0, 1.0);
    let mat       = svk();
    let assembler = Assembler::<Tet4>::new(&mesh);
    let u_zero    = DVector::zeros(3 * mesh.nodes.len());

    let t0    = std::time::Instant::now();
    let _     = assembler.assemble_tangent_sparse(&mesh, &mat, &u_zero, &SimulationContext::isotropic_static(mesh.elements.len()));
    let t_seq = t0.elapsed();

    let t0    = std::time::Instant::now();
    let _     = assembler.assemble_tangent_sparse_parallel(&mesh, &mat, &u_zero, &SimulationContext::isotropic_static(mesh.elements.len()));
    let t_par = t0.elapsed();

    println!("sequential: {:.2?}  parallel: {:.2?}  speedup: {:.2}x",
        t_seq, t_par, t_seq.as_secs_f64() / t_par.as_secs_f64());
}