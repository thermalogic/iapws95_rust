use pyo3::prelude::*;

use crate::iapws95::*;
// ==========================================================================
// Functions for (T,rho) → property calculations - direct computation
// ==========================================================================

#[pyfunction]
fn tr2p(t_c: f64, rho: f64) -> PyResult<f64> {
    let t_k = t_c + 273.15;
    Ok(calc_pressure(t_k, rho))
}

#[pyfunction]
fn tr2u(t_c: f64, rho: f64) -> PyResult<f64> {
    let t_k = t_c + 273.15;
    Ok(calc_internal_energy(t_k, rho))
}

#[pyfunction]
fn tr2h(t_c: f64, rho: f64) -> PyResult<f64> {
    let t_k = t_c + 273.15;
    Ok(calc_enthalpy(t_k, rho))
}

#[pyfunction]
fn tr2s(t_c: f64, rho: f64) -> PyResult<f64> {
    let t_k = t_c + 273.15;
    Ok(calc_entropy(t_k, rho))
}

#[pyfunction]
fn tr2cv(t_c: f64, rho: f64) -> PyResult<f64> {
    let t_k = t_c + 273.15;
    Ok(calc_cv(t_k, rho))
}

#[pyfunction]
fn tr2cp(t_c: f64, rho: f64) -> PyResult<f64> {
    let t_k = t_c + 273.15;
    Ok(calc_cp(t_k, rho))
}

#[pyfunction]
fn tr2w(t_c: f64, rho: f64) -> PyResult<f64> {
    let t_k = t_c + 273.15;
    Ok(calc_speed_of_sound(t_k, rho))
}

// ==========================================================================
// Saturation properties
// ==========================================================================

#[pyfunction]
fn sat_t(t_c: f64) -> PyResult<(f64, f64, f64, f64, f64, f64, f64)> {
    let sat = iapws95_saturation::sat_t(t_c)
        .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyValueError, _>(
            format!("Temperature {} C is out of valid saturation range", t_c)
        ))?;
    Ok((sat.p_sat, sat.rho_l, sat.rho_v, sat.h_l, sat.h_v, sat.s_l, sat.s_v))
}

#[pymodule]
fn iapws95(m: &Bound<'_, PyModule>) -> PyResult<()> {
    
    // (T,rho) → property calculations
    m.add_function(wrap_pyfunction!(tr2p, m)?)?;
    m.add_function(wrap_pyfunction!(tr2u, m)?)?;
    m.add_function(wrap_pyfunction!(tr2h, m)?)?;
    m.add_function(wrap_pyfunction!(tr2s, m)?)?;
    m.add_function(wrap_pyfunction!(tr2cv, m)?)?;
    m.add_function(wrap_pyfunction!(tr2cp, m)?)?;
    m.add_function(wrap_pyfunction!(tr2w, m)?)?;
   
    // Saturation properties
    m.add_function(wrap_pyfunction!(saturation_properties, m)?)?;
    
    Ok(())
}
