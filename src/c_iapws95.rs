//! IAPWS-95 C FFI Bindings
//! 
//! Provides a C-compatible interface for the IAPWS-95 thermodynamic properties library.
//! Compile with: cargo build --release
//! Output: target/release/libiapws95.so (Linux), target/release/iapws95.dll (Windows)

use crate::iapws95::*;
use crate::iapws95_saturation::calc_saturation_properties;
use crate::iapws95_pt::solve_density;
/// C-compatible saturation properties structure
#[repr(C)]
pub struct CIAPWS95SatProps {
    pub p_sat: f64,
    pub rho_l: f64,
    pub rho_v: f64,
    pub h_l: f64,
    pub h_v: f64,
    pub s_l: f64,
    pub s_v: f64,
}

// ==========================================================================
// Functions for (p,T) → property calculations using numerical inversion
// ==========================================================================

#[no_mangle]
pub extern "C" fn iapws95_pt2h(p: f64, t_c: f64) -> f64 {
    let t_k = t_c + 273.15;
    
    if p <= 0.0 || t_k <= 0.0 {
        return -1.0;
    }

    match solve_density(p, t_k) {
        Some(rho) => calc_enthalpy(t_k, rho),
        None => -1.0,
    }
}

#[no_mangle]
pub extern "C" fn iapws95_pt2s(p: f64, t_c: f64) -> f64 {
    let t_k = t_c + 273.15;
    
    if p <= 0.0 || t_k <= 0.0 {
        return -1.0;
    }

    match solve_density(p, t_k) {
        Some(rho) => calc_entropy(t_k, rho),
        None => -1.0,
    }
}

// ==========================================================================
// Functions for (T,rho) → property calculations - direct computation
// ==========================================================================

#[no_mangle]
pub extern "C" fn iapws95_tr2p(t_c: f64, rho: f64) -> f64 {
    tr2p(t_c, rho)
}

#[no_mangle]
pub extern "C" fn iapws95_tr2u(t_c: f64, rho: f64) -> f64 {
    tr2u(t_c, rho)
}

#[no_mangle]
pub extern "C" fn iapws95_tr2h(t_c: f64, rho: f64) -> f64 {
    tr2h(t_c, rho)
}

#[no_mangle]
pub extern "C" fn iapws95_tr2s(t_c: f64, rho: f64) -> f64 {
    tr2s(t_c, rho)
}

#[no_mangle]
pub extern "C" fn iapws95_tr2cv(t_c: f64, rho: f64) -> f64 {
    tr2cv(t_c, rho)
}

#[no_mangle]
pub extern "C" fn iapws95_tr2cp(t_c: f64, rho: f64) -> f64 {
    tr2cp(t_c, rho)
}

#[no_mangle]
pub extern "C" fn iapws95_tr2w(t_c: f64, rho: f64) -> f64 {
    tr2w(t_c, rho)
}

// ==========================================================================
// Functions for (T,x) → property calculations in two-phase region
// ==========================================================================

#[no_mangle]
pub extern "C" fn iapws95_tx2h(t_c: f64, x: f64) -> f64 {
    let t_k = t_c + 273.15;
    
    if t_k < 273.16 || t_k > IAPWS95_TCRIT {
        return -1.0;
    }

    if let Some(sat) = calc_saturation_properties(t_k) {
        let h = sat.h_l + x * (sat.h_v - sat.h_l);
        return h;
    }

    -1.0
}

#[no_mangle]
pub extern "C" fn iapws95_tx2s(t_c: f64, x: f64) -> f64 {
    let t_k = t_c + 273.15;
    
    if t_k < 273.16 || t_k > IAPWS95_TCRIT {
        return -1.0;
    }

    if let Some(sat) = calc_saturation_properties(t_k) {
        let s = sat.s_l + x * (sat.s_v - sat.s_l);
        return s;
    }

    -1.0
}

// ==========================================================================
// Saturation properties
// ==========================================================================

#[no_mangle]
pub extern "C" fn iapws95_saturation_properties(t_c: f64, props: *mut CIAPWS95SatProps) -> i32 {
    if props.is_null() {
        return -1;
    }

    let t_k = t_c + 273.15;
    
    if t_k < 273.16 || t_k > IAPWS95_TCRIT {
        return -1;
    }

    if let Some(sat) = calc_saturation_properties(t_k) {
        unsafe {
            (*props).p_sat = sat.p_sat;
            (*props).rho_l = sat.rho_l;
            (*props).rho_v = sat.rho_v;
            (*props).h_l = sat.h_l;
            (*props).h_v = sat.h_v;
            (*props).s_l = sat.s_l;
            (*props).s_v = sat.s_v;
        }
        0
    } else {
        -1
    }
}

#[no_mangle]
pub extern "C" fn iapws95_version() -> *const std::os::raw::c_char {
    static VERSION: &[u8] = b"0.1.0\0";
    VERSION.as_ptr() as *const std::os::raw::c_char
}
