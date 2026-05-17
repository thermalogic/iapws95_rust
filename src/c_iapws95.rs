//! IAPWS-95 C FFI Bindings
//! 
//! Provides a C-compatible interface for the IAPWS-95 thermodynamic properties library.
//! Compile with: cargo build --release
//! Output: target/release/libiapws95.so (Linux), target/release/iapws95.dll (Windows)

use crate::iapws95::*;
use crate::iapws95_saturation::calc_saturation_properties;
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
