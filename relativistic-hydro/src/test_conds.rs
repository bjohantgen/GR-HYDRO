/*
This file contains the initial conditions of many hydrodynamic test
problems

Author: Brayden JoHantgen
Last Update: 11/17/26
*/

//////////////
// Importing 
//////////////
use crate::simulation;

/////////////////////////////////////
// Defining Functions of Test Cases
/////////////////////////////////////
// Sod Shock Tube
pub fn sod_shock() -> simulation::InitConds {
    let sod = simulation::InitConds{
        adiabatic_index: 1.4,
        discontinuity: 0.5,
        bh_mass: 0.0,
        softness: 0.0,
        length: 0.0,
        p: (1.0, 0.125),
        rho: (1.0, 0.1),
        v1: (0.0, 0.0), 
        v2: (0.0, 0.0),
        v3: (0.0, 0.0),
    };
    sod
}

// Harder Shock Tube
pub fn hard_shock() -> simulation::InitConds {
    let hard = simulation::InitConds{
        adiabatic_index: 1.4,
        discontinuity: 0.5,
        bh_mass: 0.0,
        softness: 0.0,
        length: 0.0,
        p: (100.0, 1.0),
        rho: (10.0, 1.0),
        v1: (0.0, 0.0), 
        v2: (0.0, 0.0),
        v3: (0.0, 0.0),
    };
    hard
}

// Bondi-Accretion
pub fn bondi_accretion_state_5() -> simulation::InitConds {
    let bondi = simulation::InitConds{
        adiabatic_index: 5.0 / 3.0,
        discontinuity: 0.5,
        bh_mass: 1.0,
        softness: 2.0,
        length: 16.0,
        p: (3.0 / 5.0, 3.0 / 5.0),
        rho: (1.0, 1.0),
        v1: (0.0, 0.0), 
        v2: (0.0, 0.0),
        v3: (0.0, 0.0),
    };
    bondi
}