/*
This file contains all the physics equations used in the general
relativistic hydrocode. The characteristic values are:

Author: Brayden JoHantgen
Last Update: 11/16/26
*/

////////////
// Imports
////////////
use std::f64::consts::PI;

use crate::simulation::Cell5;
use crate::simulation::Primitives;
use crate::simulation::InitConds;
use crate::simulation::Driver;
use crate::simulation::Conserved;

/////////////////////////////////////
// Defining Position Free Functions
/////////////////////////////////////

/// Input:
///     prim (Cell5): the primitive variables of a cell
///     a_ind (f64): the adiabatic index of the gas
/// Output:
///     c (f64): the sound speed of the gas
/// Description:
///     Equation: c^2 = gamma * P / rho
///     This function calculates the sound speed of the gas.
///     It does no error checking.
pub fn sound_speed(prim: Cell5, a_ind: f64) -> f64 {
    let c = (a_ind * prim.0 / prim.1).sqrt();
    c
}

/// Input:
///     prim (Cell5): the primitive variables of a cell
///     a_ind (f64): the adiabatic index of the gas
/// Output:
///     epsilon (f64): the internal energy density
/// Description:
///     Equation: epsilon = P / (rho * (gamma - 1))
///     This equation calcualtes the internal energy density and
///     does not check for errors.
pub fn internal_energy(prim: Cell5, a_ind: f64) -> f64 {
    let epsilon = prim.0 / (prim.1 * (a_ind - 1.0));
    epsilon
}

/// Input:
///     prim (Cell5): the primitive variables of a cell
///     a_ind (f64): the adiabatic index of the gas
/// Output:
///     energy (f64): the total energy of the cell
/// Description:
///     Equation: E = rho * (v^2/2 + epsilon)
///     This functin uses the internal_energy function to calculate 
///     the total energy of the cell. This function does no error checking.
pub fn total_energy(prim: Cell5, a_ind: f64) -> f64 {
    let epsilon = internal_energy(prim, a_ind);

    let v2 = prim.2 * prim.2 + prim.3 * prim.3 + prim.4 * prim.4;
    let energy = prim.1 * (0.5 * v2 + epsilon);
    energy
}

//////////////////////////////////////////
// Defining Position Dependent Functions
//////////////////////////////////////////

//////////////////////////////////////////////////////////////
// Conserved Variable and Primitive Variable Transformations
//////////////////////////////////////////////////////////////

/// Input:
///     prim (Cell5): the primitive variables of a cell
///     a_ind (f64): the adiabatic index of the gas
/// Output:
///     con (Cell5): the conserved variables of the cell
/// Description:
///     This function uses the total energy function to construct the conserved variables
///     using the primitive variables.
pub fn prim_to_con(prim: Cell5, a_ind: f64) -> Cell5 {
    let energy = total_energy(prim, a_ind);

    let u0 = prim.1;
    let u1 = prim.1 * prim.2;
    let u2 = prim.1 * prim.3;
    let u3 = prim.1 * prim.4;
    let u4 = energy;

    let con = (u0, u1, u2, u3, u4);
    con
}

/// Input:
///     prims (&Primitives): the grid containing all the primitive variables for each cell
///     init_conds (&InitConds): the initial conditions of the simulation
///     drive (&Driver): the simulation parameters
/// Output:
///     cons (Conserved): the conserved variables on the grid.
/// Description:
///     This function uses the cell_width and prim_to_con functions to construct the conserved variables
///     on the simulation grid.
pub fn cons_from_prims(prims: &Primitives, init_conds: &InitConds, drive: &Driver) -> Conserved {
    let n = prims.p.len() as usize;

    let mut cons_rho = vec![0.0; n];
    let mut cons_m1 = vec![0.0; n];
    let mut cons_m2 = vec![0.0; n];
    let mut cons_m3 = vec![0.0; n];
    let mut cons_e = vec![0.0; n];

    for i in 0..n {
        let prim_fill = (prims.p[i], prims.rho[i], prims.v1[i], prims.v2[i], prims.v3[i]);
        let con_result = prim_to_con(prim_fill, init_conds.adiabatic_index);

        cons_rho[i] = con_result.0;
        cons_m1[i] = con_result.1;
        cons_m2[i] = con_result.2;
        cons_m3[i] = con_result.3;
        cons_e[i] = con_result.4;
    }

    let cons = Conserved{
        rho: cons_rho,
        m1: cons_m1,
        m2: cons_m2,
        m3: cons_m3,
        e: cons_e,
    };

    cons
}

/// Input:
///     con (Cell5): the conserved variables of a cell
///     a_ind (f64): the adiabatic index of the gas
/// Output:
///     prim (Cell5): the primitive variables of a cell
/// Description:
///     This function transforms the conserved variables to the primitive variables. 
///     There is a pressure and density floor implemented to prevent these 
///     values from being negative.
pub fn con_to_prim(con: Cell5, a_ind: f64) -> Cell5 {
    let rho = con.0.max(1e-8);

    let v1 = con.1 / rho;
    let v2 = con.2 / rho;
    let v3 = con.3 / rho;

    let v_square = v1 * v1 + v2 * v2 + v3 * v3;

    let p_check = (a_ind - 1.0) * (con.4 - 0.5 * rho * v_square);
    let p = p_check.max(1e-8);

    let prim = (p, rho, v1, v2, v3);
    prim
}

/// Input:
///     cons (&Conserved): the conserved variables on the grid.
///     init_conds (&InitConds): the initial conditions of the simulation
///     drive (&Driver): the simulation parameters
/// Output:
///     prims (Primitives): the grid containing all the primitive variables for each cell
/// Description:
///     This function uses the cell width function and the con_to_prim function to transform the
///     conserved grid to the primitive grid.
pub fn prims_from_cons(cons: &Conserved, init_conds: &InitConds, drive: &Driver) -> Primitives {
    let n = cons.rho.len() as usize;

    let mut prims_p = vec![0.0; n];
    let mut prims_rho = vec![0.0; n];
    let mut prims_v1 = vec![0.0; n];
    let mut prims_v2 = vec![0.0; n];
    let mut prims_v3 = vec![0.0; n];

    for i in 0..n {
        let con_fill = (cons.rho[i], cons.m1[i], cons.m2[i], cons.m3[i], cons.e[i]);
        let prim_result = con_to_prim(con_fill, init_conds.adiabatic_index);

        prims_p[i] = prim_result.0;
        prims_rho[i] = prim_result.1;
        prims_v1[i] = prim_result.2;
        prims_v2[i] = prim_result.3;
        prims_v3[i] = prim_result.4;
    }

    let prims = Primitives{
        p: prims_p,
        rho: prims_rho,
        v1: prims_v1,
        v2: prims_v2,
        v3: prims_v3,
    };
    
    prims
}

///////////////////////
// Fluxes and Sources
///////////////////////

/// Input:
///     prim (Cell5): the primitive variables of a cell
///     a_ind (f64): the adiabatic index of the gas
/// Output:
///     flux (Cell5): the flux of a cell
/// Description:
///     This function uses the total_energy function to calculate the flux.
///     There is no error checking.
pub fn flux(prim: Cell5, a_ind: f64) -> Cell5 {
    let e = total_energy(prim, a_ind);

    let f0 = prim.1 * prim.2;
    let f1 = prim.1 * prim.2 * prim.2 + prim.0;
    let f2 = prim.1 * prim.2 * prim.3;
    let f3 = prim.1 * prim.2 * prim.4;
    let f4 = prim.2 * (e + prim.0);

    let flux = (f0, f1, f2, f3, f4);
    flux
}