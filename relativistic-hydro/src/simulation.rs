/*
This file contains all the simulation functions and datatypes used 
in the general relativistic hydrocode.

Author: Brayden JoHantgen
Last Update: 11/18/26
*/

////////////
// Imports
////////////
#![allow(unused_imports)]
use std::fs;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::f64::consts::PI;

use crate::physics;

///////////////
// Data Types
///////////////
pub type Cell5 = (f64, f64, f64, f64, f64);
pub type Cell3 = (f64, f64, f64);

pub struct InitConds {
    pub adiabatic_index: f64,
    pub discontinuity: f64,
    pub bh_mass: f64,
    pub softness: f64,
    pub length: f64,
    pub p: (f64, f64),
    pub rho: (f64, f64),
    pub v1: (f64, f64),
    pub v2: (f64, f64),
    pub v3: (f64, f64),
}

pub struct Driver {
    pub cfl: f64,
    pub tfinal: f64,
    pub checkpoint: f64,
    pub n_cell: usize,
    pub n_ghost: usize,
    pub r_range: (f64, f64),
}

#[derive(Debug)]
pub struct Primitives {
    pub p: Vec<f64>, 
    pub rho: Vec<f64>,  
    pub v1: Vec<f64>,   
    pub v2: Vec<f64>,   
    pub v3: Vec<f64>,  
}

#[derive(Debug)]
pub struct Conserved {
    pub rho: Vec<f64>,
    pub m1: Vec<f64>,
    pub m2: Vec<f64>,
    pub m3: Vec<f64>,
    pub e: Vec<f64>,
}

//////////////////////////////////////////
// Defining PLM Reconstruction Functions
//////////////////////////////////////////

/// Input:
///     tup (Cell3): a tuple with three elements
/// Output:
///     max_check (f64): the maximum element of the tuple
/// Description:
///     This function compares the elements of a tuple with negative infinity
///     and stores the variable if it is larger than the previous. The function
///     checks for NaN elements of the tuple.
pub fn tuple_max(tup: Cell3) -> f64 {
    let arr = [tup.0, tup.1, tup.2];
    let mut max_check: f64 = f64::NEG_INFINITY;
    for i in 0..3 {
        debug_assert!(!arr[i].is_nan(), "tuple_max received NaN input");
        if arr[i] > max_check {
            max_check = arr[i];
        }
    }
    max_check
}

/// Input:
///     tup (Cell3): a tuple with three elements
/// Output:
///     min_check (f64): the minimum element of the tuple
/// Description:
///     This function compares the elements of a tuple with infinity
///     and stores the variable if it is smaller than the previous. The 
///     function checks for NaN elements of the tuple.
pub fn tuple_min(tup: Cell3) -> f64 {
    let arr = [tup.0, tup.1, tup.2];
    let mut min_check: f64 = f64::INFINITY;
    for i in 0..3 {
        debug_assert!(!arr[i].is_nan(), "tuple_min received NaN input");
        if arr[i] < min_check {
            min_check = arr[i]
        }
    }
    min_check
}

/// Input:
///     num (f64): any number
/// Output:
///     sign (f64): either 1.0, -1.0, or 0.0 depending on the 
///                 sign of the number
/// Description:
///     If the input number equals 0 then the out put is 0. If the 
///     input number is positive then the number equals 1. If the
///     input number is negative then the output is -1.
pub fn sgn(num: f64) -> f64 {
    let sign: f64;
    if num > 0.0 {
        sign = 1.0;
    } else if num < 0.0 {
        sign = -1.0;
    }
    else {
        sign = 0.0;
    }
    sign
}

/// Input:
///     x (f64): argument 1
///     y (f64): argument 2
///     z (f64): argument 3
/// Output:
///     mm (f64): the value calculated by the minmod funciton
/// Description:
///     Equation: minmod(x,y,z) = 0.25 * |sgn(x) + sgn(y)| * (sgn(x) + sgn(z)) * min(|x|, |y|, |z|)
///     This equation uses the sgn function and the tuple_min function to calculate the minmod function
///     of a minmod slope limiter.
pub fn minmod(x: f64, y: f64, z: f64) -> f64 {
    let mm_1 = (sgn(x) + sgn(y)).abs();
    let mm_2 = sgn(x) + sgn(z);
    let mm_3 = tuple_min((x.abs(), y.abs(), z.abs()));
    let mm = 0.25 * mm_1 * mm_2 * mm_3;
    mm
}

/// Input:
///     ci (f64): desired prim of the ith cell on the grid
///     cip1 (f64): desired prim of the ith+1 cell on the grid
///     cip2 (f64): desired prim of the ith+2 cell on the grid
///     cim1 (f64): desired prim of the ith-1 cell on the grid
///     left (bool): if this is true the left-biased interface value of the left state will be reconstructed.
/// Output:
///     a (f64): the reconstructed prim value
/// Description:
///     This function uses the minmod function to reconstruct the primitive variables
///     of multiple cells using the piecewise linear method in order to get higher order
///     spacial calculations.
pub fn plm_reconstruction(ci: f64, cip1: f64, cip2: f64, cim1: f64, left: bool) -> f64 {
    let a: f64;
    if left == true {
        let num_1 = 1.5 * (ci - cim1);
        let num_2 = 0.5 * (cip1 - cim1);
        let num_3 = 1.5 * (cip1 - ci);
        a = ci + 0.5 * minmod(num_1, num_2, num_3);
    } else {
        let num_1 = 1.5 * (cip1 - ci);
        let num_2 = 0.5 * (cip2 - ci);
        let num_3 = 1.5 * (cip2 - cip1);
        a = cip1 - 0.5 * minmod(num_1, num_2, num_3);
    }
    a
}

///////////////////
// Grid Functions
///////////////////

/// Input:
///     init_conds (&InitConds): the initial conditions of the simulation
///     drive (&Driver): the simulation parameters
/// Output:
///     prims (Primitives): the primitive variables for each cell in the grid
/// Description:
///     This funciton uses drive to initialize the grid and then uses init_conds
///     to populate the grid with the initial primitive varibles. The grid is 
///     uniformly spaced with two values on either side of a discontinuity.
pub fn init_prims(init_conds: &InitConds, drive: &Driver) -> Primitives {
    let n_cells = drive.n_cell + drive.n_ghost;

    let mut init_p = vec![0.0; n_cells];
    let mut init_rho = vec![0.0; n_cells];
    let mut init_v1 = vec![0.0; n_cells];
    let mut init_v2 = vec![0.0; n_cells];
    let mut init_v3 = vec![0.0; n_cells];

    for i in 0..n_cells {
        let cell_check = ((n_cells as f64) * init_conds.discontinuity) as usize;
        if i < cell_check {
            init_p[i] = init_conds.p.0;
            init_rho[i] = init_conds.rho.0;
            init_v1[i] = init_conds.v1.0;
            init_v2[i] = init_conds.v2.0;
            init_v3[i] = init_conds.v3.0;
        } else {
            init_p[i] = init_conds.p.1;
            init_rho[i] = init_conds.rho.1;
            init_v1[i] = init_conds.v1.1;
            init_v2[i] = init_conds.v2.1;
            init_v3[i] = init_conds.v3.1;
        }
    }

    let prims = Primitives{
        p: init_p,
        rho: init_rho,
        v1: init_v1,
        v2: init_v2,
        v3: init_v3,
    };

    prims
}

/// Input:
///     init_conds (&InitConds): the initial conditions of the simulation
///     drive (&Driver): the simulation parameters
/// Output:
///     delta (f64): the width of the cell in units of R_b
/// Description:
///     This function calculates the width of a cell using the length
///     and cell number given by the Driver construct.
pub fn cell_width(drive: &Driver) -> f64 {
    let delta = (drive.r_range.1 - drive.r_range.0) / drive.n_cell as f64;
    delta
}

/// Input:
///     prims (&Primitives): the grid containing all the primitive variables for each cell
///     i (usize): the index of the desired cell
/// Output:
///     cell (Cell5): the primitive variables of the desired cell
/// Description:
///     This function takes the index of a cell to return the primitive variables of that cell on the grid
pub fn cell_loc(prims: &Primitives, i: usize) -> Cell5 {
    let cell = (prims.p[i], prims.rho[i], prims.v1[i], prims.v2[i], prims.v3[i]);
    cell
}

/// Input:
///     prims (&mut Primitives): the grid containing all the primitive variables for each cell
///     i_org (usize): the index of the original cell
///     i_new (usize): the index of the changed cell
/// Output:
///     the updated prims 
/// Description:
///     This function takes the primitive variables on the grid and replaces a
///     desired index with a different chosen index.
pub fn copy_cell(prims: &mut Primitives, i_org: usize, i_new: usize) {
    prims.p[i_new] = prims.p[i_org];
    prims.rho[i_new] = prims.rho[i_org];
    prims.v1[i_new] = prims.v1[i_org];
    prims.v2[i_new] = prims.v2[i_org];
    prims.v3[i_new] = prims.v3[i_org];
}

/// Input:
///     prims (&mut Primitives): the grid containing all the primitive variables for each cell
///     init_conds (&InitConds): the initial conditions of the simulation
/// Output:
///     the updated prims with adjusted boundary conditions
/// Description:
///     This function takes the prims grid and adds the boundary conditiond to the ghost cells on either
///     side of the grid. The left bounday is that of a black hole and requires the gas to fall through.
///     The right boundary pulls from an infinite gas reservoir with the same parameters as the initial conditions.
pub fn boundary_conds(prims: &mut Primitives, init_conds: &InitConds, drive: &Driver) {
    let n_cells = prims.p.len();
    let inner_active = drive.n_ghost / 2;
    let outer_active = n_cells - inner_active;

    // left boundary
    for i in 0..inner_active {
        prims.p[i] = init_conds.p.0;
        prims.rho[i] = init_conds.rho.0;
        prims.v1[i] = init_conds.v1.0;
        prims.v2[i] = init_conds.v2.0;
        prims.v3[i] = init_conds.v3.0;
    }
    // right boundary
    for i in outer_active..n_cells {
        prims.p[i] = init_conds.p.1;
        prims.rho[i] = init_conds.rho.1;
        prims.v1[i] = init_conds.v1.1;
        prims.v2[i] = init_conds.v2.1;
        prims.v3[i] = init_conds.v3.1;
    }
}


//////////////////////////////////
// Defining HLL Solver Functions
//////////////////////////////////

/// Input:
///     prim_1 (Cell5): the primitive variables in position 1
///     prim_2 (Cell5): the primitive variables in position 2
///     prim_3 (Cell5): the primitive variables in position 3
///     prim_4 (Cell5): the primitive variables in position 4
///     a_ind (f64): the adiabatic index of the gas
///     i_num1 (usize): the index of cell 2
///     i_num2 (usize): the index of cell 3
/// Output:
///     hll (Cell5): the hll plm reconstructed hll flux between cells 2 and 3
/// Description:
///     This function uses the plm_reconstruction, sound_speed, prim_to_con, flux, and tuple max functions
///     to calculate the hll flux.
pub fn hll(prim_1: Cell5, prim_2: Cell5, prim_3: Cell5, prim_4: Cell5, a_ind: f64, i_num1: usize, i_num2: usize) -> Cell5 {
    let p_l = plm_reconstruction(prim_2.0, prim_3.0, prim_4.0, prim_1.0, true);
    let rho_l = plm_reconstruction(prim_2.1, prim_3.1, prim_4.1, prim_1.1, true);
    let v1_l = plm_reconstruction(prim_2.2, prim_3.2, prim_4.2, prim_1.2, true);
    let v2_l = plm_reconstruction(prim_2.3, prim_3.3, prim_4.3, prim_1.3, true);
    let v3_l = plm_reconstruction(prim_2.4, prim_3.4, prim_4.4, prim_1.4, true);

    let p_r = plm_reconstruction(prim_2.0, prim_3.0, prim_4.0, prim_1.0, false);
    let rho_r = plm_reconstruction(prim_2.1, prim_3.1, prim_4.1, prim_1.1, false);
    let v1_r = plm_reconstruction(prim_2.2, prim_3.2, prim_4.2, prim_1.2, false);
    let v2_r = plm_reconstruction(prim_2.3, prim_3.3, prim_4.3, prim_1.3, false);
    let v3_r = plm_reconstruction(prim_2.4, prim_3.4, prim_4.4, prim_1.4, false);

    let prim_l = (p_l, rho_l, v1_l, v2_l, v3_l);
    let prim_r = (p_r, rho_r, v1_r, v2_r, v3_r);

    let cf_l = physics::sound_speed(prim_l, a_ind);
    let plus_l = v1_l + cf_l;
    let minus_l = v1_l - cf_l;
    let u_l = physics::prim_to_con(prim_l, a_ind);
    let f_l = physics::flux(prim_l, a_ind);

    let cf_r = physics::sound_speed(prim_r, a_ind);
    let plus_r = v1_r + cf_r;
    let minus_r = v1_r - cf_r;
    let u_r = physics::prim_to_con(prim_r, a_ind);
    let f_r = physics::flux(prim_r, a_ind);

    let a_plus = tuple_max((0.0, plus_l, plus_r));
    let a_minus = tuple_max((0.0, -minus_l, -minus_r));

    let hll_0 = ((a_plus * f_l.0) + (a_minus * f_r.0) - (a_plus * a_minus * (u_r.0 - u_l.0))) / (a_minus + a_plus);
    let hll_1 = ((a_plus * f_l.1) + (a_minus * f_r.1) - (a_plus * a_minus * (u_r.1 - u_l.1))) / (a_minus + a_plus);
    let hll_2 = ((a_plus * f_l.2) + (a_minus * f_r.2) - (a_plus * a_minus * (u_r.2 - u_l.2))) / (a_minus + a_plus);
    let hll_3 = ((a_plus * f_l.3) + (a_minus * f_r.3) - (a_plus * a_minus * (u_r.3 - u_l.3))) / (a_minus + a_plus);
    let hll_4 = ((a_plus * f_l.4) + (a_minus * f_r.4) - (a_plus * a_minus * (u_r.4 - u_l.4))) / (a_minus + a_plus);

    let hll = (hll_0, hll_1, hll_2, hll_3, hll_4);
    hll
}

/// Input:
///     prims (&Primitives): the grid containing all the primitive variables for each cell
///     init_conds (&InitConds): the initial conditions of the simulation
///     drive (&Driver): the simulation parameters
/// Output:
///     go (Vec<Cell5>): the vector of godonov fluxes
/// Description:
///     This function calculates all of the hll fluxes on the grid.
pub fn godonov(prims: &Primitives, init_conds: &InitConds, drive: &Driver) -> Vec<Cell5> {
    let n = prims.p.len() as usize;
    let mut go = vec![(0.0, 0.0, 0.0, 0.0, 0.0); n-1];

    for i in 1..(n-2) {
        go[i] = hll(cell_loc(prims, i-1), cell_loc(prims, i), cell_loc(prims, i+1), cell_loc(prims, i+2), init_conds.adiabatic_index, i, i+1);
    }
    go
}

/// Input:
///     prims (&Primitives): the grid containing all the primitive variables for each cell
///     init_conds (&InitConds): the initial conditions of the simulation
///     drive (&Driver): the simulation parameters
/// Output:
///     l_vec (Conserved): this is the l function used to update the conserved variables.
/// Description:
///     This function uses the cell_width, godonov, and source functions to calculate the l_functions
///     to update the conserved variables.
pub fn l_function(prims: &Primitives, init_conds: &InitConds, drive: &Driver) -> Conserved {
    let n = prims.p.len() as usize;

    let dr = cell_width(&drive);

    let mut l_rho = vec![0.0; n];
    let mut l_m1 = vec![0.0; n];
    let mut l_m2 = vec![0.0; n];
    let mut l_m3 = vec![0.0; n];
    let mut l_e = vec![0.0; n];

    let f_vec = godonov(prims, init_conds, drive);

    let left_bound = drive.n_ghost / 2;
    let right_bound = n - (drive.n_ghost / 2);

    for i in 0..left_bound {
        l_rho[i] = 0.0;
        l_m1[i] = 0.0;
        l_m2[i] = 0.0;
        l_m3[i] = 0.0;
        l_e[i] = 0.0;
    }

    for i in right_bound..n {
        l_rho[i] = 0.0;
        l_m1[i] = 0.0;
        l_m2[i] = 0.0;
        l_m3[i] = 0.0;
        l_e[i] = 0.0;
    }

    for i in left_bound..right_bound {
        l_rho[i] =  - (f_vec[i].0 - f_vec[i-1].0) / dr;
        l_m1[i]  =  - (f_vec[i].1 - f_vec[i-1].1) / dr;
        l_m2[i]  =  - (f_vec[i].2 - f_vec[i-1].2) / dr;
        l_m3[i]  =  - (f_vec[i].3 - f_vec[i-1].3) / dr;
        l_e[i]   =  - (f_vec[i].4 - f_vec[i-1].4) / dr;
    }

    let l_vec = Conserved{
        rho: l_rho,
        m1: l_m1,
        m2: l_m2,
        m3: l_m3,
        e: l_e,
    };

    l_vec
}

////////////////////////////////
// Defining Timestep Functions
////////////////////////////////

/// Input:
///     prim_l (Cell5): the primitive variables for the left cell
///     prim_r (Cell5): the primitive variables for the right cell
///     a_ind (f64): the adiabatic index
///     dr (f64): the width of a cell
/// Output:
///     dt (f64): the timestep allowed for those cells
/// Description:
///     This function uses the sound_speed and tuple_max functions to determine
///     the timestep satisfied by the Courant-Friedrich-Levy condition.
pub fn compute_time_step(prim_l: Cell5, prim_r: Cell5, a_ind: f64, dr: f64) -> f64 {
    let cf_l = physics::sound_speed(prim_l, a_ind);

    let plus_l = cf_l + prim_l.2;
    let minus_l = prim_l.2 - cf_l;

    let cf_r = physics::sound_speed(prim_r, a_ind); 

    let plus_r = cf_r + prim_r.2; 
    let minus_r = prim_r.2 - cf_r;

    let a_plus = tuple_max((0.0, plus_l, plus_r));
    let a_minus = tuple_max((0.0, -minus_l, -minus_r));

    let mut dt: f64 = 0.0;
    
    if a_minus > a_plus {
        dt += dr / a_minus;
    } else {
        dt += dr / a_plus;
    }    
    dt
}

/// Input:
///     prims (&Primitives): the grid containing all the primitive variables for each cell
///     init_conds (&InitConds): the initial conditions of the simulation
///     drive (&Driver): the simulation parameters
/// Output:
///     dt (f64): the time step
/// Description:
///     This function uses the cell_width, compute_time_step, distance, and gravity_time_step functions to 
///     calculate the time step satisfying the Courant-Friedrich-Levy condition. This function checks for 
///     non-finite values from non physical primitive variables.
pub fn determine_t_step(prims: &Primitives, init_conds: &InitConds, drive: &Driver) -> f64 {
    let mut dt = f64::INFINITY;
    let n = prims.p.len();

    let dr = cell_width(drive);

        for i in 1..(n - 2) {
            let prim_fill_1 = cell_loc(prims, i);
            let prim_fill_2 = cell_loc(prims, i+1);

            let dt_check = compute_time_step(prim_fill_1, prim_fill_2, init_conds.adiabatic_index, dr);

            assert!(
                dt_check.is_finite(),
                "non-finite dt at interface {}: p=({}, {}) rho=({}, {})",
                i, prims.p[i], prims.p[i + 1], prims.rho[i], prims.rho[i + 1]
            );

            dt = dt.min(dt_check);
        }

    dt = drive.cfl * dt;
    dt
}

/// Input:
///     prims (&Primitives): the grid containing all the primitive variables for each cell
///     init_conds (&InitConds): the initial conditions of the simulation
///     drive (&Driver): the simulation parameters
///     dt (f64): the size of the timestep
/// Output:
///     cons2 (Conserved): the time stepped conserved variable.
/// Description:
///     This function uses the l_function, prims_from_cons, and boundary_conds functions to advance the
///     conserved variables through and RK3 timestep.
pub fn rk_step(prims: &Primitives, cons: &Conserved, init_conds: &InitConds, drive: &Driver, dt: f64) -> Conserved {
    let n = prims.p.len() as usize;

    let l0 = l_function(prims, init_conds, drive);
    
    let mut con0_rho = vec![0.0; n];
    let mut con0_m1 = vec![0.0; n];
    let mut con0_m2 = vec![0.0; n];
    let mut con0_m3 = vec![0.0; n];
    let mut con0_e = vec![0.0; n];

    for i in 0..n {
        con0_rho[i] = cons.rho[i] + dt * l0.rho[i];
        con0_m1[i] = cons.m1[i] + dt * l0.m1[i];
        con0_m2[i] = cons.m2[i] + dt * l0.m2[i];
        con0_m3[i] = cons.m3[i] + dt * l0.m3[i];
        con0_e[i] = cons.e[i] + dt * l0.e[i];
    }

    let cons0 = Conserved{
        rho: con0_rho,
        m1: con0_m1, 
        m2: con0_m2,
        m3: con0_m3,
        e: con0_e,
    };

    let mut prims_1 = physics::prims_from_cons(&cons0, init_conds, drive);
    boundary_conds(&mut prims_1, init_conds, drive);

    let l1 = l_function(&prims_1, init_conds, drive);

    let mut con1_rho = vec![0.0; n];
    let mut con1_m1 = vec![0.0; n];
    let mut con1_m2 = vec![0.0; n];
    let mut con1_m3 = vec![0.0; n];
    let mut con1_e = vec![0.0; n];

    let dt25 = 0.25 * dt;

    for i in 0..n{
        con1_rho[i] = 0.75 * cons.rho[i] + 0.25 * cons0.rho[i] + dt25 * l1.rho[i];
        con1_m1[i] = 0.75 * cons.m1[i] + 0.25 * cons0.m1[i] + dt25 * l1.m1[i];
        con1_m2[i] = 0.75 * cons.m2[i] + 0.25 * cons0.m2[i] + dt25 * l1.m2[i];
        con1_m3[i] = 0.75 * cons.m3[i] + 0.25 * cons0.m3[i] + dt25 * l1.m3[i];
        con1_e[i] = 0.75 * cons.e[i] + 0.25 * cons0.e[i] + dt25 * l1.e[i];
    }

    let cons1 = Conserved{
        rho: con1_rho,
        m1: con1_m1, 
        m2: con1_m2,
        m3: con1_m3,
        e: con1_e,
    };

    let mut prims_2 = physics::prims_from_cons(&cons1, init_conds, drive);
    boundary_conds(&mut prims_2, init_conds, drive);

    let l2 = l_function(&prims_2, init_conds, drive);

    let mut con2_rho = vec![0.0; n];
    let mut con2_m1 = vec![0.0; n];
    let mut con2_m2 = vec![0.0; n];
    let mut con2_m3 = vec![0.0; n];
    let mut con2_e = vec![0.0; n];

    let dt66 = (2.0 / 3.0) * dt;
    let c33 = 1.0 / 3.0;
    let c66 = 2.0 / 3.0;

    for i in 0..n{
        con2_rho[i] = c33 * cons.rho[i] + c66 * cons1.rho[i] + dt66 * l2.rho[i];
        con2_m1[i] = c33 * cons.m1[i] + c66 * cons1.m1[i] + dt66 * l2.m1[i];
        con2_m2[i] = c33 * cons.m2[i] + c66 * cons1.m2[i] + dt66 * l2.m2[i];
        con2_m3[i] = c33 * cons.m3[i] + c66 * cons1.m3[i] + dt66 * l2.m3[i];
        con2_e[i] = c33 * cons.e[i] + c66 * cons1.e[i] + dt66 * l2.e[i];
    }

    let cons2 = Conserved{
        rho: con2_rho,
        m1: con2_m1, 
        m2: con2_m2,
        m3: con2_m3,
        e: con2_e,
    };

    cons2
}

///////////////////////////////
// Checkpoint and Diagnostics
///////////////////////////////

/// Input:
///     prims (&Primitives): the grid containing all the primitive variables for each cell
///     drive (&Driver): the simulation parameters
///     t (f64): the time of the checkpoint
///     check_count (i64): the chcekpoint number
/// Output:
///     A file containing the recorded values
/// Description:
///     This function writes the number of cells, grid range, and all primitive variables into txt files for plotting.
pub fn write_checkpoint(prims: &Primitives, drive: &Driver, t: f64, check_count: i64) -> Result<(), Box<dyn std::error::Error>> {
    let file_num = check_count.to_string();
    let file_type = ".txt".to_string();
    let file_name = format!("{}{}", file_num, file_type);
    let file_path = "time_step_files/".to_string() + &file_name;
    let output_file_path = Path::new(&file_path);
    
    let t_fill = format!("{} {}", "t:".to_string(), &(t.to_string()+&" ".to_string()));
    
    let size = prims.p.len();

    let n_fill = format!("{} {}", "n_cell:".to_string(), &(drive.n_cell.to_string()+&" ".to_string()));
    let n_ghost_fill = format!("{} {}", "n_ghost_cell:".to_string(), &(drive.n_ghost.to_string()+&" ".to_string()));
    let left_bound_fill = format!("{} {}", "left_bound:".to_string(), &(drive.r_range.0.to_string()+&" ".to_string()));
    let right_bound_fill = format!("{} {}", "right_bound:".to_string(), &(drive.r_range.1.to_string()+&" ".to_string()));
    
    let mut p_string = "p: ".to_string();
    let mut rho_string = "rho: ".to_string();
    let mut v1_string = "v1: ".to_string();
    let mut v2_string = "v2: ".to_string();
    let mut v3_string = "v3: ".to_string();
        
    for i in 0..size {
        let p_fill = prims.p[i].to_string() + &" ".to_string();
        p_string.push_str(&p_fill);
    
        let rho_fill = prims.rho[i].to_string() + &" ".to_string();
        rho_string.push_str(&rho_fill);
    
        let v1_fill = prims.v1[i].to_string() + &" ".to_string();
        v1_string.push_str(&v1_fill);
    
        let v2_fill = prims.v2[i].to_string() + &" ".to_string();
        v2_string.push_str(&v2_fill);
    
        let v3_fill = prims.v3[i].to_string() + &" ".to_string();
        v3_string.push_str(&v3_fill);
    }
     
    let file = fs::OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(output_file_path)?;
    
    let mut writer = BufWriter::new(file);
    
    writeln!(writer, "{}", n_fill)?;
    writeln!(writer, "{}", n_ghost_fill)?;
    writeln!(writer, "{}", left_bound_fill)?;
    writeln!(writer, "{}", right_bound_fill)?;

    writeln!(writer, "{}", t_fill)?;

    writeln!(writer, "{}", p_string)?;
    writeln!(writer, "{}", rho_string)?;
    writeln!(writer, "{}", v1_string)?;
    writeln!(writer, "{}", v2_string)?;
    writeln!(writer, "{}", v3_string)?;

    writer.flush()?;
    
    Ok(())
}

/// Input:
///     drive (&Driver): the simulation parameters
///     t_run (f64): the runtime of the simulation in milliseconds
///     step_count (i64): the number of time steps 
/// Output:
///     print statements giving the number of time steps, the runtime
///     in seconds, the performance of zones per second
/// Description:
///     This function gives the diagnostics of the simulation
pub fn run_diagnostics(drive: &Driver, t_run: f64, step_count: i64) {
    let n = drive.n_cell + drive.n_ghost;
    
    let runtime = t_run / 1000.0;
    let performance = ((step_count as f64) * (n as f64)) / runtime;

    println!("The number of timesteps: {:?}", step_count);
    println!("The runtime in seconds: {:?}", runtime);
    println!("The performance in zones per second: {:.2?}", performance);
}