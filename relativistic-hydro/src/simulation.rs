/*
This file contains all the simulation functions and datatypes used 
in the general relativistic hydrocode.

Author: Brayden JoHantgen
Last Update: 11/13/26
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

//////////////////////////////////
// Defining HLL Solver Functions
//////////////////////////////////


////////////////////////////////
// Defining Timestep Functions
////////////////////////////////