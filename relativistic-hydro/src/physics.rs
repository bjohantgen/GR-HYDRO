/*
This file contains all the physics equations used in the general
relativistic hydrocode.

Author: Brayden JoHantgen
Last Update: 11/13/26
*/

////////////
// Imports
////////////
use std::f64::consts::PI;
use crate::simulation;

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

//////////////////////////////////////////
// Defining Position Dependent Functions
//////////////////////////////////////////

//////////////////////////////////////////////////////////////
// Conserved Variable and Primitive Variable Transformations
//////////////////////////////////////////////////////////////
