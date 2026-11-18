/*
This is the main file of the relativistic hydrocode that I wrote
for my general relativity class.

Author: Brayden JoHantgen
Last Update: 11/18/26
*/

////////////
// Imports
////////////
#![allow(unused_imports)]
use std::time::Instant;
use std::f64::consts::PI;

pub mod physics;
pub mod simulation;
pub mod test_conds;

/////////////////////////////
// Simulation Build and Run
/////////////////////////////
fn main() {
    // Initial Conditions
    let init_setup = test_conds::sod_shock();

    let drive = simulation::Driver{
                cfl: 0.3,
                tfinal: 1.001,
                checkpoint: 0.0125,
                n_cell: 1000,
                n_ghost: 4, 
                r_range: (0.0, 1.0),
    };

    // Print Error Checking
    //let test_prims = simulation::init_prims(&init_setup, &drive);
    //let test_value = physics::cons_from_prims(&test_prims, &init_setup, &drive);
    //println!("{:?}", test_prims);

    // Simulation Prerequisites
    let before = Instant::now();

    let mut t: f64 = 0.0; 
    let mut t_checkpoint = drive.checkpoint;
    let mut step_count: i64 = 0;
    let mut check_count: i64 = 0;

    let mut prims = simulation::init_prims(&init_setup, &drive);
    simulation::boundary_conds(&mut prims, &init_setup, &drive);
    let mut cons: simulation::Conserved;
    let mut dt: f64;

    let _ = simulation::write_checkpoint(&prims, &drive, t, check_count);
    check_count += 1;

    // Running the Simulation 
    while t < drive.tfinal {
        cons = physics::cons_from_prims(&prims, &init_setup, &drive);

        dt = simulation::determine_t_step(&prims, &init_setup, &drive);

        cons = simulation::rk_step(&prims, &cons, &init_setup, &drive, dt);
        prims = physics::prims_from_cons(&cons, &init_setup, &drive);
        simulation::boundary_conds(&mut prims, &init_setup, &drive);

        t += dt;
        step_count += 1;

        if t >= t_checkpoint {
            let _ = simulation::write_checkpoint(&prims, &drive, t, check_count);
            t_checkpoint += drive.checkpoint;
            check_count += 1;
        }
        
        println!("step={} t={} dt={}", step_count, t, dt);
    }

    simulation::run_diagnostics(&drive, before.elapsed().as_millis() as f64, step_count);
}