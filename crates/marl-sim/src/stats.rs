use marl_cell::cell::CellState;
use marl_field::field::Field;
use marl_field::light::LightField;
use std::time::Instant;

pub fn print_stats(
    tick: u32,
    cells: &[CellState],
    field: &Field,
    _light: &LightField,
    div: u64,
    death: u64,
    start: &Instant,
) {
    if cells.is_empty() {
        println!("t={:>5} | EXTINCT | div={} death={}", tick, div, death);
        return;
    }

    let grid = field.grid();

    let n = cells.len() as f32;
    let avg_energy: f32 = cells.iter().map(|c| c.internal[0]).sum::<f32>() / n;
    let avg_enzyme: f32 = cells.iter().map(|c| c.internal[5]).sum::<f32>() / n;
    let active_rxns: f32 = cells
        .iter()
        .map(|c| {
            c.ruleset
                .reactions
                .iter()
                .filter(|r| r.v_max.abs() > 1e-9)
                .count() as f32
        })
        .sum::<f32>()
        / n;

    // Count cells per z-third (surface / middle / deep)
    let z_third = (grid.z / 3) as u16;
    let (mut n_top, mut n_mid, mut n_bot) = (0u32, 0u32, 0u32);
    for c in cells {
        if c.pos[2] < z_third {
            n_top += 1;
        } else if c.pos[2] < z_third * 2 {
            n_mid += 1;
        } else {
            n_bot += 1;
        }
    }

    // Sample chemistry at center column
    let cx = grid.x / 2;
    let cy = grid.y / 2;
    let ox_top = field.get(cx, cy, 0, 1);
    let ox_mid = field.get(cx, cy, grid.z / 2, 1);
    let red_bot = field.get(cx, cy, grid.z - 1, 2);
    let org_mid = field.get(cx, cy, grid.z / 2, 4);

    let elapsed = start.elapsed().as_secs_f32();
    let tps = if elapsed > 0.0 {
        (tick + 1) as f32 / elapsed
    } else {
        0.0
    };

    println!(
        "t={:>5} | pop={:>6} (top:{:>5} mid:{:>5} bot:{:>5}) | E={:.2} enz={:.3} rxn={:.1} | ox={:.2}/{:.2} red={:.2} org={:.2} | {:.1} t/s",
        tick,
        cells.len(),
        n_top,
        n_mid,
        n_bot,
        avg_energy,
        avg_enzyme,
        active_rxns,
        ox_top,
        ox_mid,
        red_bot,
        org_mid,
        tps,
    );
}

pub fn print_z_profile(cells: &[CellState], field: &Field, light: &LightField) {
    let grid = field.grid();
    println!(
        "{:>3} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6}",
        "z", "cells", "light", "oxidnt", "reduct", "carbon", "organic"
    );
    let cx = grid.x / 2;
    let cy = grid.y / 2;
    for z in 0..grid.z {
        let n = cells.iter().filter(|c| c.pos[2] == z as u16).count();
        let l = light.get(cx, cy, z);
        let ox = field.get(cx, cy, z, 1);
        let re = field.get(cx, cy, z, 2);
        let ca = field.get(cx, cy, z, 3);
        let og = field.get(cx, cy, z, 4);
        println!(
            "{:>3} {:>6} {:>6.3} {:>6.3} {:>6.3} {:>6.3} {:>6.3}",
            z, n, l, ox, re, ca, og
        );
    }
}
