use marl_config::{GRID_X, GRID_Y, GRID_Z, S_EXT, SimulationConfig};
use marl_field::field::Field;
use rand::Rng;
use std::collections::{HashMap, HashSet};

/// The 6 face-neighbor offsets in 3D (±x, ±y, ±z).
/// Shared by all neighbor-scanning functions.
const FACE_OFFSETS: [(i16, i16, i16); 6] = [
    (1, 0, 0),
    (-1, 0, 0),
    (0, 1, 0),
    (0, -1, 0),
    (0, 0, 1),
    (0, 0, -1),
];

fn collect_empty_neighbors(
    pos: [u16; 3],
    cell_map: &HashMap<[u16; 3], usize>,
    reserved: Option<&HashSet<[u16; 3]>>,
) -> ([[u16; 3]; 6], usize) {
    let mut neighbors = [[0u16; 3]; 6];
    let mut count = 0;
    for &(dx, dy, dz) in &FACE_OFFSETS {
        let nx = pos[0] as i16 + dx;
        let ny = pos[1] as i16 + dy;
        let nz = pos[2] as i16 + dz;
        if nx >= 0
            && nx < GRID_X as i16
            && ny >= 0
            && ny < GRID_Y as i16
            && nz >= 0
            && nz < GRID_Z as i16
        {
            let p = [nx as u16, ny as u16, nz as u16];
            if !cell_map.contains_key(&p) && reserved.is_none_or(|r| !r.contains(&p)) {
                neighbors[count] = p;
                count += 1;
            }
        }
    }
    (neighbors, count)
}

/// Collect the positions of all in-bounds, empty face-neighbors of a voxel.
pub fn empty_neighbors(pos: [u16; 3], cell_map: &HashMap<[u16; 3], usize>) -> Vec<[u16; 3]> {
    let (neighbors, count) = collect_empty_neighbors(pos, cell_map, None);
    neighbors[..count].to_vec()
}

/// Read the average chemical environment available to a cell from
/// surrounding extracellular space.
///
/// In a real microbial column, cells access dissolved chemicals from the
/// liquid medium around them — not from inside their own body. A cell
/// surrounded by other cells has limited access to nutrients: only empty
/// (liquid-filled) neighboring voxels contribute.
///
/// Returns the mean concentration across all empty face-neighbors. If
/// the cell is completely enclosed by other cells, returns all zeros —
/// the cell is cut off from external resources and will starve.
pub fn read_neighbor_environment(
    pos: [u16; 3],
    field: &Field,
    cell_map: &HashMap<[u16; 3], usize>,
) -> [f32; S_EXT] {
    let (neighbors, count) = collect_empty_neighbors(pos, cell_map, None);
    if count == 0 {
        return [0.0; S_EXT];
    }
    let mut sum = [0.0f32; S_EXT];
    for npos in &neighbors[..count] {
        let voxel = field.read_voxel(npos[0] as usize, npos[1] as usize, npos[2] as usize);
        for s in 0..S_EXT {
            sum[s] += voxel[s];
        }
    }
    let n = count as f32;
    for value in &mut sum {
        *value /= n;
    }
    sum
}

/// Distribute a cell's secretion/consumption deltas to surrounding empty voxels.
///
/// In the real world, metabolic byproducts are released into the liquid
/// medium around the cell, and consumed substrates are depleted from that
/// same medium. The deltas are split equally among all empty face-neighbors.
///
/// If no empty neighbors exist, deltas are lost — the cell cannot exchange
/// chemicals with a fully packed environment (waste heat / trapped products).
pub fn apply_deltas_to_neighbors(
    pos: [u16; 3],
    field: &mut Field,
    cell_map: &HashMap<[u16; 3], usize>,
    deltas: &[f32; S_EXT],
) {
    let (neighbors, count) = collect_empty_neighbors(pos, cell_map, None);
    if count == 0 {
        return; // enclosed cell — deltas are lost
    }
    for s in 0..S_EXT {
        let delta = deltas[s];
        if delta == 0.0 {
            continue;
        }

        if delta > 0.0 {
            let mut split_deltas = [0.0f32; S_EXT];
            split_deltas[s] = delta / count as f32;
            for npos in &neighbors[..count] {
                field.apply_deltas(
                    npos[0] as usize,
                    npos[1] as usize,
                    npos[2] as usize,
                    &split_deltas,
                );
            }
            continue;
        }

        // Consumption is distributed in proportion to available mass so
        // heterogeneous neighbor concentrations do not silently clamp away
        // requested uptake in low-concentration voxels.
        let mut available = [0.0f32; 6];
        let mut total_available = 0.0;
        for (i, npos) in neighbors[..count].iter().enumerate() {
            let value = field.get(npos[0] as usize, npos[1] as usize, npos[2] as usize, s);
            available[i] = value.max(0.0);
            total_available += available[i];
        }
        if total_available <= f32::EPSILON {
            continue;
        }

        let requested = (-delta).min(total_available);
        for (i, npos) in neighbors[..count].iter().enumerate() {
            if available[i] <= 0.0 {
                continue;
            }
            let mut split_deltas = [0.0f32; S_EXT];
            split_deltas[s] = -requested * (available[i] / total_available);
            field.apply_deltas(
                npos[0] as usize,
                npos[1] as usize,
                npos[2] as usize,
                &split_deltas,
            );
        }
    }
}

/// Find an empty voxel `division_neighbor_distance` steps away along a face axis.
/// Daughters bud off and drift — they don't accrete like plant tissue.
/// Falls back to adjacent placement if no long-range positions are available.
pub fn find_empty_neighbor(
    pos: [u16; 3],
    occupied: &HashMap<[u16; 3], usize>,
    rng: &mut impl Rng,
    sim: &SimulationConfig,
) -> Option<[u16; 3]> {
    find_empty_neighbor_avoiding(pos, occupied, None, rng, sim)
}

pub fn find_empty_neighbor_avoiding(
    pos: [u16; 3],
    occupied: &HashMap<[u16; 3], usize>,
    reserved: Option<&HashSet<[u16; 3]>>,
    rng: &mut impl Rng,
    sim: &SimulationConfig,
) -> Option<[u16; 3]> {
    let dist = sim.division_neighbor_distance as i16;
    // Try distance-first (gap between parent and daughter)
    let mut candidates = [[0u16; 3]; 6];
    let mut count = 0;
    for &(dx, dy, dz) in &FACE_OFFSETS {
        let nx = pos[0] as i16 + dx * dist;
        let ny = pos[1] as i16 + dy * dist;
        let nz = pos[2] as i16 + dz * dist;
        if nx >= 0
            && nx < GRID_X as i16
            && ny >= 0
            && ny < GRID_Y as i16
            && nz >= 0
            && nz < GRID_Z as i16
        {
            let npos = [nx as u16, ny as u16, nz as u16];
            if !occupied.contains_key(&npos) && reserved.is_none_or(|r| !r.contains(&npos)) {
                candidates[count] = npos;
                count += 1;
            }
        }
    }
    if count > 0 {
        let idx = rng.random_range(0..count);
        return Some(candidates[idx]);
    }

    // Fallback: adjacent placement if surrounded at distance
    let (adjacent, adjacent_count) = collect_empty_neighbors(pos, occupied, reserved);
    if adjacent_count == 0 {
        None
    } else {
        let idx = rng.random_range(0..adjacent_count);
        Some(adjacent[idx])
    }
}

#[allow(dead_code)] // TODO: used by HGT when re-enabled
pub fn find_cell_neighbor(pos: [u16; 3], cell_map: &HashMap<[u16; 3], usize>) -> Option<usize> {
    let offsets: [(i16, i16, i16); 6] = [
        (1, 0, 0),
        (-1, 0, 0),
        (0, 1, 0),
        (0, -1, 0),
        (0, 0, 1),
        (0, 0, -1),
    ];
    for &(dx, dy, dz) in &offsets {
        let nx = pos[0] as i16 + dx;
        let ny = pos[1] as i16 + dy;
        let nz = pos[2] as i16 + dz;
        if nx >= 0
            && nx < GRID_X as i16
            && ny >= 0
            && ny < GRID_Y as i16
            && nz >= 0
            && nz < GRID_Z as i16
        {
            let npos = [nx as u16, ny as u16, nz as u16];
            if let Some(&idx) = cell_map.get(&npos) {
                return Some(idx);
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn consumption_is_weighted_by_available_neighbor_mass() {
        let mut field = Field::new();
        let cell_map = HashMap::from([([10, 10, 10], 0)]);
        field.set(11, 10, 10, 1, 0.0);
        field.set(9, 10, 10, 1, 1.0);

        let mut deltas = [0.0f32; S_EXT];
        deltas[1] = -0.5;
        apply_deltas_to_neighbors([10, 10, 10], &mut field, &cell_map, &deltas);

        assert_eq!(field.get(11, 10, 10, 1), 0.0);
        assert!((field.get(9, 10, 10, 1) - 0.5).abs() < 1e-6);
    }

    #[test]
    fn empty_neighbor_scan_uses_face_neighbors_only() {
        let cell_map = HashMap::from([([1, 1, 1], 0), ([2, 1, 1], 1)]);
        let neighbors = empty_neighbors([1, 1, 1], &cell_map);

        assert_eq!(neighbors.len(), 5);
        assert!(!neighbors.contains(&[2, 1, 1]));
        assert!(!neighbors.contains(&[2, 2, 1]));
    }

    #[test]
    fn division_neighbor_search_avoids_reserved_birth_positions() {
        let sim = SimulationConfig {
            division_neighbor_distance: 1,
            ..SimulationConfig::default()
        };
        let cell_map = HashMap::from([
            ([1, 1, 1], 0),
            ([0, 1, 1], 1),
            ([1, 0, 1], 2),
            ([1, 2, 1], 3),
            ([1, 1, 0], 4),
            ([1, 1, 2], 5),
        ]);
        let reserved = HashSet::from([[2, 1, 1]]);
        let mut rng = rand::rng();

        let found =
            find_empty_neighbor_avoiding([1, 1, 1], &cell_map, Some(&reserved), &mut rng, &sim);

        assert_eq!(found, None);
    }
}
