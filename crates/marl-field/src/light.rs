use crate::field::Field;
use marl_config::*;
use std::collections::HashMap;

/// Light availability field. One scalar per voxel.
pub struct LightField {
    pub data: Vec<f32>,
}

impl LightField {
    pub fn new() -> Self {
        Self {
            data: vec![0.0; GRID_X * GRID_Y * GRID_Z],
        }
    }

    #[inline]
    fn idx(x: usize, y: usize, z: usize) -> usize {
        (z * GRID_Y + y) * GRID_X + x
    }

    pub fn get(&self, x: usize, y: usize, z: usize) -> f32 {
        self.data[Self::idx(x, y, z)]
    }

    /// Beer-Lambert top-down sweep.
    /// Light enters at z=0 with intensity 1.0, attenuates with depth
    /// based on cell density and absorber concentrations.
    pub fn update(
        &mut self,
        field: &Field,
        cells: &HashMap<[u16; 3], usize>, // pos -> cell index (for density)
        sim: &SimulationConfig,
    ) {
        for y in 0..GRID_Y {
            for x in 0..GRID_X {
                let mut intensity = sim.surface_intensity.max(0.0);
                let cell_absorption = sim.cell_absorption.max(0.0);
                let chemical_absorption = sim.chemical_absorption.max(0.0);
                let light_floor = sim.light_floor.max(0.0);

                for z in 0..GRID_Z {
                    if intensity == 0.0 {
                        self.data[Self::idx(x, y, z)] = 0.0;
                        continue;
                    }

                    self.data[Self::idx(x, y, z)] = intensity;

                    // Attenuation from cells
                    let pos = [x as u16, y as u16, z as u16];
                    if cells.contains_key(&pos) {
                        intensity *= (-cell_absorption).exp();
                    }

                    // Attenuation from chemical absorbers (e.g., organic waste = species 4)
                    let absorber = field.get(x, y, z, 4).max(0.0);
                    intensity *= (-chemical_absorption * absorber).exp();

                    // Floor to prevent denormals
                    if intensity < light_floor {
                        intensity = 0.0;
                    }
                }
            }
        }
    }
}

impl Default for LightField {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn light_is_monotonic_down_column_and_floors_to_zero() {
        let sim = SimulationConfig {
            surface_intensity: 1.0,
            cell_absorption: 2.0,
            chemical_absorption: 0.0,
            light_floor: 0.2,
            ..Default::default()
        };

        let field = Field::new();
        let mut cells = HashMap::new();
        cells.insert([0, 0, 0], 0);
        cells.insert([0, 0, 1], 1);

        let mut light = LightField::new();
        light.update(&field, &cells, &sim);

        assert_eq!(light.get(0, 0, 0), 1.0);
        assert!(light.get(0, 0, 1) < light.get(0, 0, 0));
        assert_eq!(light.get(0, 0, 2), 0.0);
        assert_eq!(light.get(0, 0, GRID_Z - 1), 0.0);
    }

    #[test]
    fn negative_absorber_does_not_amplify_light() {
        let sim = SimulationConfig {
            surface_intensity: 1.0,
            cell_absorption: 0.0,
            chemical_absorption: 1.0,
            light_floor: 0.0,
            ..Default::default()
        };

        let mut field = Field::new();
        field.set(0, 0, 0, 4, -10.0);

        let mut light = LightField::new();
        light.update(&field, &HashMap::new(), &sim);

        assert_eq!(light.get(0, 0, 0), 1.0);
        assert_eq!(light.get(0, 0, 1), 1.0);
    }
}
