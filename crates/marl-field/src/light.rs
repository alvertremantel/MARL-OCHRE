use crate::field::Field;
use marl_config::*;
use std::collections::HashMap;

/// Light availability field. One scalar per voxel.
pub struct LightField {
    grid: GridDims,
    pub data: Vec<f32>,
}

impl LightField {
    pub fn new(grid: GridDims) -> Self {
        grid.validate().expect("invalid grid dimensions");
        Self {
            grid,
            data: vec![
                0.0;
                grid.voxel_count()
                    .expect("validated grid voxel count must fit usize")
            ],
        }
    }

    pub fn new_default() -> Self {
        Self::new(GridDims::default())
    }

    #[inline]
    pub fn grid(&self) -> GridDims {
        self.grid
    }

    #[inline]
    pub fn grid_x(&self) -> usize {
        self.grid.x
    }

    #[inline]
    pub fn grid_y(&self) -> usize {
        self.grid.y
    }

    #[inline]
    pub fn grid_z(&self) -> usize {
        self.grid.z
    }

    #[inline]
    pub fn voxel_count(&self) -> usize {
        self.grid
            .voxel_count()
            .expect("validated grid voxel count must fit usize")
    }

    #[inline]
    fn idx(&self, x: usize, y: usize, z: usize) -> usize {
        (z * self.grid.y + y) * self.grid.x + x
    }

    pub fn get(&self, x: usize, y: usize, z: usize) -> f32 {
        self.data[self.idx(x, y, z)]
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
        assert_eq!(
            field.grid(),
            self.grid,
            "light field grid dimensions must match chemical field grid dimensions"
        );
        for y in 0..self.grid.y {
            for x in 0..self.grid.x {
                let mut intensity = sim.surface_intensity.max(0.0);
                let cell_absorption = sim.cell_absorption.max(0.0);
                let chemical_absorption = sim.chemical_absorption.max(0.0);
                let light_floor = sim.light_floor.max(0.0);

                for z in 0..self.grid.z {
                    if intensity == 0.0 {
                        let idx = self.idx(x, y, z);
                        self.data[idx] = 0.0;
                        continue;
                    }

                    let idx = self.idx(x, y, z);
                    self.data[idx] = intensity;

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
        Self::new_default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn light_is_monotonic_down_column_and_floors_to_zero() {
        let grid = GridDims { x: 3, y: 2, z: 4 };
        let sim = SimulationConfig {
            surface_intensity: 1.0,
            cell_absorption: 2.0,
            chemical_absorption: 0.0,
            light_floor: 0.2,
            ..Default::default()
        };

        let field = Field::new(grid);
        let mut cells = HashMap::new();
        cells.insert([0, 0, 0], 0);
        cells.insert([0, 0, 1], 1);

        let mut light = LightField::new(grid);
        light.update(&field, &cells, &sim);

        assert_eq!(light.get(0, 0, 0), 1.0);
        assert!(light.get(0, 0, 1) < light.get(0, 0, 0));
        assert_eq!(light.get(0, 0, 2), 0.0);
        assert_eq!(light.get(0, 0, grid.z - 1), 0.0);
    }

    #[test]
    fn negative_absorber_does_not_amplify_light() {
        let grid = GridDims { x: 2, y: 2, z: 2 };
        let sim = SimulationConfig {
            surface_intensity: 1.0,
            cell_absorption: 0.0,
            chemical_absorption: 1.0,
            light_floor: 0.0,
            ..Default::default()
        };

        let mut field = Field::new(grid);
        field.set(0, 0, 0, 4, -10.0);

        let mut light = LightField::new(grid);
        light.update(&field, &HashMap::new(), &sim);

        assert_eq!(light.get(0, 0, 0), 1.0);
        assert_eq!(light.get(0, 0, 1), 1.0);
    }

    #[test]
    fn default_light_uses_default_grid_dimensions() {
        let light = LightField::default();

        assert_eq!(light.grid(), GridDims::default());
        assert_eq!(light.grid_x(), GRID_X);
        assert_eq!(light.grid_y(), GRID_Y);
        assert_eq!(light.grid_z(), GRID_Z);
        assert_eq!(light.voxel_count(), GRID_X * GRID_Y * GRID_Z);
        assert_eq!(light.data.len(), GRID_X * GRID_Y * GRID_Z);
    }
}
