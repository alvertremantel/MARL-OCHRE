// Naive f32 reaction-diffusion shader matching src/field.rs::diffusion_step_inner.
// Layout is AoS per voxel: [z][y][x][species].

const GRID_X: u32 = 128u;
const GRID_Y: u32 = 128u;
const GRID_Z: u32 = 64u;
const S_EXT: u32 = 12u;
const GRID_SIZE: u32 = GRID_X * GRID_Y * GRID_Z;
const STRUCTURAL_SPECIES: u32 = 7u;

struct DiffusionParams {
    dt_sub: f32,
    alpha_eps: f32,
    k_eps: f32,
    _pad0: f32,
    d_voxel: array<f32, 12>,
    lambda_decay: array<f32, 12>,
}

@group(0) @binding(0)
var<storage, read> field_in: array<f32>;

@group(0) @binding(1)
var<storage, read_write> field_out: array<f32>;

@group(0) @binding(2)
var<storage, read> occupancy: array<u32>;

@group(0) @binding(3)
var<storage, read> params: DiffusionParams;

fn field_idx(x: u32, y: u32, z: u32, s: u32) -> u32 {
    return ((z * GRID_Y + y) * GRID_X + x) * S_EXT + s;
}

fn voxel_idx(x: u32, y: u32, z: u32) -> u32 {
    return (z * GRID_Y + y) * GRID_X + x;
}

fn niche_factor(voxel: u32) -> f32 {
    let structural = max(field_in[voxel * S_EXT + STRUCTURAL_SPECIES], 0.0);
    let denom = max(params.k_eps + structural, 0.00000011920929);
    return clamp(1.0 - max(params.alpha_eps, 0.0) * structural / denom, 0.0, 1.0);
}

struct Neighbor {
    voxel: u32,
    niche: f32,
    present: u32,
}

fn neighbor_info(nx: i32, ny: i32, nz: i32) -> Neighbor {
    if (nx < 0 || nx >= i32(GRID_X) || ny < 0 || ny >= i32(GRID_Y) || nz < 0 || nz >= i32(GRID_Z)) {
        return Neighbor(0u, 0.0, 0u);
    }

    let x = u32(nx);
    let y = u32(ny);
    let z = u32(nz);
    let neighbor_voxel = voxel_idx(x, y, z);
    if (occupancy[neighbor_voxel] != 0u) {
        return Neighbor(0u, 0.0, 0u);
    }

    return Neighbor(neighbor_voxel, niche_factor(neighbor_voxel), 1u);
}

fn neighbor_flux(neighbor: Neighbor, s: u32, c: f32, base_d: f32, niche_here: f32) -> f32 {
    if (neighbor.present == 0u) {
        return 0.0;
    }

    let d_face = 0.5 * base_d * (niche_here + neighbor.niche);
    return d_face * (field_in[neighbor.voxel * S_EXT + s] - c);
}

@compute @workgroup_size(256)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let voxel = gid.x;
    if (voxel >= GRID_SIZE) {
        return;
    }

    let x = voxel % GRID_X;
    let y = (voxel / GRID_X) % GRID_Y;
    let z = voxel / (GRID_X * GRID_Y);
    let base = voxel * S_EXT;

    if (occupancy[voxel] != 0u) {
        for (var s = 0u; s < S_EXT; s = s + 1u) {
            field_out[base + s] = field_in[base + s];
        }
        return;
    }

    let niche_here = niche_factor(voxel);
    let xi = i32(x);
    let yi = i32(y);
    let zi = i32(z);
    let xm = neighbor_info(xi - 1, yi, zi);
    let xp = neighbor_info(xi + 1, yi, zi);
    let ym = neighbor_info(xi, yi - 1, zi);
    let yp = neighbor_info(xi, yi + 1, zi);
    let zm = neighbor_info(xi, yi, zi - 1);
    let zp = neighbor_info(xi, yi, zi + 1);

    for (var s = 0u; s < S_EXT; s = s + 1u) {
        let c = field_in[base + s];
        let base_d = max(params.d_voxel[s], 0.0);
        let diffusion =
            neighbor_flux(xm, s, c, base_d, niche_here) +
            neighbor_flux(xp, s, c, base_d, niche_here) +
            neighbor_flux(ym, s, c, base_d, niche_here) +
            neighbor_flux(yp, s, c, base_d, niche_here) +
            neighbor_flux(zm, s, c, base_d, niche_here) +
            neighbor_flux(zp, s, c, base_d, niche_here);
        let decay = params.lambda_decay[s] * c;
        let new_c = c + params.dt_sub * (diffusion - decay);
        field_out[base + s] = max(new_c, 0.0);
    }
}
