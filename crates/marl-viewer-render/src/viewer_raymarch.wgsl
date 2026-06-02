// MARL Isometric Viewer Raymarch Shader
//
// Orthographic ray/AABB traversal through a normalized 3D simulation volume.
// Binds:
//   @binding(0) field_tex: 3D R32Float — packed chemical field (tex_x = voxel_x * s_ext + species)
//   @binding(1) params:    uniform — ViewerParams struct
//   @binding(2) cell_tex:  3D Rgba8Uint — microbe occupancy/identity texture

struct Params {
    grid: vec4<u32>,        // grid_x, grid_y, grid_z, s_ext
    render: vec4<u32>,      // width, height, species, steps
    transfer: vec4<f32>,    // exposure, density_scale, cell_alpha, _unused
    axis_scale: vec4<f32>,  // grid / max_dim, 0
    cam_right: vec4<f32>,   // right.xyz, right.w = zoom
    cam_up: vec4<f32>,      // up.xyz, 0
    cam_dir: vec4<f32>,     // dir.xyz, 0
    options: vec4<u32>,     // options.x = cells_enabled
};

@group(0) @binding(0) var field_tex: texture_3d<f32>;
@group(0) @binding(1) var<uniform> params: Params;
@group(0) @binding(2) var cell_tex: texture_3d<u32>;

// ---------------------------------------------------------------------------
// Vertex shader — full-screen triangle from Phase 1
// ---------------------------------------------------------------------------

struct VertexOut {
    @builtin(position) position: vec4<f32>,
};

@vertex
fn vs_main(@builtin(vertex_index) vertex_index: u32) -> VertexOut {
    var positions = array<vec2<f32>, 3>(
        vec2<f32>(-1.0, -3.0),
        vec2<f32>(3.0, 1.0),
        vec2<f32>(-1.0, 1.0),
    );

    var out: VertexOut;
    out.position = vec4<f32>(positions[vertex_index], 0.0, 1.0);
    return out;
}

// ---------------------------------------------------------------------------
// Ray / AABB intersection (slab method)
// ---------------------------------------------------------------------------

fn intersect_box(
    origin: vec3<f32>,
    dir: vec3<f32>,
    box_min: vec3<f32>,
    box_max: vec3<f32>,
) -> bool {
    var t_near = -1e38;
    var t_far = 1e38;

    for (var i = 0u; i < 3u; i = i + 1u) {
        let d = dir[i];
        let o = origin[i];
        let mn = box_min[i];
        let mx = box_max[i];

        if (abs(d) < 1e-10) {
            // Ray is parallel to this slab; must be inside
            if (o < mn || o > mx) {
                return false;
            }
        } else {
            let inv_d = 1.0 / d;
            var t1 = (mn - o) * inv_d;
            var t2 = (mx - o) * inv_d;
            if (t1 > t2) {
                let tmp = t1;
                t1 = t2;
                t2 = tmp;
            }
            t_near = max(t_near, t1);
            t_far = min(t_far, t2);
        }
    }

    return t_near <= t_far && t_far > 0.0;
}

// ---------------------------------------------------------------------------
// Convert world position to voxel index
// ---------------------------------------------------------------------------

fn world_to_voxel(world: vec3<f32>) -> vec3<i32> {
    let ascl = params.axis_scale.xyz;
    let upper = 0.99999994;
    // world x in [-half.x, +half.x] → voxel x in [0, grid_x-1]
    let nx = clamp(world.x / ascl.x + 0.5, 0.0, upper);
    let vx = clamp(i32(floor(nx * f32(params.grid.x))), 0, i32(params.grid.x) - 1);
    // world y in [-half.y, +half.y] → voxel y in [0, grid_y-1]
    let ny = clamp(world.y / ascl.y + 0.5, 0.0, upper);
    let vy = clamp(i32(floor(ny * f32(params.grid.y))), 0, i32(params.grid.y) - 1);
    // world z in [+half.z, -half.z] → voxel z in [0, grid_z-1]
    //   top surface (tex_z=0) corresponds to world_z = +half.z
    let nz = clamp(0.5 - world.z / ascl.z, 0.0, upper);
    let vz = clamp(i32(floor(nz * f32(params.grid.z))), 0, i32(params.grid.z) - 1);
    return vec3<i32>(vx, vy, vz);
}

// The same voxel mapping, with x expanded for species-packed field sampling.
fn world_to_field_texel(world: vec3<f32>) -> vec3<i32> {
    let v = world_to_voxel(world);
    let species = i32(params.render.z);
    let s_ext = i32(params.grid.w);
    let tex_x = v.x * s_ext + min(species, s_ext - 1);
    return vec3<i32>(tex_x, v.y, v.z);
}

// ---------------------------------------------------------------------------
// Colour palette for chemical density (similar to Phase 1)
// ---------------------------------------------------------------------------

fn palette(value: f32, z_frac: f32) -> vec3<f32> {
    let cold = vec3<f32>(0.03, 0.18, 0.45);
    let warm = vec3<f32>(0.95, 0.62, 0.18);
    let hot = vec3<f32>(0.85, 0.95, 0.92);
    let mid = mix(cold, warm, clamp(value * 1.8, 0.0, 1.0));
    return mix(mid, hot, clamp(value - 0.55, 0.0, 1.0)) * (0.65 + 0.35 * z_frac);
}

// ---------------------------------------------------------------------------
// Fragment shader — isometric ray/AABB traversal + cell compositing
// ---------------------------------------------------------------------------

@fragment
fn fs_main(@builtin(position) frag_coord: vec4<f32>) -> @location(0) vec4<f32> {
    let viewport = vec2<f32>(f32(params.render.x), f32(params.render.y));
    // Normalised screen coordinates in [-1, 1] centred, aspect-correct
    // wgpu: y=0 is top, y=viewport.y is bottom, so we flip y
    var screen = (frag_coord.xy / max(viewport, vec2<f32>(1.0, 1.0))) * 2.0 - 1.0;
    // Correct aspect: scale x by (width/height)
    let aspect = viewport.x / max(viewport.y, 1.0);
    screen.x *= aspect;
    // Flip y so +y = up
    screen.y = -screen.y;

    // Unpack camera
    let zoom = params.cam_right.w;
    let right = params.cam_right.xyz;
    let up = params.cam_up.xyz;
    let dir = params.cam_dir.xyz;

    // Orthographic ray
    let origin = right * screen.x * zoom + up * screen.y * zoom - dir * 2.0;
    let ray_dir = dir;

    // Normalised simulation box centred at origin
    let half_box = 0.5 * params.axis_scale.xyz;
    let box_min = -half_box;
    let box_max = half_box;

    // Early out: pixels outside the volume projection
    if (!intersect_box(origin, ray_dir, box_min, box_max)) {
        let background = vec3<f32>(0.008, 0.01, 0.014);
        return vec4<f32>(background, 1.0);
    }

    // Effective step count: at least enough to cover the longest box diagonal
    // without skipping one-voxel cell markers.
    let max_grid = max(max(params.grid.x, params.grid.y), params.grid.z);
    let max_dim_f = max(f32(max_grid), params.axis_scale.x * 2.0);
    let max_dim = u32(max_dim_f);
    let effective_steps = max(params.render.w, 2u * max_dim);

    let species_count = params.grid.w;
    let density_scale = params.transfer.y;
    let exposure = params.transfer.x;
    let cell_alpha = params.transfer.z;
    let cells_enabled = params.options.x != 0u;

    // March distance: from box entry to exit
    // We'll use a fixed step size that covers the full box traversal.
    // Box extent along ray direction:
    // For orthographic dir, the total travel through the box is the distance
    // along the ray from the near intersection to far intersection.
    // We just march from a point behind the box to a point beyond it.
    let step_count = max(effective_steps, 1u);
    // Compute approximate near/far along ray
    var near = -1e38;
    var far = 1e38;
    for (var i = 0u; i < 3u; i = i + 1u) {
        let d = ray_dir[i];
        let o = origin[i];
        if (abs(d) < 1e-10) { continue; }
        let inv_d = 1.0 / d;
        var t1 = (box_min[i] - o) * inv_d;
        var t2 = (box_max[i] - o) * inv_d;
        if (t1 > t2) {
            let tmp = t1; t1 = t2; t2 = tmp;
        }
        near = max(near, t1);
        far = min(far, t2);
    }
    near = max(near, 0.0);
    if (near >= far || far <= 0.0) {
        let background = vec3<f32>(0.008, 0.01, 0.014);
        return vec4<f32>(background, 1.0);
    }

    let total_dist = far - near;
    let step_dist = total_dist / f32(step_count);

    var rgb = vec3<f32>(0.0, 0.0, 0.0);
    var alpha = 0.0;

    // Track previous voxel for cell deduplication
    var prev_voxel = vec3<i32>(-999, -999, -999);

    for (var step = 0u; step < step_count; step = step + 1u) {
        let t = near + (f32(step) + 0.5) * step_dist;
        let world_p = origin + ray_dir * t;

        // Sample field density
        let field_texel = world_to_field_texel(world_p);
        let raw = textureLoad(field_tex, field_texel, 0).r;
        let density = clamp(raw * density_scale, 0.0, 1.0);

        // Composite field density
        let sample_alpha = 1.0 - exp(-density * exposure / f32(step_count));
        let z_frac = clamp((world_p.z - box_min.z) / (box_max.z - box_min.z + 1e-10), 0.0, 1.0);
        let field_rgb = palette(density, 1.0 - z_frac);

        rgb = rgb + (1.0 - alpha) * sample_alpha * field_rgb;
        alpha = alpha + (1.0 - alpha) * sample_alpha;

        // Composite cell voxel (only when voxel changes)
        if (cells_enabled) {
            let voxel = world_to_voxel(world_p);
            if (voxel.x != prev_voxel.x || voxel.y != prev_voxel.y || voxel.z != prev_voxel.z) {
                prev_voxel = voxel;
                let gx = i32(params.grid.x);
                let gy = i32(params.grid.y);
                let gz = i32(params.grid.z);
                if (voxel.x >= 0 && voxel.x < gx &&
                    voxel.y >= 0 && voxel.y < gy &&
                    voxel.z >= 0 && voxel.z < gz)
                {
                    let cell_raw = textureLoad(cell_tex, voxel, 0);
                    let cell_a = f32(cell_raw.a); // 0-255

                    if (cell_a > 0.0) {
                        let cr = f32(cell_raw.r) / 255.0;
                        let cg = f32(cell_raw.g) / 255.0;
                        let cb = f32(cell_raw.b) / 255.0;
                        // User-controlled alpha (0-1) multiplied by texture alpha (0-255/255)
                        let composite_a = cell_a / 255.0 * cell_alpha;

                        let cell_rgb = vec3<f32>(cr, cg, cb);
                        rgb = rgb + (1.0 - alpha) * composite_a * cell_rgb;
                        alpha = alpha + (1.0 - alpha) * composite_a;
                    }
                }
            }
        }

        if (alpha > 0.985) {
            break;
        }
    }

    let background = vec3<f32>(0.008, 0.01, 0.014);
    return vec4<f32>(rgb + (1.0 - alpha) * background, 1.0);
}
