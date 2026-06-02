use crate::context::{GpuContext, GpuError};
use marl_config::{GRID_X, GRID_Y, GRID_Z, S_EXT, SimulationConfig};
use marl_field::field::{Field, stable_diffusion_substeps, validate_diffusion_config};

const GRID_SIZE: usize = GRID_X * GRID_Y * GRID_Z;
const FIELD_FLOATS: usize = GRID_SIZE * S_EXT;
const WORKGROUP_SIZE: u32 = 256;

#[repr(C)]
#[derive(Clone, Copy, Debug, bytemuck::Pod, bytemuck::Zeroable)]
struct DiffusionParams {
    dt_sub: f32,
    alpha_eps: f32,
    k_eps: f32,
    _pad0: f32,
    d_voxel: [f32; S_EXT],
    lambda_decay: [f32; S_EXT],
}

pub struct GpuFieldDiffuser {
    context: GpuContext,
    pipeline: wgpu::ComputePipeline,
    bind_group_a_to_b: wgpu::BindGroup,
    bind_group_b_to_a: wgpu::BindGroup,
    field_buffer_a: wgpu::Buffer,
    field_buffer_b: wgpu::Buffer,
    occupancy_buffer: wgpu::Buffer,
    params_buffer: wgpu::Buffer,
    staging_buffer: wgpu::Buffer,
    occupancy_upload: Vec<u32>,
}

impl GpuFieldDiffuser {
    pub fn new() -> Result<Self, GpuError> {
        if GRID_X != 128 || GRID_Y != 128 || GRID_Z != 64 || S_EXT != 12 {
            return Err(GpuError::InvalidInput(format!(
                "GPU shader constants require 128x128x64 with 12 species, got {GRID_X}x{GRID_Y}x{GRID_Z} with {S_EXT} species"
            )));
        }

        let context = GpuContext::new()?;
        let device = &context.device;

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("MARL Field Diffusion Shader"),
            source: wgpu::ShaderSource::Wgsl(include_str!("shaders/field_diffuse.wgsl").into()),
        });

        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("MARL Field Diffusion Bind Group Layout"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: false },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 3,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
            ],
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("MARL Field Diffusion Pipeline Layout"),
            bind_group_layouts: &[Some(&bind_group_layout)],
            immediate_size: 0,
        });

        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("MARL Field Diffusion Pipeline"),
            layout: Some(&pipeline_layout),
            module: &shader,
            entry_point: Some("main"),
            compilation_options: wgpu::PipelineCompilationOptions::default(),
            cache: None,
        });

        let field_bytes = field_buffer_bytes();
        let occupancy_bytes = (GRID_SIZE * std::mem::size_of::<u32>()) as u64;

        let field_usage = wgpu::BufferUsages::STORAGE
            | wgpu::BufferUsages::COPY_DST
            | wgpu::BufferUsages::COPY_SRC;

        let field_buffer_a = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("MARL Field Buffer A"),
            size: field_bytes,
            usage: field_usage,
            mapped_at_creation: false,
        });
        let field_buffer_b = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("MARL Field Buffer B"),
            size: field_bytes,
            usage: field_usage,
            mapped_at_creation: false,
        });
        let occupancy_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("MARL Occupancy Buffer"),
            size: occupancy_bytes,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let params_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("MARL Diffusion Params Buffer"),
            size: std::mem::size_of::<DiffusionParams>() as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let staging_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("MARL Field Readback Buffer"),
            size: field_bytes,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });
        let bind_group_a_to_b = create_diffusion_bind_group(
            device,
            &bind_group_layout,
            &field_buffer_a,
            &field_buffer_b,
            &occupancy_buffer,
            &params_buffer,
            "MARL Field Diffusion Bind Group A to B",
        );
        let bind_group_b_to_a = create_diffusion_bind_group(
            device,
            &bind_group_layout,
            &field_buffer_b,
            &field_buffer_a,
            &occupancy_buffer,
            &params_buffer,
            "MARL Field Diffusion Bind Group B to A",
        );

        Ok(Self {
            context,
            pipeline,
            bind_group_a_to_b,
            bind_group_b_to_a,
            field_buffer_a,
            field_buffer_b,
            occupancy_buffer,
            params_buffer,
            staging_buffer,
            occupancy_upload: vec![0; GRID_SIZE],
        })
    }

    pub fn diffuse_tick_with_cells(
        &mut self,
        field: &mut Field,
        occupancy: &[bool],
        sim: &SimulationConfig,
    ) -> Result<(), GpuError> {
        if field.data.len() != FIELD_FLOATS {
            return Err(GpuError::InvalidInput(format!(
                "field has {} floats, expected {FIELD_FLOATS}",
                field.data.len()
            )));
        }
        if occupancy.len() != GRID_SIZE {
            return Err(GpuError::InvalidInput(format!(
                "occupancy has {} voxels, expected {GRID_SIZE}",
                occupancy.len()
            )));
        }
        validate_diffusion_config(sim).map_err(GpuError::InvalidInput)?;
        let substeps = stable_diffusion_substeps(sim);
        if substeps == 0 {
            return Ok(());
        }

        let params = DiffusionParams {
            dt_sub: sim.dt / substeps as f32,
            alpha_eps: sim.alpha_eps,
            k_eps: sim.k_eps,
            _pad0: 0.0,
            d_voxel: sim.d_voxel,
            lambda_decay: sim.lambda_decay,
        };
        for (dst, &occupied) in self.occupancy_upload.iter_mut().zip(occupancy) {
            *dst = u32::from(occupied);
        }

        self.context
            .queue
            .write_buffer(&self.field_buffer_a, 0, bytemuck::cast_slice(&field.data));
        self.context.queue.write_buffer(
            &self.occupancy_buffer,
            0,
            bytemuck::cast_slice(&self.occupancy_upload),
        );
        self.context
            .queue
            .write_buffer(&self.params_buffer, 0, bytemuck::bytes_of(&params));

        let mut encoder =
            self.context
                .device
                .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                    label: Some("MARL Field Diffusion Encoder"),
                });

        for substep in 0..substeps {
            let bind_group = if substep % 2 == 0 {
                &self.bind_group_a_to_b
            } else {
                &self.bind_group_b_to_a
            };

            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("MARL Field Diffusion Pass"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&self.pipeline);
            pass.set_bind_group(0, bind_group, &[]);
            pass.dispatch_workgroups((GRID_SIZE as u32).div_ceil(WORKGROUP_SIZE), 1, 1);
        }

        let final_buffer = if substeps.is_multiple_of(2) {
            &self.field_buffer_a
        } else {
            &self.field_buffer_b
        };
        encoder.copy_buffer_to_buffer(
            final_buffer,
            0,
            &self.staging_buffer,
            0,
            field_buffer_bytes(),
        );
        self.context.queue.submit(Some(encoder.finish()));

        self.readback_field(field)
    }

    fn readback_field(&self, field: &mut Field) -> Result<(), GpuError> {
        let slice = self.staging_buffer.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |result| {
            let _ = tx.send(result.map_err(|e| e.to_string()));
        });

        self.context
            .device
            .poll(wgpu::PollType::wait_indefinitely())
            .map_err(|e| GpuError::BufferMap(e.to_string()))?;
        rx.recv()
            .map_err(|e| GpuError::BufferMap(e.to_string()))?
            .map_err(GpuError::BufferMap)?;

        {
            let view = slice.get_mapped_range();
            let data: &[f32] = bytemuck::cast_slice(&view);
            field.data.copy_from_slice(data);
        }
        self.staging_buffer.unmap();
        Ok(())
    }
}

fn create_diffusion_bind_group(
    device: &wgpu::Device,
    layout: &wgpu::BindGroupLayout,
    input: &wgpu::Buffer,
    output: &wgpu::Buffer,
    occupancy: &wgpu::Buffer,
    params: &wgpu::Buffer,
    label: &'static str,
) -> wgpu::BindGroup {
    device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some(label),
        layout,
        entries: &[
            wgpu::BindGroupEntry {
                binding: 0,
                resource: input.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 1,
                resource: output.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 2,
                resource: occupancy.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 3,
                resource: params.as_entire_binding(),
            },
        ],
    })
}

fn field_buffer_bytes() -> u64 {
    (FIELD_FLOATS * std::mem::size_of::<f32>()) as u64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diffusion_params_layout_is_stable() {
        assert_eq!(std::mem::size_of::<DiffusionParams>(), 112);
        assert_eq!(std::mem::offset_of!(DiffusionParams, dt_sub), 0);
        assert_eq!(std::mem::offset_of!(DiffusionParams, alpha_eps), 4);
        assert_eq!(std::mem::offset_of!(DiffusionParams, k_eps), 8);
        assert_eq!(std::mem::offset_of!(DiffusionParams, _pad0), 12);
        assert_eq!(std::mem::offset_of!(DiffusionParams, d_voxel), 16);
        assert_eq!(std::mem::offset_of!(DiffusionParams, lambda_decay), 64);
    }
}
