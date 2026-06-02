use std::borrow::Cow;
use std::error::Error;
use std::path::PathBuf;
use std::sync::Arc;

use wgpu::util::DeviceExt;
use winit::dpi::PhysicalSize;
use winit::event::WindowEvent;
use winit::window::Window;

#[cfg(test)]
use marl_viewer_core::args::ViewMode;
use marl_viewer_core::args::{CellMode, ViewerArgs};
use marl_viewer_core::camera::{CameraBasis, camera_basis};
use marl_viewer_core::io::{
    LoadedCell, SnapshotPayload, discover_field_ticks, load_run_meta, load_snapshot,
};
use marl_viewer_core::types::{GuiAction, SnapshotInfo, choose_initial_tick, neighbor_tick};

use crate::gui::GuiState;

// ---------------------------------------------------------------------------
// ViewerParams — must be #[repr(C)], Pod, and mirror WGSL Params exactly
// ---------------------------------------------------------------------------

#[repr(C)]
#[derive(Clone, Copy, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub(crate) struct ViewerParams {
    pub(crate) grid: [u32; 4],       // grid_x, grid_y, grid_z, s_ext
    pub(crate) render: [u32; 4],     // width, height, species, steps
    pub(crate) transfer: [f32; 4],   // exposure, density_scale, cell_alpha, _unused
    pub(crate) axis_scale: [f32; 4], // grid / max_dim, 0
    pub(crate) cam_right: [f32; 4],  // right.xyz, right.w = zoom
    pub(crate) cam_up: [f32; 4],     // up.xyz, 0
    pub(crate) cam_dir: [f32; 4],    // dir.xyz, 0
    pub(crate) options: [u32; 4],    // options.x = cells_enabled
}

// ---------------------------------------------------------------------------
// SnapshotGpuResources — reloadable GPU state for the active snapshot
// ---------------------------------------------------------------------------

struct SnapshotGpuResources {
    bind_group: wgpu::BindGroup,
    params_buffer: wgpu::Buffer,
    params: ViewerParams,
    _field_texture: wgpu::Texture,
    _field_view: wgpu::TextureView,
    _cell_texture: wgpu::Texture,
    _cell_view: wgpu::TextureView,
}

// ---------------------------------------------------------------------------
// Renderer
// ---------------------------------------------------------------------------

pub struct Renderer {
    pub window: Arc<Window>,
    surface: wgpu::Surface<'static>,
    device: wgpu::Device,
    queue: wgpu::Queue,
    config: wgpu::SurfaceConfiguration,
    pipeline: wgpu::RenderPipeline,
    bind_group_layout: wgpu::BindGroupLayout,
    snapshot: SnapshotGpuResources,
    pub(crate) loaded_info: Option<SnapshotInfo>,
    pub(crate) args: ViewerArgs,
    // egui integration
    egui_ctx: egui::Context,
    egui_state: egui_winit::State,
    egui_renderer: egui_wgpu::Renderer,
    gui: GuiState,
}

impl Renderer {
    pub async fn new(window: Arc<Window>, args: ViewerArgs) -> Result<Self, Box<dyn Error>> {
        let size = window.inner_size();
        let instance = wgpu::Instance::default();
        let surface = instance.create_surface(window.clone())?;
        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                compatible_surface: Some(&surface),
                force_fallback_adapter: false,
            })
            .await
            .map_err(|_| "no compatible GPU adapter found")?;

        let info = adapter.get_info();
        eprintln!("[viewer] using adapter: {} ({:?})", info.name, info.backend);

        let (device, queue) = adapter
            .request_device(&wgpu::DeviceDescriptor {
                label: Some("MARL Viewer Device"),
                required_features: wgpu::Features::empty(),
                required_limits: wgpu::Limits::default(),
                memory_hints: wgpu::MemoryHints::Performance,
                trace: wgpu::Trace::Off,
                ..Default::default()
            })
            .await?;

        let surface_caps = surface.get_capabilities(&adapter);
        let surface_format = surface_caps
            .formats
            .iter()
            .copied()
            .find(|format| format.is_srgb())
            .unwrap_or(surface_caps.formats[0]);
        let present_mode = surface_caps
            .present_modes
            .iter()
            .copied()
            .find(|mode| *mode == wgpu::PresentMode::Fifo)
            .unwrap_or(surface_caps.present_modes[0]);
        let alpha_mode = surface_caps.alpha_modes[0];

        let config = wgpu::SurfaceConfiguration {
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            format: surface_format,
            width: size.width.max(1),
            height: size.height.max(1),
            present_mode,
            alpha_mode,
            view_formats: vec![],
            desired_maximum_frame_latency: 2,
        };
        surface.configure(&device, &config);

        // Bind group layout (stored for reloading)
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("MARL Viewer Bind Group Layout"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Texture {
                        sample_type: wgpu::TextureSampleType::Float { filterable: false },
                        view_dimension: wgpu::TextureViewDimension::D3,
                        multisampled: false,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Texture {
                        sample_type: wgpu::TextureSampleType::Uint,
                        view_dimension: wgpu::TextureViewDimension::D3,
                        multisampled: false,
                    },
                    count: None,
                },
            ],
        });

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("MARL Viewer Raymarch Shader"),
            source: wgpu::ShaderSource::Wgsl(Cow::Borrowed(include_str!("viewer_raymarch.wgsl"))),
        });
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("MARL Viewer Pipeline Layout"),
            bind_group_layouts: &[Some(&bind_group_layout)],
            immediate_size: 0,
        });
        let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("MARL Viewer Raymarch Pipeline"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: Some("vs_main"),
                buffers: &[],
                compilation_options: wgpu::PipelineCompilationOptions::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader,
                entry_point: Some("fs_main"),
                targets: &[Some(wgpu::ColorTargetState {
                    format: surface_format,
                    blend: Some(wgpu::BlendState::REPLACE),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
                compilation_options: wgpu::PipelineCompilationOptions::default(),
            }),
            primitive: wgpu::PrimitiveState::default(),
            depth_stencil: None,
            multisample: wgpu::MultisampleState::default(),
            multiview_mask: None,
            cache: None,
        });

        // -------------------------------------------------------------------
        // Initialize egui
        // -------------------------------------------------------------------
        let egui_ctx = egui::Context::default();
        let egui_state = egui_winit::State::new(
            egui_ctx.clone(),
            egui::ViewportId::ROOT,
            window.as_ref(),
            Some(window.scale_factor() as f32),
            window.theme(),
            Some(device.limits().max_texture_dimension_2d as usize),
        );
        let egui_renderer = egui_wgpu::Renderer::new(
            &device,
            config.format,
            egui_wgpu::RendererOptions::default(),
        );
        let gui = GuiState::new(&args);

        // Try to load the initial snapshot; fall back to placeholder on failure.
        let (snapshot, loaded_info) = match try_load_snapshot_resources(
            &device,
            &queue,
            &bind_group_layout,
            &args,
            config.width,
            config.height,
        ) {
            Ok((res, info)) => {
                eprintln!(
                    "[viewer] loaded tick {} species {} view {:?} cells {:?} ({} field bytes, {} cells)",
                    info.tick,
                    info.species,
                    info.view_mode,
                    info.cell_mode,
                    info.field_bytes,
                    info.cell_count
                );
                (res, Some(info))
            }
            Err(e) => {
                eprintln!("[viewer] initial snapshot load failed: {e}");
                let snap = create_placeholder_resources(
                    &device,
                    &queue,
                    &bind_group_layout,
                    &args,
                    config.width,
                    config.height,
                )?;
                (snap, None)
            }
        };

        // Set window title
        if let Some(ref info) = loaded_info {
            window.set_title(&format!(
                "MARL Viewer - tick {} species {} view {:?} cells {:?} ({} cells, {}x{}x{})",
                info.tick,
                info.species,
                info.view_mode,
                info.cell_mode,
                info.cell_count,
                info.grid[0],
                info.grid[1],
                info.grid[2]
            ));
        } else {
            window.set_title("MARL Viewer - no snapshot loaded");
        }

        Ok(Self {
            window,
            surface,
            device,
            queue,
            config,
            pipeline,
            bind_group_layout,
            snapshot,
            loaded_info,
            args,
            egui_ctx,
            egui_state,
            egui_renderer,
            gui,
        })
    }

    /// Forward a window event to egui. Returns `true` if egui consumed the event
    /// and a repaint is needed.
    pub fn handle_window_event(&mut self, event: &WindowEvent) -> bool {
        let response = self.egui_state.on_window_event(self.window.as_ref(), event);
        response.repaint
    }

    pub fn resize(&mut self, size: PhysicalSize<u32>) {
        if size.width == 0 || size.height == 0 {
            return;
        }

        self.config.width = size.width;
        self.config.height = size.height;
        self.surface.configure(&self.device, &self.config);
        self.snapshot.params.render[0] = size.width;
        self.snapshot.params.render[1] = size.height;
        self.queue.write_buffer(
            &self.snapshot.params_buffer,
            0,
            bytemuck::bytes_of(&self.snapshot.params),
        );
    }

    pub fn render(&mut self) -> RenderResult {
        // -------------------------------------------------------------------
        // 1. Collect egui input and build GUI
        // -------------------------------------------------------------------
        let raw_input = self.egui_state.take_egui_input(self.window.as_ref());
        let egui_ctx = self.egui_ctx.clone();
        egui_ctx.begin_pass(raw_input);

        let loaded_ref = self.loaded_info.as_ref();
        let actions = self.gui.show(&egui_ctx, loaded_ref, &self.args);

        let full_output = egui_ctx.end_pass();

        // Handle platform output (cursor changes, etc.)
        self.egui_state
            .handle_platform_output(self.window.as_ref(), full_output.platform_output.clone());

        // Clipped primitives for rendering
        let paint_jobs = egui_ctx.tessellate(full_output.shapes, full_output.pixels_per_point);
        let screen_descriptor = egui_wgpu::ScreenDescriptor {
            size_in_pixels: [self.config.width, self.config.height],
            pixels_per_point: full_output.pixels_per_point,
        };

        // Update textures
        for (id, image_delta) in &full_output.textures_delta.set {
            self.egui_renderer
                .update_texture(&self.device, &self.queue, *id, image_delta);
        }

        // -------------------------------------------------------------------
        // 2. Process GUI actions
        // -------------------------------------------------------------------
        for action in actions {
            match action {
                GuiAction::OpenDirectoryDialog => {
                    let current_dir = if self.args.output_dir.is_dir() {
                        Some(self.args.output_dir.clone())
                    } else {
                        None
                    };
                    let picked = rfd::FileDialog::new()
                        .set_directory(current_dir.unwrap_or_else(|| PathBuf::from(".")))
                        .pick_folder();
                    if let Some(dir) = picked {
                        self.gui.set_info(format!("loading {}…", dir.display()));
                        match self.load_directory_from_gui(dir) {
                            Ok(()) => {}
                            Err(e) => self.gui.set_error(format!("load failed: {e}")),
                        }
                    }
                }
                GuiAction::LoadDirectory(dir) => {
                    self.gui.set_info(format!("loading {}…", dir.display()));
                    match self.load_directory_from_gui(dir) {
                        Ok(()) => {}
                        Err(e) => self.gui.set_error(format!("load failed: {e}")),
                    }
                }
                GuiAction::LoadTick(tick) => {
                    let new_args = ViewerArgs {
                        tick,
                        ..self.args.clone()
                    };
                    match self.apply_args(new_args) {
                        Ok(info) => {
                            let ticks =
                                discover_field_ticks(&self.args.output_dir).unwrap_or_default();
                            self.gui.sync_loaded(&info, &self.args, ticks);
                        }
                        Err(e) => {
                            self.gui
                                .set_error(format!("failed to load tick {tick}: {e}"));
                        }
                    }
                }
                GuiAction::ReloadCurrent => {
                    self.gui.set_info("reloading…".to_string());
                    match self.apply_args(self.args.clone()) {
                        Ok(info) => {
                            let ticks =
                                discover_field_ticks(&self.args.output_dir).unwrap_or_default();
                            self.gui.sync_loaded(&info, &self.args, ticks);
                        }
                        Err(e) => {
                            self.gui.set_error(format!("reload failed: {e}"));
                        }
                    }
                }
                GuiAction::FirstTick => {
                    if let Some(&tick) = self.gui.available_ticks.first() {
                        self.gui.tick_text = tick.to_string();
                        let new_args = ViewerArgs {
                            tick,
                            ..self.args.clone()
                        };
                        match self.apply_args(new_args) {
                            Ok(info) => {
                                let ticks =
                                    discover_field_ticks(&self.args.output_dir).unwrap_or_default();
                                self.gui.sync_loaded(&info, &self.args, ticks);
                            }
                            Err(e) => {
                                self.gui
                                    .set_error(format!("failed to load tick {tick}: {e}"));
                            }
                        }
                    }
                }
                GuiAction::LastTick => {
                    if let Some(&tick) = self.gui.available_ticks.last() {
                        self.gui.tick_text = tick.to_string();
                        let new_args = ViewerArgs {
                            tick,
                            ..self.args.clone()
                        };
                        match self.apply_args(new_args) {
                            Ok(info) => {
                                let ticks =
                                    discover_field_ticks(&self.args.output_dir).unwrap_or_default();
                                self.gui.sync_loaded(&info, &self.args, ticks);
                            }
                            Err(e) => {
                                self.gui
                                    .set_error(format!("failed to load tick {tick}: {e}"));
                            }
                        }
                    }
                }
                GuiAction::PrevTick => {
                    let current = self.args.tick;
                    if let Some(tick) = neighbor_tick(current, &self.gui.available_ticks, -1) {
                        self.gui.tick_text = tick.to_string();
                        let new_args = ViewerArgs {
                            tick,
                            ..self.args.clone()
                        };
                        match self.apply_args(new_args) {
                            Ok(info) => {
                                let ticks =
                                    discover_field_ticks(&self.args.output_dir).unwrap_or_default();
                                self.gui.sync_loaded(&info, &self.args, ticks);
                            }
                            Err(e) => {
                                self.gui
                                    .set_error(format!("failed to load tick {tick}: {e}"));
                            }
                        }
                    }
                }
                GuiAction::NextTick => {
                    let current = self.args.tick;
                    if let Some(tick) = neighbor_tick(current, &self.gui.available_ticks, 1) {
                        self.gui.tick_text = tick.to_string();
                        let new_args = ViewerArgs {
                            tick,
                            ..self.args.clone()
                        };
                        match self.apply_args(new_args) {
                            Ok(info) => {
                                let ticks =
                                    discover_field_ticks(&self.args.output_dir).unwrap_or_default();
                                self.gui.sync_loaded(&info, &self.args, ticks);
                            }
                            Err(e) => {
                                self.gui
                                    .set_error(format!("failed to load tick {tick}: {e}"));
                            }
                        }
                    }
                }
                GuiAction::ApplyViewSettings => {
                    let s_ext = self.loaded_info.as_ref().map(|i| i.s_ext);
                    match self.gui.build_view_args_from_drafts(&self.args, s_ext) {
                        Ok(new_args) => {
                            self.gui.set_info("applying view settings…".to_string());
                            match self.apply_args(new_args) {
                                Ok(info) => {
                                    let ticks = discover_field_ticks(&self.args.output_dir)
                                        .unwrap_or_default();
                                    self.gui.sync_loaded(&info, &self.args, ticks);
                                }
                                Err(e) => {
                                    self.gui.set_error(format!("apply failed: {e}"));
                                }
                            }
                        }
                        Err(msg) => {
                            self.gui.set_error(msg);
                        }
                    }
                }
                GuiAction::ResetDraftFromLoaded => {
                    self.gui.reset_drafts_from_args(&self.args);
                }
            }
        }

        // -------------------------------------------------------------------
        // 3. Acquire surface frame
        // -------------------------------------------------------------------
        let (frame, should_reconfigure) = match self.surface.get_current_texture() {
            wgpu::CurrentSurfaceTexture::Success(frame) => (frame, false),
            wgpu::CurrentSurfaceTexture::Suboptimal(frame) => (frame, true),
            wgpu::CurrentSurfaceTexture::Timeout | wgpu::CurrentSurfaceTexture::Occluded => {
                return RenderResult::Skip;
            }
            wgpu::CurrentSurfaceTexture::Outdated
            | wgpu::CurrentSurfaceTexture::Lost
            | wgpu::CurrentSurfaceTexture::Validation => return RenderResult::Reconfigure,
        };
        let view = frame
            .texture
            .create_view(&wgpu::TextureViewDescriptor::default());
        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("MARL Viewer Render Encoder"),
            });

        // Upload egui vertex/index data before render passes
        self.egui_renderer.update_buffers(
            &self.device,
            &self.queue,
            &mut encoder,
            &paint_jobs,
            &screen_descriptor,
        );

        // -------------------------------------------------------------------
        // 4. Raymarch pass (first, with Clear)
        // -------------------------------------------------------------------
        {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("MARL Viewer Raymarch Pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view,
                    depth_slice: None,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color {
                            r: 0.008,
                            g: 0.01,
                            b: 0.014,
                            a: 1.0,
                        }),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: None,
                timestamp_writes: None,
                occlusion_query_set: None,
                multiview_mask: None,
            });
            pass.set_pipeline(&self.pipeline);
            pass.set_bind_group(0, &self.snapshot.bind_group, &[]);
            pass.draw(0..3, 0..1);
        }

        // -------------------------------------------------------------------
        // 5. Egui pass (second, with Load)
        // -------------------------------------------------------------------
        {
            let pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("MARL Viewer Egui Pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view,
                    depth_slice: None,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Load,
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: None,
                timestamp_writes: None,
                occlusion_query_set: None,
                multiview_mask: None,
            });
            // egui_wgpu::Renderer::render expects RenderPass<'static>
            let mut static_pass = pass.forget_lifetime();
            self.egui_renderer
                .render(&mut static_pass, &paint_jobs, &screen_descriptor);
        }

        // -------------------------------------------------------------------
        // 6. Submit, present, free textures
        // -------------------------------------------------------------------
        self.queue.submit(Some(encoder.finish()));
        frame.present();

        for id in &full_output.textures_delta.free {
            self.egui_renderer.free_texture(id);
        }

        if should_reconfigure {
            RenderResult::Reconfigure
        } else {
            RenderResult::Drawn
        }
    }

    // -------------------------------------------------------------------
    // Snapshot reloading
    // -------------------------------------------------------------------

    /// Try to load a new snapshot from disk and replace current GPU resources
    /// atomically. On failure, the old snapshot/placeholder is preserved.
    pub(crate) fn apply_args(
        &mut self,
        new_args: ViewerArgs,
    ) -> Result<SnapshotInfo, Box<dyn Error>> {
        let (snapshot, info) = try_load_snapshot_resources(
            &self.device,
            &self.queue,
            &self.bind_group_layout,
            &new_args,
            self.config.width,
            self.config.height,
        )?;

        self.snapshot = snapshot;
        self.loaded_info = Some(info.clone());
        self.args = new_args;

        self.window.set_title(&format!(
            "MARL Viewer - tick {} species {} view {:?} cells {:?} ({} cells, {}x{}x{})",
            info.tick,
            info.species,
            info.view_mode,
            info.cell_mode,
            info.cell_count,
            info.grid[0],
            info.grid[1],
            info.grid[2]
        ));

        Ok(info)
    }

    /// Load a directory: discover ticks, choose initial tick, load snapshot.
    pub(crate) fn load_directory_from_gui(&mut self, dir: PathBuf) -> Result<(), Box<dyn Error>> {
        let ticks = discover_field_ticks(&dir)?;
        if ticks.is_empty() {
            return Err(format!("no tick_*.field.bin snapshots found in {}", dir.display()).into());
        }

        let tick = choose_initial_tick(self.args.tick, &ticks).ok_or("no valid tick found")?;

        // Validate run_meta.json exists and optionally clamp species
        let mut species = self.args.species;
        if let Ok(meta) = load_run_meta(&dir)
            && species >= meta.s_ext
        {
            species = meta.s_ext.saturating_sub(1);
        }

        let new_args = ViewerArgs {
            output_dir: dir,
            tick,
            species,
            ..self.args.clone()
        };

        let info = self.apply_args(new_args)?;
        self.gui.sync_loaded(&info, &self.args, ticks);
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// RenderResult
// ---------------------------------------------------------------------------

pub enum RenderResult {
    Drawn,
    Skip,
    Reconfigure,
}

// ---------------------------------------------------------------------------
// Snapshot GPU resource construction
// ---------------------------------------------------------------------------

/// Build `ViewerParams` from a loaded snapshot and viewer args.
fn build_viewer_params(
    payload: &SnapshotPayload,
    args: &ViewerArgs,
    width: u32,
    height: u32,
) -> Result<ViewerParams, Box<dyn Error>> {
    let cam: CameraBasis = camera_basis(args.view_mode);

    let max_dim = payload
        .meta
        .grid_x
        .max(payload.meta.grid_y)
        .max(payload.meta.grid_z) as f32;
    if max_dim <= 0.0 {
        return Err("zero max grid dimension".into());
    }
    let axis_scale = [
        payload.meta.grid_x as f32 / max_dim,
        payload.meta.grid_y as f32 / max_dim,
        payload.meta.grid_z as f32 / max_dim,
        0.0,
    ];

    let cells_enabled = if payload.cell_mode == CellMode::Off {
        0u32
    } else {
        1u32
    };

    Ok(ViewerParams {
        grid: [
            payload.meta.grid_x,
            payload.meta.grid_y,
            payload.meta.grid_z,
            payload.meta.s_ext,
        ],
        render: [width, height, payload.species, payload.steps],
        transfer: [
            payload.exposure,
            payload.density_scale,
            payload.cell_alpha,
            0.0,
        ],
        axis_scale,
        cam_right: [cam.right[0], cam.right[1], cam.right[2], cam.zoom],
        cam_up: [cam.up[0], cam.up[1], cam.up[2], 0.0],
        cam_dir: [cam.dir[0], cam.dir[1], cam.dir[2], 0.0],
        options: [cells_enabled, 0, 0, 0],
    })
}

/// Create a bind group using the shared layout and given resources.
fn create_snapshot_bind_group(
    device: &wgpu::Device,
    bind_group_layout: &wgpu::BindGroupLayout,
    field_view: &wgpu::TextureView,
    params_buffer: &wgpu::Buffer,
    cell_view: &wgpu::TextureView,
) -> wgpu::BindGroup {
    device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("MARL Viewer Bind Group"),
        layout: bind_group_layout,
        entries: &[
            wgpu::BindGroupEntry {
                binding: 0,
                resource: wgpu::BindingResource::TextureView(field_view),
            },
            wgpu::BindGroupEntry {
                binding: 1,
                resource: params_buffer.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 2,
                resource: wgpu::BindingResource::TextureView(cell_view),
            },
        ],
    })
}

/// Load snapshot from disk and build all GPU resources.
/// Returns both the GPU resources and display metadata, or an error.
fn try_load_snapshot_resources(
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    bind_group_layout: &wgpu::BindGroupLayout,
    args: &ViewerArgs,
    width: u32,
    height: u32,
) -> Result<(SnapshotGpuResources, SnapshotInfo), Box<dyn Error>> {
    let payload = load_snapshot(args)?;

    let params = build_viewer_params(&payload, args, width, height)?;

    // Field texture
    let (field_texture, field_view) = create_field_texture(device, queue, &payload)?;

    // Cell texture
    let cells_enabled = params.options[0] != 0;
    let (cell_texture, cell_view) = if cells_enabled {
        let (tex, view) = create_cell_texture(device, queue, &payload, args)?;
        (tex, view)
    } else {
        let (tex, view) = create_empty_cell_texture(device, queue, &payload)?;
        (tex, view)
    };

    // Params buffer
    let params_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("MARL Viewer Params"),
        contents: bytemuck::bytes_of(&params),
        usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
    });

    // Bind group
    let bind_group = create_snapshot_bind_group(
        device,
        bind_group_layout,
        &field_view,
        &params_buffer,
        &cell_view,
    );

    let info = SnapshotInfo {
        output_dir: args.output_dir.clone(),
        tick: payload.tick,
        species: payload.species,
        view_mode: args.view_mode,
        cell_mode: args.cell_mode,
        cell_count: payload.cells.len(),
        field_bytes: payload.field_bytes.len(),
        grid: [
            payload.meta.grid_x,
            payload.meta.grid_y,
            payload.meta.grid_z,
        ],
        s_ext: payload.meta.s_ext,
    };

    Ok((
        SnapshotGpuResources {
            bind_group,
            params_buffer,
            params,
            _field_texture: field_texture,
            _field_view: field_view,
            _cell_texture: cell_texture,
            _cell_view: cell_view,
        },
        info,
    ))
}

// ---------------------------------------------------------------------------
// Placeholder resources (1×1×1 zero textures)
// ---------------------------------------------------------------------------

/// Create a minimal 1×1×1 placeholder `SnapshotPayload`.
fn placeholder_payload(args: &ViewerArgs) -> SnapshotPayload {
    use marl_format::RunMeta;

    let meta = RunMeta::new(1, 1, 1, 1, 0, true, false);
    let field_bytes = vec![0u8; 4]; // 1×1×1×1×4 bytes
    SnapshotPayload {
        meta,
        field_bytes,
        cells: Vec::new(),
        tick: args.tick,
        species: 0,
        exposure: args.exposure,
        density_scale: args.density_scale,
        steps: args.steps,
        cell_mode: CellMode::Off,
        cell_alpha: args.cell_alpha,
    }
}

/// Build SnapshotGpuResources from a placeholder payload.
fn create_placeholder_resources(
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    bind_group_layout: &wgpu::BindGroupLayout,
    args: &ViewerArgs,
    width: u32,
    height: u32,
) -> Result<SnapshotGpuResources, Box<dyn Error>> {
    let payload = placeholder_payload(args);
    let params = build_viewer_params(&payload, args, width, height)?;
    let (field_texture, field_view) = create_field_texture(device, queue, &payload)?;
    let (cell_texture, cell_view) = create_empty_cell_texture(device, queue, &payload)?;
    let params_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("MARL Viewer Params (placeholder)"),
        contents: bytemuck::bytes_of(&params),
        usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
    });
    let bind_group = create_snapshot_bind_group(
        device,
        bind_group_layout,
        &field_view,
        &params_buffer,
        &cell_view,
    );
    Ok(SnapshotGpuResources {
        bind_group,
        params_buffer,
        params,
        _field_texture: field_texture,
        _field_view: field_view,
        _cell_texture: cell_texture,
        _cell_view: cell_view,
    })
}

// ---------------------------------------------------------------------------
// Field texture creation (preserved from Phase 1)
// ---------------------------------------------------------------------------

pub(crate) fn create_field_texture(
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    payload: &SnapshotPayload,
) -> Result<(wgpu::Texture, wgpu::TextureView), Box<dyn Error>> {
    let texture_width = payload
        .meta
        .grid_x
        .checked_mul(payload.meta.s_ext)
        .ok_or("field texture width overflows u32")?;
    let max_3d = device.limits().max_texture_dimension_3d;
    if texture_width > max_3d || payload.meta.grid_y > max_3d || payload.meta.grid_z > max_3d {
        return Err(format!(
            "field texture dimensions {}x{}x{} exceed this adapter's 3D texture limit of {}",
            texture_width, payload.meta.grid_y, payload.meta.grid_z, max_3d
        )
        .into());
    }
    let size = wgpu::Extent3d {
        width: texture_width,
        height: payload.meta.grid_y,
        depth_or_array_layers: payload.meta.grid_z,
    };
    let texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("MARL Viewer Field Texture"),
        size,
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D3,
        format: wgpu::TextureFormat::R32Float,
        usage: wgpu::TextureUsages::COPY_DST | wgpu::TextureUsages::TEXTURE_BINDING,
        view_formats: &[],
    });

    queue.write_texture(
        wgpu::TexelCopyTextureInfo {
            texture: &texture,
            mip_level: 0,
            origin: wgpu::Origin3d::ZERO,
            aspect: wgpu::TextureAspect::All,
        },
        &payload.field_bytes,
        wgpu::TexelCopyBufferLayout {
            offset: 0,
            bytes_per_row: Some(texture_width * 4),
            rows_per_image: Some(payload.meta.grid_y),
        },
        size,
    );

    let view = texture.create_view(&wgpu::TextureViewDescriptor {
        label: Some("MARL Viewer Field Texture View"),
        dimension: Some(wgpu::TextureViewDimension::D3),
        ..Default::default()
    });
    Ok((texture, view))
}

// ---------------------------------------------------------------------------
// Cell texture creation
// ---------------------------------------------------------------------------

/// Build a 3D `Rgba8Uint` texture sized `(grid_x, grid_y, grid_z)`.
///
/// Each occupied voxel stores RGBA bytes encoding the cell marker color and
/// alpha. Unoccupied voxels are zeroed.
fn create_cell_texture(
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    payload: &SnapshotPayload,
    args: &ViewerArgs,
) -> Result<(wgpu::Texture, wgpu::TextureView), Box<dyn Error>> {
    let gx = payload.meta.grid_x as usize;
    let gy = payload.meta.grid_y as usize;
    let gz = payload.meta.grid_z as usize;

    let max_3d = device.limits().max_texture_dimension_3d;
    if gx as u32 > max_3d || gy as u32 > max_3d || gz as u32 > max_3d {
        return Err(format!(
            "cell texture dimensions {}x{}x{} exceed adapter 3D limit {}",
            gx, gy, gz, max_3d
        )
        .into());
    }

    let size = wgpu::Extent3d {
        width: gx as u32,
        height: gy as u32,
        depth_or_array_layers: gz as u32,
    };

    let texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("MARL Viewer Cell Texture"),
        size,
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D3,
        format: wgpu::TextureFormat::Rgba8Uint,
        usage: wgpu::TextureUsages::COPY_DST | wgpu::TextureUsages::TEXTURE_BINDING,
        view_formats: &[],
    });

    let alpha_byte = (args.cell_alpha * 255.0).round() as u8;
    checked_cell_texture_len(gx, gy, gz)?;
    let data = build_cell_texture_data(gx, gy, gz, &payload.cells, payload.cell_mode, alpha_byte);

    let bytes_per_row = (gx as u32)
        .checked_mul(4)
        .ok_or("cell texture bytes_per_row overflow")?;

    queue.write_texture(
        wgpu::TexelCopyTextureInfo {
            texture: &texture,
            mip_level: 0,
            origin: wgpu::Origin3d::ZERO,
            aspect: wgpu::TextureAspect::All,
        },
        &data,
        wgpu::TexelCopyBufferLayout {
            offset: 0,
            bytes_per_row: Some(bytes_per_row),
            rows_per_image: Some(gy as u32),
        },
        size,
    );

    let view = texture.create_view(&wgpu::TextureViewDescriptor {
        label: Some("MARL Viewer Cell Texture View"),
        dimension: Some(wgpu::TextureViewDimension::D3),
        ..Default::default()
    });

    Ok((texture, view))
}

/// Create a minimal all-zeros cell texture for when cells are disabled.
fn create_empty_cell_texture(
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    payload: &SnapshotPayload,
) -> Result<(wgpu::Texture, wgpu::TextureView), Box<dyn Error>> {
    let gx = payload.meta.grid_x;
    let gy = payload.meta.grid_y;
    let gz = payload.meta.grid_z;

    let size = wgpu::Extent3d {
        width: gx,
        height: gy,
        depth_or_array_layers: gz,
    };
    let texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("MARL Viewer Empty Cell Texture"),
        size,
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D3,
        format: wgpu::TextureFormat::Rgba8Uint,
        usage: wgpu::TextureUsages::COPY_DST | wgpu::TextureUsages::TEXTURE_BINDING,
        view_formats: &[],
    });

    let data = vec![0u8; checked_cell_texture_len(gx as usize, gy as usize, gz as usize)?];
    let bytes_per_row = gx
        .checked_mul(4)
        .ok_or("empty cell texture bytes_per_row overflow")?;
    queue.write_texture(
        wgpu::TexelCopyTextureInfo {
            texture: &texture,
            mip_level: 0,
            origin: wgpu::Origin3d::ZERO,
            aspect: wgpu::TextureAspect::All,
        },
        &data,
        wgpu::TexelCopyBufferLayout {
            offset: 0,
            bytes_per_row: Some(bytes_per_row),
            rows_per_image: Some(gy),
        },
        size,
    );

    let view = texture.create_view(&wgpu::TextureViewDescriptor {
        label: Some("MARL Viewer Empty Cell Texture View"),
        dimension: Some(wgpu::TextureViewDimension::D3),
        ..Default::default()
    });
    Ok((texture, view))
}

// ---------------------------------------------------------------------------
// Cell texture data builder (pure, testable)
// ---------------------------------------------------------------------------

pub(crate) fn build_cell_texture_data(
    gx: usize,
    gy: usize,
    gz: usize,
    cells: &[LoadedCell],
    mode: CellMode,
    alpha: u8,
) -> Vec<u8> {
    let mut data = vec![
        0u8;
        checked_cell_texture_len(gx, gy, gz)
            .expect("cell texture dimensions overflow addressable memory")
    ];

    match mode {
        CellMode::Off => {
            // already all zeros
        }
        CellMode::Starter => {
            for cell in cells {
                let x = cell.pos[0] as usize;
                let y = cell.pos[1] as usize;
                let z = cell.pos[2] as usize;
                if x >= gx || y >= gy || z >= gz {
                    continue;
                }
                let idx = ((z * gy + y) * gx + x) * 4;
                // If already populated (duplicate), keep first; log warning
                if data[idx + 3] != 0 {
                    // duplicate: keep the first, skip this one
                    continue;
                }
                let (r, g, b) = starter_color(cell.starter_type);
                data[idx] = r;
                data[idx + 1] = g;
                data[idx + 2] = b;
                data[idx + 3] = alpha;
            }
        }
        CellMode::Energy => {
            let max_energy = cells
                .iter()
                .map(|c| c.energy)
                .fold(0.0f32, |a, b| a.max(b))
                .max(1.0); // safe fallback
            for cell in cells {
                let x = cell.pos[0] as usize;
                let y = cell.pos[1] as usize;
                let z = cell.pos[2] as usize;
                if x >= gx || y >= gy || z >= gz {
                    continue;
                }
                let idx = ((z * gy + y) * gx + x) * 4;
                // For energy mode, keep the higher-energy cell on duplicate
                if data[idx + 3] != 0 {
                    continue; // first wins for simplicity
                }
                let t = (cell.energy / max_energy).clamp(0.0, 1.0);
                let (r, g, b) = energy_color(t);
                data[idx] = r;
                data[idx + 1] = g;
                data[idx + 2] = b;
                data[idx + 3] = alpha;
            }
        }
    }

    data
}

fn checked_cell_texture_len(gx: usize, gy: usize, gz: usize) -> Result<usize, Box<dyn Error>> {
    gx.checked_mul(gy)
        .and_then(|v| v.checked_mul(gz))
        .and_then(|v| v.checked_mul(4))
        .ok_or_else(|| {
            format!("cell texture dimensions {gx}x{gy}x{gz} overflow addressable memory").into()
        })
}

fn starter_color(starter_type: u8) -> (u8, u8, u8) {
    match starter_type {
        0 => (230, 60, 55),  // phototroph: red
        1 => (70, 220, 90),  // chemolithotroph: green
        2 => (75, 130, 255), // anaerobe: blue
        _ => (230, 75, 230), // unknown: magenta
    }
}

fn energy_color(t: f32) -> (u8, u8, u8) {
    // Dark purple (low energy) → yellow (high energy)
    let t = t.clamp(0.0, 1.0);
    let r = (30.0 + t * 225.0).round() as u8;
    let g = (10.0 + t * 245.0).round() as u8;
    let b = (60.0 * (1.0 - t)).round() as u8;
    (r, g, b)
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn viewer_params_size_multiple_of_16() {
        assert_eq!(std::mem::size_of::<ViewerParams>() % 16, 0);
    }

    #[test]
    fn empty_cell_data_all_zeros() {
        let data = build_cell_texture_data(2, 2, 2, &[], CellMode::Starter, 200);
        assert_eq!(data.len(), 2 * 2 * 2 * 4);
        assert!(data.iter().all(|&b| b == 0));
    }

    #[test]
    fn starter_voxel_at_origin() {
        let cell = LoadedCell {
            pos: [0, 0, 0],
            lineage_id: 1,
            starter_type: 0,
            energy: 1.0,
        };
        let data = build_cell_texture_data(2, 2, 2, &[cell], CellMode::Starter, 200);
        //  4 bytes per voxel = RGBA
        assert_eq!(data[0], 230); // phototroph red
        assert_eq!(data[1], 60);
        assert_eq!(data[2], 55);
        assert_eq!(data[3], 200); // alpha
        // Rest zero
        for (i, byte) in data.iter().enumerate().skip(4) {
            assert_eq!(*byte, 0, "byte {i} should be zero");
        }
    }

    #[test]
    fn starter_voxel_different_positions() {
        let cells = vec![
            LoadedCell {
                pos: [0, 0, 0],
                lineage_id: 1,
                starter_type: 1,
                energy: 1.0,
            },
            LoadedCell {
                pos: [1, 0, 0],
                lineage_id: 2,
                starter_type: 2,
                energy: 1.0,
            },
        ];
        let data = build_cell_texture_data(2, 1, 1, &cells, CellMode::Starter, 255);
        // cell 0 at idx 0: green
        assert_eq!(&data[0..4], &[70, 220, 90, 255]);
        // cell 1 at idx 4: blue
        assert_eq!(&data[4..8], &[75, 130, 255, 255]);
    }

    #[test]
    fn energy_voxel_color() {
        let cell = LoadedCell {
            pos: [0, 0, 0],
            lineage_id: 1,
            starter_type: 0,
            energy: 5.0,
        };
        let cells = vec![cell];
        let data = build_cell_texture_data(2, 2, 2, &cells, CellMode::Energy, 128);
        // max_energy = 5.0, so t = 1.0 => full energy color
        assert_eq!(data[3], 128); // alpha
        // R,G,B should be near (255, 255, 0) for t=1
        assert!(data[0] > 200);
        assert!(data[1] > 200);
        assert!(data[2] < 30);
    }

    #[test]
    fn duplicate_position_first_wins_starter() {
        let cells = vec![
            LoadedCell {
                pos: [0, 0, 0],
                lineage_id: 1,
                starter_type: 0,
                energy: 1.0,
            },
            LoadedCell {
                pos: [0, 0, 0],
                lineage_id: 2,
                starter_type: 1,
                energy: 2.0,
            },
        ];
        let data = build_cell_texture_data(2, 2, 2, &cells, CellMode::Starter, 200);
        // First cell (phototroph) should win
        assert_eq!(data[0], 230);
        assert_eq!(data[1], 60);
        assert_eq!(data[2], 55);
    }

    #[test]
    fn out_of_bounds_pos_ignored() {
        let cell = LoadedCell {
            pos: [2, 0, 0], // x >= gx=2
            lineage_id: 1,
            starter_type: 0,
            energy: 1.0,
        };
        let data = build_cell_texture_data(2, 2, 2, &[cell], CellMode::Starter, 200);
        assert!(data.iter().all(|&b| b == 0));
    }

    #[test]
    fn cell_mode_off_always_empty() {
        let cell = LoadedCell {
            pos: [0, 0, 0],
            lineage_id: 1,
            starter_type: 0,
            energy: 1.0,
        };
        let data = build_cell_texture_data(2, 2, 2, &[cell], CellMode::Off, 200);
        assert!(data.iter().all(|&b| b == 0));
    }

    #[test]
    fn checked_cell_texture_len_rejects_overflow() {
        let err = checked_cell_texture_len(usize::MAX, 2, 1).unwrap_err();
        assert!(err.to_string().contains("overflow"));
    }

    #[test]
    fn placeholder_payload_dimensions() {
        let args = ViewerArgs {
            output_dir: std::path::PathBuf::from("/nonexistent"),
            tick: 0,
            species: 1,
            exposure: 18.0,
            density_scale: 2.0,
            steps: 160,
            view_mode: ViewMode::Iso,
            cell_mode: CellMode::Starter,
            cell_alpha: 0.95,
        };
        let payload = placeholder_payload(&args);
        assert_eq!(payload.meta.grid_x, 1);
        assert_eq!(payload.meta.grid_y, 1);
        assert_eq!(payload.meta.grid_z, 1);
        assert_eq!(payload.meta.s_ext, 1);
        assert_eq!(payload.field_bytes.len(), 4); // 1×1×1×1 × 4 bytes
        assert!(payload.cells.is_empty());
    }
}
