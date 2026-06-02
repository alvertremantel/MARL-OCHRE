//! GUI state and action helpers for the MARL viewer.
//!
//! `GuiState` tracks the draft form fields, available ticks, and status
//! messages.  The `show()` method draws the toolbar and sidebar, returning
//! `GuiAction` values for the renderer to process.

use std::path::PathBuf;

use marl_viewer_core::args::{CellMode, ViewMode, ViewerArgs};
use marl_viewer_core::types::{GuiAction, SnapshotInfo};

// ---------------------------------------------------------------------------
// GuiState
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub(crate) struct GuiState {
    pub(crate) directory_text: String,
    pub(crate) tick_text: String,
    pub(crate) available_ticks: Vec<u64>,
    // Draft view settings (may differ from loaded args)
    pub(crate) draft_species: String,
    pub(crate) draft_view_mode: ViewMode,
    pub(crate) draft_cell_mode: CellMode,
    pub(crate) draft_cell_alpha: String,
    pub(crate) draft_density_scale: String,
    pub(crate) draft_exposure: String,
    pub(crate) draft_steps: String,
    pub(crate) status: String,
    pub(crate) status_is_error: bool,
}

impl GuiState {
    pub(crate) fn new(args: &ViewerArgs) -> Self {
        Self {
            directory_text: args.output_dir.to_str().unwrap_or("").to_string(),
            tick_text: args.tick.to_string(),
            available_ticks: Vec::new(),
            draft_species: args.species.to_string(),
            draft_view_mode: args.view_mode,
            draft_cell_mode: args.cell_mode,
            draft_cell_alpha: format!("{:.2}", args.cell_alpha),
            draft_density_scale: format!("{:.2}", args.density_scale),
            draft_exposure: format!("{:.1}", args.exposure),
            draft_steps: args.steps.to_string(),
            status: String::new(),
            status_is_error: false,
        }
    }

    /// Update display fields to match a successfully loaded snapshot.
    pub(crate) fn sync_loaded(
        &mut self,
        info: &SnapshotInfo,
        args: &ViewerArgs,
        available_ticks: Vec<u64>,
    ) {
        self.directory_text = info.output_dir.to_str().unwrap_or("").to_string();
        self.tick_text = info.tick.to_string();
        self.available_ticks = available_ticks;
        self.draft_species = args.species.to_string();
        self.draft_view_mode = args.view_mode;
        self.draft_cell_mode = args.cell_mode;
        self.draft_cell_alpha = format!("{:.2}", args.cell_alpha);
        self.draft_density_scale = format!("{:.2}", args.density_scale);
        self.draft_exposure = format!("{:.1}", args.exposure);
        self.draft_steps = args.steps.to_string();
        self.clear_status();
    }

    pub(crate) fn set_error(&mut self, msg: impl Into<String>) {
        self.status = msg.into();
        self.status_is_error = true;
    }

    pub(crate) fn set_info(&mut self, msg: impl Into<String>) {
        self.status = msg.into();
        self.status_is_error = false;
    }

    pub(crate) fn clear_status(&mut self) {
        self.status.clear();
        self.status_is_error = false;
    }

    /// Parse draft view settings into a new `ViewerArgs`, keeping the
    /// `output_dir` and `tick` from the provided base. Returns an error
    /// message string on parse/validation failure.
    pub(crate) fn build_view_args_from_drafts(
        &self,
        base: &ViewerArgs,
        s_ext: Option<u32>,
    ) -> Result<ViewerArgs, String> {
        let species: u32 = self
            .draft_species
            .parse()
            .map_err(|_| format!("invalid species number: {}", self.draft_species))?;
        if let Some(max_ext) = s_ext
            && species >= max_ext
        {
            return Err(format!(
                "species {species} is out of range for {max_ext} external species"
            ));
        }

        let exposure: f32 = self
            .draft_exposure
            .parse()
            .map_err(|_| format!("invalid exposure: {}", self.draft_exposure))?;
        if !exposure.is_finite() || exposure <= 0.0 {
            return Err(format!(
                "exposure must be a positive finite number, got {}",
                self.draft_exposure
            ));
        }

        let density_scale: f32 = self
            .draft_density_scale
            .parse()
            .map_err(|_| format!("invalid scale: {}", self.draft_density_scale))?;
        if !density_scale.is_finite() || density_scale <= 0.0 {
            return Err(format!(
                "scale must be a positive finite number, got {}",
                self.draft_density_scale
            ));
        }

        let steps: u32 = self
            .draft_steps
            .parse()
            .map_err(|_| format!("invalid steps: {}", self.draft_steps))?;
        if steps == 0 {
            return Err(format!(
                "steps must be greater than zero, got {}",
                self.draft_steps
            ));
        }

        let cell_alpha: f32 = self
            .draft_cell_alpha
            .parse()
            .map_err(|_| format!("invalid cell alpha: {}", self.draft_cell_alpha))?;
        if !cell_alpha.is_finite() || cell_alpha <= 0.0 || cell_alpha > 1.0 {
            return Err(format!(
                "cell alpha must be in (0, 1], got {}",
                self.draft_cell_alpha
            ));
        }

        Ok(ViewerArgs {
            output_dir: base.output_dir.clone(),
            tick: base.tick,
            species,
            exposure,
            density_scale,
            steps,
            view_mode: self.draft_view_mode,
            cell_mode: self.draft_cell_mode,
            cell_alpha,
        })
    }

    /// Reset all draft fields to match the currently loaded args.
    pub(crate) fn reset_drafts_from_args(&mut self, args: &ViewerArgs) {
        self.draft_species = args.species.to_string();
        self.draft_view_mode = args.view_mode;
        self.draft_cell_mode = args.cell_mode;
        self.draft_cell_alpha = format!("{:.2}", args.cell_alpha);
        self.draft_density_scale = format!("{:.2}", args.density_scale);
        self.draft_exposure = format!("{:.1}", args.exposure);
        self.draft_steps = args.steps.to_string();
    }
}

// ---------------------------------------------------------------------------
// GUI drawing
// ---------------------------------------------------------------------------

impl GuiState {
    /// Draw the top toolbar and return any user actions.
    pub(crate) fn show(
        &mut self,
        ctx: &egui::Context,
        loaded: Option<&SnapshotInfo>,
        _args: &ViewerArgs,
    ) -> Vec<GuiAction> {
        let mut actions = Vec::new();

        // Top toolbar — drawn directly at root level (no CentralPanel wrapper)
        #[allow(deprecated)]
        egui::TopBottomPanel::top("marl_toolbar").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label("Directory:");
                if ui
                    .add(
                        egui::TextEdit::singleline(&mut self.directory_text)
                            .hint_text("path/to/output")
                            .desired_width(300.0),
                    )
                    .lost_focus()
                    && ui.input(|i| i.key_pressed(egui::Key::Enter))
                {
                    actions.push(GuiAction::LoadDirectory(PathBuf::from(
                        self.directory_text.clone(),
                    )));
                }
                if ui.button("Open…").clicked() {
                    actions.push(GuiAction::OpenDirectoryDialog);
                }
                if ui.button("Load Dir").clicked() {
                    actions.push(GuiAction::LoadDirectory(PathBuf::from(
                        self.directory_text.clone(),
                    )));
                }

                ui.separator();

                ui.label("Tick:");
                if ui
                    .add(egui::TextEdit::singleline(&mut self.tick_text).desired_width(60.0))
                    .lost_focus()
                    && ui.input(|i| i.key_pressed(egui::Key::Enter))
                    && let Ok(tick) = self.tick_text.parse::<u64>()
                {
                    actions.push(GuiAction::LoadTick(tick));
                }
                if ui.button("Go").clicked()
                    && let Ok(tick) = self.tick_text.parse::<u64>()
                {
                    actions.push(GuiAction::LoadTick(tick));
                }

                // Navigation buttons
                ui.add_enabled_ui(!self.available_ticks.is_empty(), |ui| {
                    if ui.button("|<").clicked() {
                        actions.push(GuiAction::FirstTick);
                    }
                    if ui.button("<").clicked() {
                        actions.push(GuiAction::PrevTick);
                    }
                    if ui.button(">").clicked() {
                        actions.push(GuiAction::NextTick);
                    }
                    if ui.button(">|").clicked() {
                        actions.push(GuiAction::LastTick);
                    }
                });
                if ui.button("Reload").clicked() {
                    actions.push(GuiAction::ReloadCurrent);
                }

                // Summary of loaded state
                if let Some(info) = loaded {
                    ui.label(format!(
                        "  |  {}×{}×{}  s_ext={}  {} cells",
                        info.grid[0], info.grid[1], info.grid[2], info.s_ext, info.cell_count,
                    ));
                    if !self.available_ticks.is_empty() {
                        ui.label(format!(
                            "  ticks: {}",
                            tick_list_summary(&self.available_ticks)
                        ));
                    }
                } else {
                    ui.label("  |  No snapshot loaded");
                }
            });
            // Status bar (second row inside toolbar panel)
            if !self.status.is_empty() {
                ui.horizontal(|ui| {
                    if self.status_is_error {
                        ui.colored_label(egui::Color32::RED, &self.status);
                    } else {
                        ui.label(&self.status);
                    }
                });
            }
        });

        // Left view settings panel
        #[allow(deprecated)]
        egui::SidePanel::left("marl_view_settings")
            .resizable(false)
            .default_width(220.0)
            .show(ctx, |ui| {
                egui::CollapsingHeader::new("View Settings")
                    .default_open(false)
                    .show(ui, |ui| {
                        ui.horizontal(|ui| {
                            ui.label("Species:");
                            ui.add(
                                egui::TextEdit::singleline(&mut self.draft_species)
                                    .desired_width(40.0),
                            );
                        });
                        ui.horizontal(|ui| {
                            ui.label("View:");
                            egui::ComboBox::from_id_salt("view_mode")
                                .selected_text(self.draft_view_mode.as_str())
                                .show_ui(ui, |ui| {
                                    for mode in ViewMode::all() {
                                        ui.selectable_value(
                                            &mut self.draft_view_mode,
                                            mode,
                                            mode.as_str(),
                                        );
                                    }
                                });
                        });
                        ui.horizontal(|ui| {
                            ui.label("Cells:");
                            egui::ComboBox::from_id_salt("cell_mode")
                                .selected_text(self.draft_cell_mode.as_str())
                                .show_ui(ui, |ui| {
                                    for mode in CellMode::all() {
                                        ui.selectable_value(
                                            &mut self.draft_cell_mode,
                                            mode,
                                            mode.as_str(),
                                        );
                                    }
                                });
                        });
                        ui.horizontal(|ui| {
                            ui.label("Alpha:");
                            ui.add(
                                egui::TextEdit::singleline(&mut self.draft_cell_alpha)
                                    .desired_width(50.0),
                            );
                        });
                        ui.horizontal(|ui| {
                            ui.label("Scale:");
                            ui.add(
                                egui::TextEdit::singleline(&mut self.draft_density_scale)
                                    .desired_width(50.0),
                            );
                        });
                        ui.horizontal(|ui| {
                            ui.label("Exposure:");
                            ui.add(
                                egui::TextEdit::singleline(&mut self.draft_exposure)
                                    .desired_width(50.0),
                            );
                        });
                        ui.horizontal(|ui| {
                            ui.label("Steps:");
                            ui.add(
                                egui::TextEdit::singleline(&mut self.draft_steps)
                                    .desired_width(50.0),
                            );
                        });
                        ui.horizontal(|ui| {
                            if ui.button("Apply").clicked() {
                                actions.push(GuiAction::ApplyViewSettings);
                            }
                            if ui.button("Reset").clicked() {
                                actions.push(GuiAction::ResetDraftFromLoaded);
                            }
                        });
                    });
            });

        actions
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Make a compact summary of the available tick list.
fn tick_list_summary(ticks: &[u64]) -> String {
    if ticks.is_empty() {
        return "(none)".to_string();
    }
    if ticks.len() <= 5 {
        let items: Vec<String> = ticks.iter().map(|t| t.to_string()).collect();
        return items.join(", ");
    }
    format!(
        "{} … {}  ({} total)",
        ticks.first().unwrap(),
        ticks.last().unwrap(),
        ticks.len()
    )
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tick_list_summary_small() {
        assert_eq!(tick_list_summary(&[0, 500]), "0, 500");
    }

    #[test]
    fn tick_list_summary_large() {
        let ticks: Vec<u64> = (0..10).collect();
        let s = tick_list_summary(&ticks);
        assert!(s.contains("…"));
        assert!(s.contains("total"));
    }

    #[test]
    fn tick_list_summary_empty() {
        assert_eq!(tick_list_summary(&[]), "(none)");
    }
}
