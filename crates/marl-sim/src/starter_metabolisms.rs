use marl_cell::cell::{
    CellState, EffectorParams, FateParams, Reaction, ReceptorParams, Ruleset, TransportParams,
};
use marl_config::{M_INT, R_MAX, S_EFFECTORS, S_RECEPTORS, S_TRANSPORTERS};

fn inactive_receptor() -> ReceptorParams {
    ReceptorParams {
        k_half: 1.0,
        n_hill: 2.0,
        gain: 0.0,
    }
}

fn inactive_transport() -> TransportParams {
    TransportParams {
        uptake_rate: 0.0,
        secrete_rate: 0.0,
        ext_species: 0,
        int_species: 0,
        gate_receptor: 0,
        gate_weight: 0.0,
    }
}

fn inactive_reaction() -> Reaction {
    Reaction {
        substrate: 0,
        product: 0,
        catalyst: 0,
        cofactor: 0xFF,
        k_m: 1.0,
        v_max: 0.0,
        k_cat: 0.5,
    }
}

fn inactive_effector() -> EffectorParams {
    EffectorParams {
        threshold: 10.0,
        rate: 0.0,
        int_species: 0,
        ext_species: 0,
    }
}

/// Phototroph: uses light + reductant to produce energy and oxidant.
/// Surface dweller. Spec Section 4.1.
pub fn make_phototroph(pos: [u16; 3], lineage_id: u64) -> CellState {
    let mut receptors: [ReceptorParams; S_RECEPTORS] = std::array::from_fn(|_| inactive_receptor());
    let mut transport: [TransportParams; S_TRANSPORTERS] =
        std::array::from_fn(|_| inactive_transport());
    let mut reactions: [Reaction; R_MAX] = std::array::from_fn(|_| inactive_reaction());
    let mut effectors: [EffectorParams; S_EFFECTORS] = std::array::from_fn(|_| inactive_effector());

    // Sense reductant and carbon
    receptors[2] = ReceptorParams {
        k_half: 0.3,
        n_hill: 2.0,
        gain: 1.0,
    };
    receptors[3] = ReceptorParams {
        k_half: 0.5,
        n_hill: 2.0,
        gain: 1.0,
    };

    // Transport: carbon is the primary fuel, secretes oxidant + organic waste
    transport[0] = TransportParams {
        uptake_rate: 0.6,
        secrete_rate: 0.0,
        ext_species: 3,
        int_species: 3,
        gate_receptor: 0,
        gate_weight: 0.0,
    }; // carbon in (primary)
    transport[1] = TransportParams {
        uptake_rate: 0.0,
        secrete_rate: 0.6,
        ext_species: 1,
        int_species: 1,
        gate_receptor: 0,
        gate_weight: 0.0,
    }; // oxidant out
    transport[2] = TransportParams {
        uptake_rate: 0.0,
        secrete_rate: 0.3,
        ext_species: 4,
        int_species: 4,
        gate_receptor: 0,
        gate_weight: 0.0,
    }; // organic waste out
    transport[3] = TransportParams {
        uptake_rate: 0.1,
        secrete_rate: 0.0,
        ext_species: 2,
        int_species: 2,
        gate_receptor: 0,
        gate_weight: 0.0,
    }; // some reductant in (secondary)

    // Rxn 0: carbon(3) -> energy(0), cat=LIGHT(15)  — photosynthesis: CO2 + light -> energy
    reactions[0] = Reaction {
        substrate: 3,
        product: 0,
        catalyst: 15,
        cofactor: 0xFF,
        k_m: 0.2,
        v_max: 0.8,
        k_cat: 0.1,
    };
    // Rxn 1: carbon(3) -> oxidant(1), cat=LIGHT(15)  — O2 production (water splitting analog)
    reactions[1] = Reaction {
        substrate: 3,
        product: 1,
        catalyst: 15,
        cofactor: 0xFF,
        k_m: 0.2,
        v_max: 0.4,
        k_cat: 0.1,
    };
    // Rxn 2: carbon(3) -> organic(4), cat=energy(0)  — carbon fixation into biomass
    reactions[2] = Reaction {
        substrate: 3,
        product: 4,
        catalyst: 0,
        cofactor: 0xFF,
        k_m: 0.3,
        v_max: 0.3,
        k_cat: 0.3,
    };
    // Rxn 3: carbon(3) -> enzyme-A(5), cat=enzyme-B(6)  — autocatalytic pair
    reactions[3] = Reaction {
        substrate: 3,
        product: 5,
        catalyst: 6,
        cofactor: 0xFF,
        k_m: 0.3,
        v_max: 0.15,
        k_cat: 0.1,
    };
    // Rxn 4: carbon(3) -> enzyme-B(6), cat=enzyme-A(5)
    reactions[4] = Reaction {
        substrate: 3,
        product: 6,
        catalyst: 5,
        cofactor: 0xFF,
        k_m: 0.3,
        v_max: 0.1,
        k_cat: 0.1,
    };

    // Effectors: secrete oxidant and organic waste
    effectors[0] = EffectorParams {
        threshold: 0.5,
        rate: 0.8,
        int_species: 1,
        ext_species: 1,
    }; // oxidant out
    effectors[1] = EffectorParams {
        threshold: 1.0,
        rate: 0.2,
        int_species: 4,
        ext_species: 4,
    }; // waste out

    let ruleset = Ruleset {
        receptors,
        transport,
        reactions,
        effectors,
        fate: FateParams {
            division_energy: 1.2,
            death_energy: 0.05,
            quiescence_energy: 0.15,
            division_prep_ticks: 20.0,
        },
        hgt_propensity: 0.1,
        mutation_rate: 0.05,
    };

    let mut internal = [0.0f32; M_INT];
    internal[0] = 0.3; // starting energy
    internal[3] = 0.3; // starting carbon (photosynthesis substrate)
    internal[5] = 0.05; // enzyme-A seed
    internal[6] = 0.05; // enzyme-B seed

    CellState {
        pos,
        lineage_id,
        age: 0,
        internal,
        ruleset,
        quiescent: false,
        starter_type: 0,
        prep_remaining: 0,
    }
}

/// Chemolithotroph: oxidizes reductant using oxidant at the chemocline.
/// Models sulfur-oxidizing bacteria at the oxic-anoxic interface.
pub fn make_chemolithotroph(pos: [u16; 3], lineage_id: u64) -> CellState {
    let mut receptors: [ReceptorParams; S_RECEPTORS] = std::array::from_fn(|_| inactive_receptor());
    let mut transport: [TransportParams; S_TRANSPORTERS] =
        std::array::from_fn(|_| inactive_transport());
    let mut reactions: [Reaction; R_MAX] = std::array::from_fn(|_| inactive_reaction());
    let mut effectors: [EffectorParams; S_EFFECTORS] = std::array::from_fn(|_| inactive_effector());

    receptors[1] = ReceptorParams {
        k_half: 0.3,
        n_hill: 2.0,
        gain: 1.0,
    }; // oxidant
    receptors[2] = ReceptorParams {
        k_half: 0.3,
        n_hill: 2.0,
        gain: 1.0,
    }; // reductant

    // Transport: take in both oxidant and reductant, secrete organic waste
    transport[0] = TransportParams {
        uptake_rate: 0.7,
        secrete_rate: 0.0,
        ext_species: 1,
        int_species: 1,
        gate_receptor: 0,
        gate_weight: 0.0,
    }; // oxidant in
    transport[1] = TransportParams {
        uptake_rate: 0.7,
        secrete_rate: 0.0,
        ext_species: 2,
        int_species: 2,
        gate_receptor: 0,
        gate_weight: 0.0,
    }; // reductant in
    transport[2] = TransportParams {
        uptake_rate: 0.2,
        secrete_rate: 0.0,
        ext_species: 3,
        int_species: 3,
        gate_receptor: 0,
        gate_weight: 0.0,
    }; // carbon in (for enzymes)
    transport[3] = TransportParams {
        uptake_rate: 0.0,
        secrete_rate: 0.3,
        ext_species: 4,
        int_species: 4,
        gate_receptor: 0,
        gate_weight: 0.0,
    }; // organic waste out

    // Rxn 0: reductant(2) -> energy(0), cat=enzyme-A(5), cofactor=oxidant(1)  — sulfur oxidation
    reactions[0] = Reaction {
        substrate: 2,
        product: 0,
        catalyst: 5,
        cofactor: 1,
        k_m: 0.15,
        v_max: 0.7,
        k_cat: 0.1,
    };
    // Rxn 1: oxidant(1) -> organic(4), cat=energy(0)  — oxidant processing byproduct
    reactions[1] = Reaction {
        substrate: 1,
        product: 4,
        catalyst: 0,
        cofactor: 0xFF,
        k_m: 0.3,
        v_max: 0.3,
        k_cat: 0.2,
    };
    // Rxn 2: carb_reserve(7) -> energy(0), cat=enzyme-A(5)  — slow burn of internal carbon store
    reactions[2] = Reaction {
        substrate: 7,
        product: 0,
        catalyst: 5,
        cofactor: 0xFF,
        k_m: 0.3,
        v_max: 0.15,
        k_cat: 0.1,
    };
    // Rxn 3-4: autocatalytic enzyme loop
    reactions[3] = Reaction {
        substrate: 3,
        product: 5,
        catalyst: 6,
        cofactor: 0xFF,
        k_m: 0.3,
        v_max: 0.15,
        k_cat: 0.1,
    };
    reactions[4] = Reaction {
        substrate: 3,
        product: 6,
        catalyst: 5,
        cofactor: 0xFF,
        k_m: 0.3,
        v_max: 0.1,
        k_cat: 0.1,
    };

    // Effectors: secrete organic waste
    effectors[0] = EffectorParams {
        threshold: 0.5,
        rate: 0.3,
        int_species: 4,
        ext_species: 4,
    };

    let ruleset = Ruleset {
        receptors,
        transport,
        reactions,
        effectors,
        fate: FateParams {
            division_energy: 1.0,
            death_energy: 0.05,
            quiescence_energy: 0.12,
            division_prep_ticks: 20.0,
        },
        hgt_propensity: 0.1,
        mutation_rate: 0.05,
    };

    let mut internal = [0.0f32; M_INT];
    internal[0] = 1.5; // substantial starting energy
    internal[1] = 0.5; // starting oxidant
    internal[2] = 0.5; // starting reductant
    internal[5] = 0.05;
    internal[6] = 0.05;
    internal[7] = 5.0; // carb reserve — slow-burn fuel while waiting for gradients to form

    CellState {
        pos,
        lineage_id,
        age: 0,
        internal,
        ruleset,
        quiescent: false,
        starter_type: 1,
        prep_remaining: 0,
    }
}

/// Anaerobe: uses reductant for energy. Killed by oxidant. Deep zone.
/// Spec Section 4.3.
pub fn make_anaerobe(pos: [u16; 3], lineage_id: u64) -> CellState {
    let mut receptors: [ReceptorParams; S_RECEPTORS] = std::array::from_fn(|_| inactive_receptor());
    let mut transport: [TransportParams; S_TRANSPORTERS] =
        std::array::from_fn(|_| inactive_transport());
    let mut reactions: [Reaction; R_MAX] = std::array::from_fn(|_| inactive_reaction());
    let mut effectors: [EffectorParams; S_EFFECTORS] = std::array::from_fn(|_| inactive_effector());

    receptors[2] = ReceptorParams {
        k_half: 0.3,
        n_hill: 2.0,
        gain: 1.0,
    }; // reductant

    // Transport — higher uptake rates to match the strong vent chemistry
    transport[0] = TransportParams {
        uptake_rate: 0.9,
        secrete_rate: 0.0,
        ext_species: 2,
        int_species: 2,
        gate_receptor: 0,
        gate_weight: 0.0,
    }; // reductant in (primary fuel)
    transport[1] = TransportParams {
        uptake_rate: 0.3,
        secrete_rate: 0.0,
        ext_species: 3,
        int_species: 3,
        gate_receptor: 0,
        gate_weight: 0.0,
    }; // carbon in
    transport[2] = TransportParams {
        uptake_rate: 0.0,
        secrete_rate: 0.5,
        ext_species: 4,
        int_species: 4,
        gate_receptor: 0,
        gate_weight: 0.0,
    }; // organic waste out
    transport[3] = TransportParams {
        uptake_rate: 0.1,
        secrete_rate: 0.0,
        ext_species: 1,
        int_species: 1,
        gate_receptor: 0,
        gate_weight: 0.0,
    }; // oxidant in (inadvertent!)

    // Rxn 0: reductant(2) -> energy(0), cat=enzyme-A(5)  — anaerobic respiration (BUFFED v_max)
    reactions[0] = Reaction {
        substrate: 2,
        product: 0,
        catalyst: 5,
        cofactor: 0xFF,
        k_m: 0.15,
        v_max: 0.8,
        k_cat: 0.1,
    };
    // Rxn 1: carbon(3) -> organic(4), cat=energy(0)  — fermentation
    reactions[1] = Reaction {
        substrate: 3,
        product: 4,
        catalyst: 0,
        cofactor: 0xFF,
        k_m: 0.2,
        v_max: 0.3,
        k_cat: 0.2,
    };
    // Rxn 2: OXIDANT TOXICITY — energy(0) -> carbon(3), cat=oxidant(1)
    //   High v_max + low k_m = even trace oxidant is lethal
    reactions[2] = Reaction {
        substrate: 0,
        product: 3,
        catalyst: 1,
        cofactor: 0xFF,
        k_m: 0.01,
        v_max: 2.0,
        k_cat: 0.01,
    };
    // Rxn 3-4: autocatalytic enzyme loop
    reactions[3] = Reaction {
        substrate: 3,
        product: 5,
        catalyst: 6,
        cofactor: 0xFF,
        k_m: 0.3,
        v_max: 0.15,
        k_cat: 0.1,
    };
    reactions[4] = Reaction {
        substrate: 3,
        product: 6,
        catalyst: 5,
        cofactor: 0xFF,
        k_m: 0.3,
        v_max: 0.1,
        k_cat: 0.1,
    };

    // Effectors: secrete organic waste
    effectors[0] = EffectorParams {
        threshold: 0.5,
        rate: 0.5,
        int_species: 4,
        ext_species: 4,
    };

    let ruleset = Ruleset {
        receptors,
        transport,
        reactions,
        effectors,
        fate: FateParams {
            division_energy: 0.8,
            death_energy: 0.05,
            quiescence_energy: 0.1,
            division_prep_ticks: 20.0,
        },
        hgt_propensity: 0.1,
        mutation_rate: 0.05,
    };

    let mut internal = [0.0f32; M_INT];
    internal[0] = 0.5; // more starting energy (needs to survive the 20-tick prep phase)
    internal[2] = 0.7; // more starting reductant (vents are chemically rich, help bootstrap)
    internal[5] = 0.05;
    internal[6] = 0.05;

    CellState {
        pos,
        lineage_id,
        age: 0,
        internal,
        ruleset,
        quiescent: false,
        starter_type: 2,
        prep_remaining: 0,
    }
}
