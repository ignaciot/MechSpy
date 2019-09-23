import open_tg_gates_canonical_mechs
import open_tg_gates_all_other_chemicals
import dixa_heparg_experiments
import dixa_hepg2_experiments
import dixa_kidney_experiments
import dixa_lung_experiments
import dixa078_HepG2_experiments
import tobacco_nasal_experiments
import tobacco_buccal_experiments
import tobacco_bronchial_experiments

EXPERIMENTS = { **open_tg_gates_canonical_mechs.EXPERIMENTS,
                **open_tg_gates_all_other_chemicals.EXPERIMENTS,
                **dixa_heparg_experiments.EXPERIMENTS,
                **dixa_hepg2_experiments.EXPERIMENTS,
                **dixa_kidney_experiments.EXPERIMENTS,
                **dixa_lung_experiments.EXPERIMENTS,
                **dixa078_HepG2_experiments.EXPERIMENTS,
                **tobacco_nasal_experiments.EXPERIMENTS,
                **tobacco_buccal_experiments.EXPERIMENTS,
                **tobacco_bronchial_experiments.EXPERIMENTS,
                }

EXP_SET_ID = 'all-experiments'
