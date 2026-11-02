# demo

An alternate copy of [`Initialization/`](../Initialization/README.md) and
[`Analysis/`](../Analysis/README.md). This copy matches the [original instructions](../DOCUMENTATION.md#original-instructions): real-cell inputs
in the RainbowSTORM format with z in column 12, `_C1`/`_C2` file names for the analysis, and the
mask values 0.17 / 500000. Function behaviour is documented in those two folders' READMEs. This
page lists only what differs. Known issues are in [../DOCUMENTATION.md](../DOCUMENTATION.md).

## Contents

| File | Compared with the main copy |
|---|---|
| `Initialization/calc_decay.m`, `compile_function_in.m`, `create_variables.m`, `epsilon_cost.m`, `exclude_bg.m`, `find_area.m`, `find_gt_width.m`, `find_within.m`, `find_within2.m`, `param_opt2_5.m`, `sim_histones_3d.m`, `sim_nucleus_3d.m` | identical |
| `Initialization/denstest_1.m` | input cell `OVCAR5_ALDH+_021` instead of `well_2_2` |
| `Initialization/simulation_optimization_2.m` | input cell `OVCAR5_ALDH-_k4_002`; **reads z from column 12** (main: column 5) |
| `Initialization/sim_densities_3.m` | **first removes `nmeth, nac, nseedm, nseeda, name, ratcool*` from `sim_param2.mat`**, so `sim_histones_3d_2` keeps the ×n density scaling (main copy lacks this) |
| `Initialization/sim_histones_3d_2.m` | **lacks** the main copy's fix that keeps the output name after `load('sim_param2')`. `sim_densities_3` here clears `name` from that file but not `filename`, so every density level is written to the saved name (`Simulation.csv`), each overwriting the previous one |
| `Initialization/split_me_ac.m`, `split_script_4.m` | **not present**; `simulation_optimization_2.m` calls `split_me_ac`, so copy it from the main folder |
| `Analysis/analyze_STORM.m` | reads `<cell>_C1.csv` / `_C2.csv`, z from column 12; mask threshold 0.17 (line 181) and minimum object 500 000 px (line 184); many extra commented-out plotting blocks; ring ratio stored as `ratio.n`. **Its final `save` (line 509) names `C1_within`, `C2_within`, … which are never created, so it errors**; rename them to the `meth_*` / `ac_*` names the main copy saves |
| `Analysis/compile_data.m` | single cell `k27_k27_thaw000.mat` → `k27_k27_thaw.mat`; expects the main copy's `meth_*` / `ac_*` variable names |
| `Analysis/compile_function.m` | identical |
| `Analysis/acetyl_epsilon_surround2.mat` | an example ac calibration: `coefa = 51.10`, `expa = −0.0766` (ε ≈ 51 nm · density^−0.077, density in localizations / µm²). The matching `methyl_epsilon_surround2.mat` is not included |

The six `.mat` inputs `Initialization/` needs (`initial_variables`, `mean_ratio`, `mean_mdens`,
`mean_adens`, `ratio`, `mdens2`), plus a `gt_values.mat`, were also deleted from this folder in
commit `ce89e0c` and can be restored from git (see the main documentation).
