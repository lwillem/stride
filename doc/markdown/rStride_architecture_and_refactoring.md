# rStride: As-Built Architecture and Refactoring Analysis

**Date:** 2026-09-25
**Branch analysed:** `measles_usa` (48 commits ahead of `master`, 1 behind)
**Scope:** the R workbench (`main/r/`) and its coupling to the C++ kernel build/install flow.
**Working checkout:** `~/Documents/university/research/stride/repo/stride_2023`. Findings
below were gathered from an identical checkout at the same commit; everything they rest on
— the git history, the shared install root `~/opt/stride-<N>`, and the machine's
toolchain — is common to both.

This document describes how the system **actually works today**, the reasoning
behind the problems observed, and a proposed refactoring sequence. It complements
`rStride_readme.Rmd`, which describes how rStride is *intended* to be used.

---

## 1. Current flow (as built)

### 1.1 From repository to a runnable workbench

```
repo/stride_2023/
  main/cpp/          C++ kernel            ──┐
  main/r/*.R         experiment scripts    ──┤
  main/r/rstride/    core R functions      ──┼── make install ──▶  ~/opt/stride-<N>/
  main/resources/    data, config          ──┘                        bin/stride      (binary)
                                                                      bin/*.R         (experiment scripts)
                                                                      bin/rstride/    (core R, copied)
                                                                      config/
                                                                      data/
                                                                      tests/          (regression .rds)
                                                                      sim_output/     (created at runtime)
```

The install prefix is derived in `Makefile:35-36`:

```make
LABEL = $(shell git rev-list HEAD --count)
CMAKE_INSTALL_PREFIX = $(HOME)/opt/stride-$(LABEL)
```

**The install root is a function of the commit count.** Every commit produces a new
install directory. Observed on this machine: `stride-745, 746, 748, 798, 806, 815, 820`.

Experiment scripts are installed by an explicit file list in `main/r/CMakeLists.txt`;
the `rstride/` library is installed wholesale as a directory (excluding `*.Rproj*`
and `*.Rhistory`).

### 1.2 Running an experiment

The user changes directory to the install root and executes a script from there,
because all internal paths are relative to that root:

```
cd ~/opt/stride-820
./bin/rStride_contacts.R
```

`.rstride$set_wd()` (`Misc.R:588`) automates this: it scans `$HOME/opt`, parses the
numeric suffix of `stride-<N>`, and `setwd()`s to the highest one. It falls back to
`$VSC_SCRATCH` on the UA cluster.

### 1.3 The R load sequence

Every experiment script begins with the same two lines:

```r
rm(list=ls())
source('./bin/rstride/rStride.R')
```

`rStride.R` then:

1. Installs/loads `simid.rtools` from GitHub (version >= 0.1.43), installing it if absent.
2. Force-loads 20 CRAN packages via `smd_load_packages()`, including `sf`, `tigris`,
   `usmap`, `haven`, `VGAM` (needed only by the USA population factory).
3. Sources `Misc.R`, which creates the `.rstride` environment.
4. Sources every `.R` file under `./bin/rstride` found by `dir(recursive=TRUE)`,
   minus a hardcoded exclusion list.
5. Sets `options(scipen=999)` to avoid scientific notation reaching the C++ layer.
6. Defines the public controller functions.
7. As its **final statement**, captures the whole global environment:
   `rStride_functions <- ls(all.names = TRUE)`.

### 1.4 The experiment pipeline

```
exp_param_list   (named list of parameter vectors, set in the experiment script)
      │
      ▼  .rstride$get_full_grid_exp_design()
exp_design       (data.frame: one row per run, full-factorial × rng seeds)
      │
      ▼  run_rStride()
      ├─ validation gate: data_files_exist, log_levels_exist, valid_r0_values,
      │                   valid_immunity_profiles, valid_seed_infected, valid_cnt_param
      ├─ create project dir: sim_output/<timestamp><dir_postfix>/
      ├─ smd_start_cluster()
      ├─ foreach (%dopar%) over rows of exp_design:
      │     ├─ create_config_exp()             merge defaults + design row
      │     ├─ integrate_parameters_in_calendar()
      │     ├─ save_config_xml()               write <exp_tag>.xml
      │     ├─ system('./bin/stride -c <xml>') run the C++ kernel
      │     ├─ read summary.csv, cbind with config
      │     └─ parse_log_file()                event_log.txt ──▶ <exp_tag>_parsed.rds
      ├─ write <run_tag>_summary.csv
      ├─ .rstride$aggregate_compressed_output()
      └─ smd_stop_cluster()
      │
      ▼  returns project_dir
inspect_*(project_dir)   8 post-processing entry points, uniform signature:
      inspect_summary, inspect_participant_data, inspect_contact_data,
      inspect_transmission_data, inspect_transmission_dynamics,
      inspect_incidence_data, inspect_prevalence_data, inspect_tracing_data
```

Parallel workers receive the library via
`foreach(..., .export = c('par_nodes_info', rStride_functions))` — i.e. by exporting
the captured global environment from step 7 above.

### 1.5 The regression test flow

The load-bearing test suite is **`rStride_gtester_covid19.R`** (717 lines), not the
C++ gtester (426 lines in `test/cpp/gtester/`).

```
22 scenario designs (covid_base, covid_hhcl, covid_tracing, covid_collectivity,
   covid_airborne, covid_subpools, covid_fitting, ...) × 5 rng seeds
      │
      ▼ run_rStride()
      ▼ compare against reference .rds in ./tests/
        summary / incidence / prevalence / contacts / participants / out_abc
      ├─ exact equality diff per column
      ├─ order-of-magnitude reporting for floating-point drift
      ├─ reports which gtester_label diverged
      └─ run-time comparison (performance regression)

rrv()       resets the reference .rds in ./tests/
rrv_repo()  additionally writes them back into the repository
```

Reference files live in the repo at `main/resources/rstride_test/`.

### 1.6 The USA population / contact-matrix flow

A separate, manually-run pipeline that is **not** part of the loaded library:

```
~/opt/FRED_population_usa/   (manual download, RTI synthetic population)
      │
      ▼  PopulationFactory_USA.R :: getFREDdata(state, county, com_target_size, rng_seed)
pop_usa   (population matrix)
      │
      ▼  social_contacts_usa2026.R
      ├─ contactdata::contact_matrix()  Prem et al. 2020, by location
      ├─ convert to *conditional* rates (conditional on school enrolment / employment)
      ├─ estimate cluster-size adjustment factors by uniroot()
      └─ write: <run_tag>.csv, contact_matrix_<run_tag>.xml, *_METADATA.txt, *.pdf
                into sim_output/<timestamp>_<run_tag>/
      │
      ▼  C++ AgeContactProfile.cpp:53
         reads matrices.adjustment_factor.<type>.value and scales the whole age profile
```

Note that `social_contacts_usa2026.R` is an *experiment script* that lives inside the
*library* directory `main/r/rstride/`, and is therefore explicitly excluded from the
library load.

---

## 2. Observed state

### 2.1 Code volume

| Layer | Files | Lines |
|---|---:|---:|
| C++ kernel (`main/cpp`) | 89 | 11,341 |
| R workbench (`main/r`) | ~37 | 13,668 |
| C++ tests (`test/cpp/gtester`) | 3 | 426 |
| R regression suite (`rStride_gtester_covid19.R`) | 1 | 717 |

The C++ figure includes **898 lines under `main/cpp/mdp/` that are not compiled** — no
`mdp/*.cpp` appears in `STRIDE_SRC`. Together with `main/python/` and `rStride_MDP.R`,
roughly 1,200 lines of the totals above are dormant; see F11.

### 2.2 Function inventory (R)

- 108 public functions defined into the global environment
- 40 private functions in the `.rstride` environment (40 defined, 40 distinct call sites — no dead weight)
- 6 total occurrences of `stop()` / `warning()` / `tryCatch()` across the whole R layer
- 16 `if(0==1){ attach(...) }` interactive-debug blocks

### 2.3 Longest functions

| Lines | Function |
|---:|---|
| 334 | `CalendarFactory_USA.R :: create_calendar_file` |
| 310 | `TransmissionAnalyst.R :: analyse_transmission_data_for_r0` |
| 308 | `CalendarFactory.R :: create_calendar_file` |
| 262 | `ParameterEstimator.R :: select_ensemble_and_plot` |
| 238 | `rStride.R :: run_rStride` (with a ~100-line `foreach` body inline) |
| 237 | `CalendarFactory_testing.R :: create_calenders_universal_testing` |
| 232 | `HealthEconomist.R :: calculate_cost_effectiveness` |
| 232 | `HealthEconomist_USA.R :: calculate_cost_effectiveness` |

### 2.4 Duplicated ("twin") files

| Original | Fork(s) | Lines |
|---|---|---|
| `HealthEconomist.R` | `HealthEconomist_USA.R` | 557 / 559 |
| `CalendarFactory.R` | `CalendarFactory_USA.R`, `CalendarFactory_testing.R` | 751 / 777 / 283 |
| `ImmunityProfileFactory.R` | `ImmunityProfileFactory_USA.R` | 127 / 86 |
| `rStride_default_param.R` | `rStride_covid19_default_param.R`, `rStride_measles_default_param.R` | 210 / 131 / 131 |

`calculate_cost_effectiveness` being 232 lines in *both* HealthEconomist variants
indicates the same function with differing constants.

---

## 3. Findings and reasoning

### F1. The install directory moves on every commit

`LABEL = git rev-list HEAD --count` means one commit changes the install root.
Because `sim_output/` is created *inside* the install root, every commit strands the
previous outputs in the old directory.

**This is the root cause of several downstream symptoms:**

- Experiment scripts accumulate commented-out absolute-ish paths such as
  `sim_output/20260827_234016_pop_usa_tx_gaines_c1000/...`, which are valid only in
  `stride-820` and break on the next commit.
- Generated artefacts get committed into `main/resources/data/`
  (e.g. `pop_usa_wisconsin_dane474k_c1000.csv`, `contact_matrix_usa_conditional.xml`)
  purely so they survive a commit. Large generated data therefore enters git.
- `.rstride$set_wd()` exists solely to cope with the moving target.

Secondary risk: `git rev-list --count` is **branch-dependent and not unique**. Two
branches, or a worktree, can produce the same count and silently share an install root.

### F2. rStride is a package that was never allowed to become one

148 functions, a hand-rolled private namespace (`.rstride`), a clean library/script
separation and a uniform public API — all the structure of a package, loaded by
`source()`-ing a `dir()` listing filtered through a hardcoded blacklist.

The blacklist has already rotted: `TransmissionInspector_old.R` is listed for exclusion
but **no longer exists**. It exists at all because an experiment script
(`social_contacts_usa2026.R`) lives in the library directory.

Consequences:

- **No namespace.** 108 names occupy the user's global environment.
- **No dependency declaration.** 20 packages load unconditionally, including heavy
  geospatial ones (`sf`, `tigris`, `usmap`) required only by a file that is *excluded
  from loading*.
- **No test framework.** `testthat` is unavailable, which is precisely what the
  regression work needs.
- **Load order is alphabetical** by `dir()`; any load-time dependency between files
  works only by luck.

### F3. The parallel layer is load-bearing on a clean global environment

`rStride_functions <- ls(all.names=TRUE)` captures globalenv at load time and feeds it
to `foreach(.export=)`. This is why every experiment script must open with
`rm(list=ls())`: otherwise user objects are serialised to every worker.

The convention works but couples three unrelated things — user workspace hygiene,
library loading, and parallel dispatch. As a package this line disappears entirely
(`.packages='rStride'`).

### F4. Kernel failures are silent

`rStride.R:385`:

```r
system(cmd, ignore.stdout = ignore_stdout)
run_summary <- read.table(summary_filename, header=T, sep=',')
```

The exit status is discarded. A kernel crash surfaces as an opaque `foreach` error
about a missing `summary.csv`, from inside a worker, rather than
"stride exited 139 on experiment 7". With 6 error-handling constructs in 13,668 lines,
this is the highest-friction defect in the workbench for long overnight grids.

### F5. Paths are hardcoded inside the library

`'./bin/rstride/Misc.R'`, `'./bin/stride'`, `'./config/run_default.xml'`, `'sim_output'`
are embedded in the library rather than injected. The library can only execute from one
directory. This is what makes F1 painful rather than merely untidy.

### F6. The regression harness is sound but mis-scoped

The design is correct: 22 scenarios, 6 output streams, exact diffs plus float-magnitude
reporting, reference reset, performance tracking. Three problems:

1. **Coverage is COVID/Belgium only.** No scenario exercises USA populations, measles,
   or the contact adjustment factors — i.e. none of the code currently being changed.
2. **References are internally inconsistent.** In `~/opt/stride-820/tests/`:
   `contacts`, `out_abc`, `participants` are dated Jul 6; `incidence`, `prevalence`,
   `summary` are dated Aug 25. A partial `rrv()` was run, so different streams encode
   different code states.
3. **`rrv_repo()` hardcodes an absolute personal path** —
   `~/Documents/university/research/stride/repo/stride_2023/main/resources/rstride_test`.
   The path is correct for the author's own checkout and the function is labelled "local
   function for LW" in the source, so this is a fragility rather than a defect: the
   reference-promotion step only works for one person on one machine, and silently writes
   nothing useful for anyone else. It should derive the repository root rather than
   assume it.

The comparison engine is also entangled with the 22 scenario definitions in one script,
so adding a measles suite by the obvious route means copy-pasting ~700 lines — creating
another twin.

### F7. The install file list has drifted

`main/r/CMakeLists.txt` enumerates experiment scripts by hand. Three have fallen off and
are absent from `~/opt/stride-820/bin/`:

- `rStride_measles_default_param.R`
- `rStride_r0_measles.R`
- `rStride_param.R`

Relevant because `rStride_gtester_covid19.R` sources `./bin/rStride_covid19_default_param.R`;
a measles gtester would need `rStride_measles_default_param.R`, which does not install.

`rStride_readme.Rmd` already documents the manual maintenance of this list, and the
"edit in the repo, not in the install dir, or your changes are overwritten" workaround —
confirming the friction is known.

### F8. Uncommitted model change in the working tree

`Infector.cpp:265-271` switches the contact probability from the **minimum** of the two
age-specific probabilities to their **average**, with the old block commented out:

```cpp
// double contact_probability = individual_contact_probability_p1 ;
// if(individual_contact_probability_p2 < individual_contact_probability_p1){
//     contact_probability = individual_contact_probability_p2;
// }
double contact_probability = (individual_contact_probability_p1 + individual_contact_probability_p2 ) / 2;
```

This is a substantive change to core transmission. It interacts directly with the
cluster-size adjustment factors being calibrated in `social_contacts_usa2026.R`
(averaging raises realised contacts relative to the minimum rule) and it will alter
every one of the 22 regression scenarios. It must be resolved before any reference
reset, because a golden master cannot be established against a moving target.

### F9. Minor items

- `social_contacts_usa2026.R:193` computes `workplace_threshold`, which is never used.
- Vendored third-party code under `main/r/rstride/lib/` (`npsp`, `cascsim`, `socrates`)
  carries its own `DESCRIPTION` files; licence compatibility with GPL-3 is worth checking.
- Commit `87fe899` records an unresolved issue: the school adjustment factor's
  "age-specific behaviour". Likely cause: `adj_fctr_school` is a single scalar applied to
  all ages by `AgeContactProfile`, while the target
  (`social_contacts_usa2026.R:140`) is an unweighted mean over enrolled ages and the
  school-size distribution is pooled across all ages — primary and secondary schools
  differ substantially in size.

### F10. The contact rule is a calibration dependency for every disease file

The change in F8 is not local to `Infector.cpp`. Each disease configuration carries a
fitted regression of R0 on the transmission probability:

```xml
<transmission>
  <b0>0.820990685751851</b0>
  <b1>18.91204185976</b1>
  <b2>0</b2>
</transmission>
```

`TransmissionProfile.cpp:57-77` inverts `E(R0) = b0 + b1*p + b2*p^2` to derive the
transmission probability from a requested R0. All ten disease files in
`main/resources/data/` carry such a fit:

| File | b0 | b1 | b2 |
|---|---|---|---|
| `disease_covid19_age.xml` | 0.14744 | 43.960 | 0 |
| `disease_covid19_age_15min.xml` | 0.21936 | 29.283 | 0 |
| `disease_covid19_child.xml` | 0.04610 | 38.610 | 0 |
| `disease_covid19_lognorm.xml` | 0.12449 | 39.646 | 0 |
| `disease_covid19_lognorm_child.xml` | 0.06896 | 34.571 | 0 |
| `disease_influenza.xml` | 0 | 37.481 | -27.444 |
| `disease_influenza_15touch.xml` | 0 | 15.990 | -7.738 |
| `disease_measles_adaptive_behavior.xml` | 0 | 40.948 | -14.680 |
| `disease_measles_adaptive_behavior_15min.xml` | 0 | 34.044 | -14.232 |
| `disease_measles_usa.xml` | 0.82099 | 18.912 | 0 |

These fits were produced against the **current `min` implementation**. Switching to
`mean` raises realised contacts per person, so the same transmission probability yields
a higher R0 and every fit becomes silently wrong: a run requesting `r0 = 12` would
receive a different effective R0, with no warning.

**The calibration provenance is already broken.** The `<label>` block of
`disease_measles_usa.xml` records the fitting run — `num_rng_seeds=45`,
`dim_exp_design=81` (approximately 3,645 runs), `fit_r0_range=0-24`,
`run_tag=20260727_171147_r0` — and references

```xml
<population_file>sim_output/20260710_111118_WI-Dane_conditional_social_contacts/...csv</population_file>
```

which lies in a superseded install directory and no longer exists. This is F1
materialising as a reproducibility failure in a committed artefact.

**Options:**

| | Action | Cost | Consequence |
|---|---|---|---|
| A | Adopt `mean`, refit all ten files | compute-bound (thousands of runs per file); scripted, but see the write-back caveat below | refits performed with the tooling described below |
| B | Revert to `min` | none | loses whatever motivated the change |
| C | **Make the rule configurable, default `min`** | small | existing calibrations remain valid; `mean` becomes opt-in |

**Option C is recommended.** It converts a globally-invalidating change into a local,
opt-in one: existing disease files and regression references stay valid, so the
refactoring retains its safety net, and the two rules become comparable on the same
population — which is impossible today, since `min` and `mean` are mutually exclusive
states of the working tree. Implementation is a configuration value read once and stored
as a member, branched on in `GetContactProbability`. The branch is loop-invariant; if
zero overhead is required, `InfectorExec.h` / `InfectorMap.h` already establish a
compile-time dispatch pattern.

**How refitting actually works.** The fit is produced by `bin/rStride_r0.R`, which runs
an R0 sweep and calls `analyse_transmission_data_for_r0()`
(`TransmissionAnalyst.R`). That function regresses secondary cases on the transmission
probability, writes the coefficients into the parsed disease configuration
(lines 280-282) together with a `<label>` provenance block, and saves the result as

```r
disease_config_update_file <- paste0(run_tag,'_',basename(disease_config_file))
.rstride$save_config_xml(config_disease,'disease',file.path(project_dir,disease_config_update_file))
```

Two consequences:

1. **The original file is never overwritten.** The refit lands as
   `sim_output/<timestamp>_r0/<run_tag>_disease_<name>.xml` and must be **copied by
   hand** into `main/resources/data/`. Since `sim_output` sits inside the install root,
   that output is orphaned by the next commit (F1). This is the mechanism that produced
   the stale provenance in `disease_measles_usa.xml`: it carries
   `run_tag=20260727_171147_r0` and a `population_file` path that no longer exists.
   Relocating `sim_output` (Phase 3 step 6) and scripting the copy-back are prerequisites
   for a defensible calibration.

2. **The quadratic term is currently disabled.** `TransmissionAnalyst.R:96` hardcodes
   `fit_b2 <- 0`, with `mod$coefficients[3,1]` commented out on the following line. Four
   committed files carry a non-zero `b2` (`disease_influenza.xml`,
   `disease_influenza_15touch.xml`, `disease_measles_adaptive_behavior.xml`,
   `disease_measles_adaptive_behavior_15min.xml`) and therefore originate from an earlier
   quadratic version of this function. Refitting them with the current code would
   silently change their functional form from quadratic to linear — and
   `TransmissionProfile.cpp:62` branches on `b2 == 0` to choose between the linear and
   quadratic inversion. Whether the quadratic term should be restored is a modelling
   decision that must be settled *before* any recalibration campaign, independently of
   the min/mean question.

**Ordering.** The *mechanism* must land before the §6.4 consolidation tag, because that
tag is the baseline against which the whole refactoring is verified, and regression
references cannot be reset against a moving target. The *recalibration campaign* should
not precede the refactoring, because refitting depends on
`rStride_r0.R` / `rStride_r0_measles.R` -> `TransmissionAnalyst::analyse_transmission_data_for_r0`
-> `ParameterEstimator`, which is currently among the least reliable code in the
workbench: `rStride_r0_measles.R` does not install at all (F7); the analysis function is
310 lines and untestable; `system()` swallows kernel failures (F4), so in a 3,645-run
fitting grid a subset can die silently and the fit is taken over the survivors; and fit
provenance points at vanished directories (F1). Refitting under those conditions produces
numbers that cannot be reproduced or defended — `disease_measles_usa.xml` already
demonstrates the failure mode.

### F11. A dormant Python / MDP layer, disabled in place

Roughly 1,200 lines implement a pybind11 interface exposing Stride as a step-wise,
controllable simulator for reinforcement-learning agents — an RL agent advances the
simulation a day at a time and applies vaccination actions between steps. It spans three
layers: `main/r/rStride_MDP.R` (configuration generator, 108 lines, still installed),
`main/python/pybind/` (bindings, ~155 lines plus build files), and `main/cpp/mdp/`
(step-wise control API, 898 lines).

**It is disabled in two independent places.** `main/python/CMakeLists.txt` contains
`#add_subdirectory(pybind)`, and no `mdp/*.cpp` appears in `STRIDE_SRC`, so the MDP core
is compiled only by the disabled pybind target. The code is present in the tree and
absent from every build artefact.

**Verification of its current state** (2026-09-25) found the opposite of what disabled
code usually looks like:

| Check | Result |
|---|---|
| `mdp/MDP.cpp`, `mdp/MDPRunner.cpp` compile against the current core | clean, 0 errors |
| CMake variables the pybind build needs, defined outside `pybind/` | all present |
| `libstride` target available for linking | yes |
| pybind wrappers consistent with `MDP.h` | yes |
| `test.py` consistent with the bound API | **no** — calls `Create()` with 1 of 3 required arguments |
| pybind11 dependency | pinned to `v2.6.0` (October 2020) |

So the C++ has not drifted from the core API despite two years outside the build. The
obstruction is environmental: the pybind11 pin predates Python 3.12 support, declares
CMake compatibility below 3.5 (rejected by CMake 4.0), uses the single-argument
`FetchContent_Populate` form removed in CMake 4.0, and requires network access at
configure time. On the development machine at the time of writing (CMake 3.24.4,
Python 3.8.3) it remained inside the supported window, which is why the build failed only
when a different interpreter was selected rather than unconditionally.

**Why this matters beyond the line count.** Disabling a subsystem in place, rather than
removing it with a specification, is the failure mode this plan exists to correct: the
code decays invisibly, nobody can tell whether it still works without reconstructing the
build, and the knowledge of what it was for lives only in the author's memory. Detail and
removal requirements are in `removed_features/python_mdp_interface.md`; the excision is
Phase 0b.

One residual coupling reaches the live core: `ContactPoolSys::ClearContactPools()` is
compiled into `libstride`, and its only caller is `MDP.cpp:454` (open decision 4).

### F12. The venue extension is intrusive, and the population format cannot absorb it

Four contact-pool types were added alongside the original set: `OtherHouse`, `RestoCafe`,
`OtherPlace` and `Transport`. They are **structurally different** from the classic pools,
and that difference is legitimate:

| | Classic pools | New venues |
|---|---|---|
| Pool id per person | one (`m_pool_ids[type][0]`) | one **per day of week** |
| Source | population CSV column | separate `subpools_community_file` |
| Contact rate | `AgeContactProfile`, by age | per-person `CPoolContacts(type)[day]` |
| Duration | not tracked | tracked per day |

The problem is not the extension. It is that the difference is expressed by **naming the
four types at every decision point** instead of declaring what makes them different.

#### F12.1 Identity used as a proxy for properties

The literal four-type list appears **nine times across five files**:

```
PopBuilder.cpp:155,165,188,240   SimBuilder.cpp:93   Sim.cpp:141   Infector.cpp:222
Person.cpp:106-109,126-134       ContactDivider.cpp:48-61,117-120 (unrolled by hand)
```

These lists encode **three different predicates that do not coincide**, which is invisible
because they look alike:

| Predicate | Sites | Members |
|---|---|---|
| loaded from subpool file / day-indexed | PopBuilder x4, `SimBuilder.cpp:93`, `Sim.cpp:141`, `Infector.cpp:222` | the 4 venues |
| has individual contact variation | `Infector.cpp:275` | 4 venues **+ Workplace + both Community**, **minus OtherHouse** |
| presence toggled per day | `Person.cpp:106-134` | the 4 venues |
| has physical venue characteristics (airborne) | `PoolCharacteristicsSeeder.cpp:94` | **School, Workplace, Collectivity + the 4 venues** |

`Infector.cpp:275` and `PoolCharacteristicsSeeder.cpp:94` are genuinely different sets and
nothing states so. Adding a fifth venue requires locating all the sites and knowing which
of the four groups it joins.

#### F12.2 The population file format is positional and count-inferred

`PopBuilder.cpp:82-88` determines the layout by probing a value and a column count:

```cpp
bool bool_profession = Trim(headers[2]) == "worker";           // layout from a VALUE
unsigned int profession_adj = bool_profession ? 2 : 0;
bool has_extra_column = headers.size() == (7+profession_adj);  // and from COLUMN COUNT
if (has_extra_column) extra_id = Trim(headers[6+profession_adj]);
bool household_cluster_id = extra_id == "household_cluster_id";
bool collectivity_id      = extra_id == "collectivity_id";
```

Five liabilities:

1. **Every field is read positionally** — `values[1+profession_adj]` and so on. Column
   order is load-bearing.
2. **The layout is inferred from the column count using `==`.** An eight-column file makes
   `has_extra_column` false and column 7 is **silently ignored**. Adding a column anywhere
   breaks detection with no error.
3. **`household_cluster_id` and `collectivity_id` are mutually exclusive.** Only one extra
   column is representable; both can never be present.
4. **Ragged rows degrade silently** — the per-row `if (values.size() == 7+profession_adj)`
   applies defaults instead of raising an error.
5. **No header validation** beyond those two probes.

The four new venues consequently could not be added to this file at all, which is why they
arrive through a separate `subpools_community_file`. The side-file is a reasonable
isolation, but it was forced by the format rather than chosen.

#### F12.3 Memory cost, measured

Compiled against the real headers on 2026-09-25:

```
NumOfTypes()   = 12   (the enum has 11 values)
sizeof(Person) = 1192 bytes
  IdSubscriptArray<array<unsigned int,7>> = 336 bytes each
  three of them (ids + durations + contacts) = 1008 bytes
```

**85% of every `Person` is those three arrays**, and most of it is never used:

| Waste | Per person |
|---|---:|
| `NumOfTypes()` returns 12 for 11 types — a phantom slot in all three arrays | 84 B |
| `m_pool_durations` + `m_pool_contacts` sized for 12 types, used by 4 | 448 B |
| `m_pool_ids` gives 7 days to classic types that use only day 0 | 168 B |
| **Total** | **~700 B of 1192 (59%)** |

At population scale: Dane WI (474k) holds roughly 565 MB of `Person` objects, about
**332 MB of it waste**; a 3M Belgian population roughly 3.6 GB with about 2.1 GB waste.
This is the probable reason `.rstride$print_system_memory_info()` exists in the R layer.

The `NumOfTypes()` discrepancy is a one-character fix worth about 40 MB at 474k, and is
almost certainly a remnant of the `College` type still visible commented out at
`ContactType.cpp:41`.

#### F12.4 Two defects found while investigating

1. **`ContactType::IsId()` can never return true** (`ContactType.cpp:52-53`): it
   upper-cases its argument and then looks it up in a mixed-case map. It currently has no
   callers, so the defect is dormant — but it is a trap for the next user.
2. **`ContactType::ToId()` ships debug output** — on an unrecognised string it dumps the
   entire map to `cerr` before throwing (`ContactType.cpp:80-86`). It also strips `"`
   characters from its input, a symptom of CSV quoting being handled ad hoc rather than by
   the parser.

#### F12.5 The venue characteristics are already generalized; the day array is redundant

The extension is less special than it appears. `ContactPool` **already** carries the venue
attributes, and `PoolCharacteristicsSeeder` **already** applies them to every type:

```cpp
double       GetVentilation() / SetVentilation(double)
double       GetAirMass()     / SetAirMass(double)
unsigned int GetDayWeek()     / SetDayWeek(unsigned int)
bool         IsNonComplier()  / SetNonComplier()
```

`pool_characteristics.xml` declares `ventilation_reduction` for Household, School,
Workplace, HouseholdCluster, Collectivity and the four venues on equal footing, and the
seeder loops `for (ContactType::Id typ : ContactType::IdList)`.

Duration is likewise **not** a venue-only field. `PoolCharacteristicsSeeder.cpp:148-161`
already populates it for classic types:

```cpp
if (typ == Id::School || typ == Id::Workplace) {
        for each member: for dayWeek 1..5: p->PoolDurations(typ)[dayWeek] = pool_duration;
} else if (typ == Id::Collectivity) {
        for each member: for dayWeek 0..6: p->PoolDurations(typ)[dayWeek] = pool_duration;
}
```

So `PoolDurations` has two loaders writing the same field: a type-level value from the XML
for School/Workplace/Collectivity, and a per-person-per-day value from the subpools file
for the venues.

**What remains genuinely venue-specific is therefore only two things:** per-day pool
membership, and a contact rate taken per-person rather than from the age profile. The
separate data file is a consequence of F12.2, not a modelling requirement.

**And the per-day membership array is itself redundant.** `Sim.cpp:141-146` already filters
on the *pool* side:

```cpp
const auto day_week_pool = poolSys.RefPools(typ)[i].GetDayWeek();
if (day_week_pool != dayWeek) { continue; }
```

Each venue pool carries its own day, so pool 17 *is* a Tuesday pool and `Sim` runs it only
on Tuesdays. A person's membership of pool 17 is therefore already a Tuesday fact:
`m_pool_ids[RestoCafe][2] == 17` duplicates `pool17.m_day_week` together with
`pool17.m_members`. That duplication is what costs the ~700 bytes per person in F12.3.

Per-individual non-attendance is currently encoded as **pool id 0** — `PopBuilder.cpp:257`
adds a member only when `subpool_id > 0`, while storing the zero in the person's array
regardless, and `Sim` iterates pools from index 1, so pool 0 is a null pool. Attendance
therefore already varies per individual per day, but through the id array rather than
through `m_in_pools`, which `UpdatePresence()` sets uniformly to `true` for every
non-isolated person.

**Three day-gating mechanisms consequently overlap:**

| Mechanism | Location | Applies to |
|---|---|---|
| type-level skip | `Sim.cpp:131-137` | School, Workplace, Community, HouseholdCluster |
| pool-level day match | `Sim.cpp:141-146` | the 4 venues |
| per-person presence | `Person::UpdatePresence()` | symptoms / isolation only |

The first two do the same job at different granularities. Phase 5b collapses them.

#### F12.6 Probable defect: one random draw shared by every pool

`PoolCharacteristicsSeeder::Seed()` samples once at function scope:

```cpp
auto uniform01Number = m_rn_man->at(0U).SampleUniform01();    // once, at the top
...
double variatie = uniform01Number * variability_area;          // inside the per-pool loop
double grootte  = average_area_per_person + variatie;
```

The single draw is reused for **every pool of every type**, so `variability_area` produces
no variation between pools — only a constant offset applied uniformly. If per-pool area
variability was intended, it is not occurring. This should be confirmed with the author
before changing, because correcting it will move results and therefore requires a
deliberate reference reset under the rule in §7.5.

### F13. Build configuration is dated, and two settings are actively harmful

The CMake files carry copyright dates of 2017-2019 and predate "modern CMake" throughout.
Findings below were verified against what the build **actually produces**, not only what
the files say. The flags reaching the compiler in a Release build are:

```
-g -fvisibility=hidden -std=c++17 -Wall -Wextra -pedantic -Weffc++ \
-Wno-unknown-pragmas -std=c++1z -O3 -DNDEBUG -ffast-math -arch arm64 ...
```

#### F13.1 `-ffast-math` undermines the regression strategy

`CMakeCPP.cmake:45` adds `-ffast-math` to every Release build. It permits floating-point
reassociation, flushes denormals and assumes no NaN or Inf.

The project's entire test strategy is **byte-exact comparison against stored `.rds`
references** (§1.5). Under `-ffast-math`, results may differ across compiler versions,
optimisation levels and architectures, so a colleague on a different machine can obtain
different numbers from identical source with no way to distinguish that from a genuine
regression. It also acts directly on the `std::exp()` calls in the airborne transmission
path.

Recommendation: remove it. If the speed is wanted back, `-fno-math-errno` alone is safe
for IEEE semantics and captures most of the benefit in `exp`/`log`-heavy code. Removal
**changes results** and therefore requires its own pull request with a deliberate
reference reset, and must precede the reset in Phase 4.

#### F13.2 OpenMP is not linked, and its absence is silent

Verified on the development machine:

| | |
|---|---|
| `~/opt/stride-820/bin/stride` | Mach-O **arm64**; `otool -L` shows only `libSystem`, `libc++` |
| `/usr/local/Cellar/libomp/17.0.6/lib/libomp.dylib` | Mach-O **x86_64** |
| `/opt/homebrew` (Apple Silicon prefix) | does not exist |
| `CMakeCache.txt` | `HAVE_FOUND_OpenMP:BOOL=FALSE`, `OpenMP_CXX_FLAGS:STRING=NOTFOUND` |

An arm64 binary cannot link an x86_64 library, so detection failed and the build fell back
to the dummy-OpenMP stubs in `main/resources/lib/domp`. Every `#pragma omp parallel`
compiles to nothing.

**The performance impact today is nil**, because `rStride.R:87` hardcodes
`config_default$num_threads <- 1` and no live script raises it — only `rStride_MDP.R:41`
sets 16, and that belongs to the dormant workflow of F11. rStride parallelises at the
*experiment* level through `foreach %dopar%`, spawning independent stride processes, which
is unaffected by OpenMP.

**The correctness impact is not nil, and is the reason to fix this.** Without OpenMP
locally the parallel code cannot be compiled, exercised or tested at all:

- **Half the C++ suite is a duplicate.** `main/resources/lib/domp` stubs return
  `omp_get_num_threads() = 1` and `omp_get_thread_num() = 0`, so
  `ConfigInfo::NumberAvailableThreads()` returns 1 and the gtester's two instantiations —
  `RunTest(tag, data, 1U)` and `RunTest(tag, data, ConfigInfo::NumberAvailableThreads())` —
  run identical single-threaded configurations. Of 44 test instances, 22 are redundant and
  the multi-threaded path has **zero coverage**.
- **A malformed `#pragma omp` cannot be detected.** `-Wno-unknown-pragmas`
  (`CMakeCPP.cmake:44`) suppresses the only diagnostic that would report a dropped or
  misspelled directive.
- **Data races and incorrect sharing clauses are invisible**, because the region never
  executes in parallel.
- **Thread-indexed code never leaves index 0.** `m_rn_man_ptr->at(thread_num)`
  (`Sim.cpp:108,127`) is only ever called with 0. If a configuration sets `num_threads`
  above the size the `RnMan` was constructed with, `.at()` throws — a failure that cannot
  occur locally but can on a machine where OpenMP is present.
- Anyone setting `num_threads > 1` gets no parallelism and no warning.

This matters because **OpenMP is live on a deployment target in actual use**: the VSC
cluster, for which `.rstride$is_ua_cluster()` exists, and `rStride_MDP.R` with
`num_threads = 16`. A defect in the parallel path would surface there, or in CI once the
workflows of §9 land, rather than on the developer's machine.

**A previous attempt to fix this is still in the tree and has never had any effect.**
`CMakeLocal.cmake` — the file included by `CMakeLists.txt:23` as the designated local
override hook — contains:

```cmake
set(LDFLAGS  -L/usr/local/opt/libomp/lib)
set(CPPFLAGS -I/usr/local/opt/libomp/include)
```

`LDFLAGS` and `CPPFLAGS` are autotools and make environment names; **CMake ignores both
entirely**. They also point at the Intel Homebrew prefix. Together with the cached
detection failure and the suppressed pragma warnings, this is why the gap went unnoticed.

**Root cause on the development machine:** the shell runs under **Rosetta 2** (`arch`
reports `i386`, `uname -m` reports `x86_64`) on Apple M2 hardware, so Homebrew and its
`libomp` are x86_64, while the build — configured from a native context — targets arm64.
A native shell is available (`arch -arm64 zsh`).

**Note a trap in the repair:** deleting `cmake-build-release/` and reconfiguring *from a
Rosetta shell* would make CMake detect `x86_64` and produce an Intel binary, which would
then link the existing x86_64 `libomp` successfully — OpenMP would appear fixed while the
whole simulator ran translated. Configuration must happen from a native arm64 shell.

Finally, enabling threads is **not** results-neutral: because the random-number manager is
indexed by `omp_get_thread_num()`, a run with `num_threads > 1` draws from different
streams than the same run with one thread. Repairing the *build* so the parallel path can
be compiled and tested is therefore a separate decision from *enabling* threads in
production, which would need its own reference set.

#### F13.3 Concrete defects

| Finding | Location |
|---|---|
| `-std=c++17` and `-std=c++1z` both passed; the `AppleClang` branch adds a second, redundant standard flag | `CMakeCPP.cmake:44, 70` |
| `CMAKE_ENABLE_COMPILE_COMMANDS` is **not a CMake variable** — the correct name is `CMAKE_EXPORT_COMPILE_COMMANDS`. The cache shows it empty, which is why no `compile_commands.json` exists for IDEs or clang-tidy | `CMakeConfig.cmake:25` |
| `-g` present in Release builds, from the cached `CMAKE_CXX_FLAGS` | build cache |
| Dead branch: `if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_HOST_APPLE)` never fires, because Xcode reports `AppleClang` | `CMakeCPP.cmake:60` |
| `USE_PYLIBSTRIDE` defaults `ON` but does nothing (the pybind subdirectory is commented out, F11), and its help text names a different variable, `USE_PYTHON` | `CMakeLists.txt:70` |
| `CMAKE_CXX_STANDARD` is never set — only `..._REQUIRED` and `..._EXTENSIONS`, with the standard supplied by a raw flag | `CMakeCPP.cmake:37-38` |
| Interprocedural optimisation disabled with the comment "issues with gcc 8.1.0 on travis"; Travis is no longer in use | `main/cpp/CMakeLists.txt:82-87` |

#### F13.4 Pre-modern CMake style

Directory-scoped `include_directories()` and `add_definitions()` rather than
`target_include_directories()` / `target_compile_definitions()`; `set(LIBS ${LIBS} x)`
accumulation rather than `target_link_libraries()`; hand-plumbed `${OpenMP_CXX_FLAGS}`
rather than the `OpenMP::OpenMP_CXX` imported target — which matters, because the current
code applies OpenMP *compile* flags but never its *link* flags.
`CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS` has been obsolete since CMake 2.x. The floor of
`cmake_minimum_required(VERSION 3.12...3.13)` is from 2018; raising it to 3.16 would make
`target_precompile_headers()` and `CMAKE_UNITY_BUILD` available, both of which would
reduce build time appreciably given the weight of the spdlog and Boost headers.

### F14. Repository history is dominated by re-committed binary data

`git count-objects -vH` reports a **270 MB** pack. The working tree's
`main/resources/data` is 213 MB, but 166 MB of that is
`pop_belgium3000k_c500_teachers_censushh_all.zip`, which is **gitignored**
(`.gitignore:123`) and therefore local-only — it has never been in a clone.

Aggregating every blob in history by path gives the actual composition.

**Deleted files that still cost every clone — about 57 MB:**

| Size | Path |
|---:|---|
| 35.8 MB | `main/resources/data/pop_belgium3000k_..._extended3_size2.zip` |
| 6.2 MB | `main/resources/data/pop_belgium600k_..._extended3.zip` |
| 5.7 MB | `main/resources/data/pop_flanders600.csv.zip` |
| 4.8 MB | `main/r/rstride/lib/simid_rtools-master.zip` |
| 4.5 MB | `main/resources/data/pop_belgium600k_c1k_teachers_censushh.zip` |

**Superseded versions of files still present at HEAD — about 50 MB:**

| Total | Versions | Path |
|---:|---:|---|
| 51.5 MB | **5** | `pop_belgium100k_c500_teachers_censushh.zip` (current version is 6.4 MB) |
| 9.9 MB | 2 | `pop_belgium600k_c500_teachers_censushh.zip` |
| 5.2 MB | 5 | `pop_belgium10k_c500_teachers_censushh.zip` |

And one that is small per version but diagnostic:
`main/resources/rstride_test/regression_rstride_incidence.rds` exists in **51 versions**
totalling 9.4 MB — every `rrv()` reference reset adds a new binary blob permanently.

So roughly **110 MB of the 270 MB is historical binary weight that nothing at HEAD
requires.**

**Root cause.** Population archives are *re-committed rather than versioned*: regenerate,
overwrite, commit, and the full size is added to history for good. This is the same
mechanism as F1 (outputs stranded in a moving install root) and F10 (calibration files
copied by hand): generated artefacts enter version control because there is nowhere else
durable to put them.

#### F14.1 A better home for population and other generated data

The project **already publishes synthetic populations on Zenodo** — the README cites
*Synthetic population data for Belgium*, `10.5281/zenodo.4485995`. The repository copies
are therefore largely duplicates of an existing, citable, versioned archive. Options, in
increasing order of effort:

1. **Fetch script plus checksums (recommended).** Keep a small manifest in the repository
   listing each population file with its DOI or URL, SHA-256 and target path, and a
   `make data` target that downloads and verifies. Version control then holds the
   *reference* to the data, not the data; a stale or corrupted download fails loudly; and
   the provenance requirement of §8.4 is satisfied by the same hash.
2. **GitHub Releases as asset storage** for files without a Zenodo record — keeps them
   near the code without entering the git object database.
3. **Git LFS** for files that genuinely must be tracked. Note that migrating *existing*
   files to LFS rewrites history and carries the cost in F14.2.

Keep only what tests need in-tree: the regression suite uses
`pop_belgium600k_...csv` (4.7 MB) and `pop_belgium10k_...csv` (1.6 MB). The 100k, 1000k
and extended variants are not exercised by any test and are the bulk of the weight.

#### F14.2 On rewriting history

`git filter-repo` would reclaim roughly 110 MB, but it rewrites every commit hash and so
simultaneously invalidates all six local branches, `origin`, and the three collaborator
remotes (`as/`, `ek/`, `ic/`). With active work on `measles_usa_rm` the coordination cost
greatly exceeds 110 MB of disk.

If it is done at all, the moment is immediately after the §6.4 consolidation tag, when
everyone is on a single branch, and it requires every collaborator to re-clone the same
day. **It is not a prerequisite for anything else in this plan**, and specifically not for
CI: `fetch-depth: 1` reduces a CI checkout to HEAD blobs only, roughly 47 MB once the
ignored archive is excluded.

---

## 4. Proposed refactoring sequence

Ordered so that each phase makes the next one safe. Each phase ends with a working system.

### Phase 0 — Settle the contact-probability rule (blocks everything)

0. Implement the contact rule as a configuration option defaulting to `min`
   (F8, F10, option C). Land it before the §6.4 consolidation tag. No disease file is
   refitted at this stage and no regression reference changes, because the default
   behaviour is unchanged.

### Phase 0b — Feature excision

Both excisions are behaviour-preserving, remove code that no build currently exercises,
and reduce the surface every later phase must carry. Each is specified in
`removed_features/` so it can be re-implemented later, and each must be a single
self-contained commit so `git revert` restores it. Sequence both **after** the §6.4
consolidation tag, so the tag remains a faithful snapshot of the last version containing
them.

| Feature | Specification | Footprint |
|---|---|---|
| `track_index_case` (TIC) | `removed_features/track_index_case.md` | template parameter + 5 branches; instantiations 10 -> 5 |
| Python / MDP interface | `removed_features/python_mdp_interface.md` | ~1,200 lines; already absent from every build artefact |

Notes:

- **TIC** was verified never to be enabled from the R workbench (`rStride.R:95` hardcodes
  `'false'`, no script overrides it). Removal requires a configuration tombstone, because
  Boost ptree silently ignores unknown keys and an old configuration setting
  `track_index_case = true` would otherwise run without the feature and produce quietly
  wrong results. Removal does not foreclose a future measurement use of the feature — the
  specification exists for that purpose — but it is explicitly **not** wanted for R0
  calibration, which must be fitted under a natural flow of infection
  (`removed_features/track_index_case.md` §2.1).
- **Python / MDP** (F11) is already disabled in two places: `main/python/CMakeLists.txt`
  comments out `add_subdirectory(pybind)`, and no `mdp/*.cpp` appears in `STRIDE_SRC`.
  The excision replaces disabled-in-place code, which decays silently, with a
  specification plus a clean deletion. Note that the C++ still **compiles cleanly** against
  the current core — this is not removal of broken code, but of working code that the
  configured build cannot produce, blocked by a pybind11 pin from 2020. Re-enabling is
  therefore not a matter of uncommenting the line; see
  `removed_features/python_mdp_interface.md` §3.1. Delete rather than comment out, since
  commenting out is precisely what produced this situation. One decision is required:
  `ContactPoolSys::ClearContactPools()` is compiled into `libstride` but its only caller
  is `MDP.cpp:454` (open decision 4).

### Phase 0c — Build configuration repair

Addresses F13.2 and F13.3. All steps here are behaviour-preserving: reference `.rds` files
must not change. Independent of the branch consolidation and safe to do at any time.

1. Correct `CMAKE_ENABLE_COMPILE_COMMANDS` to `CMAKE_EXPORT_COMPILE_COMMANDS`, producing
   `compile_commands.json` for IDEs and clang-tidy.
2. Set `CMAKE_CXX_STANDARD 17` and drop the raw `-std=` flags, removing the duplicate
   `-std=c++1z`.
3. Delete the unreachable `Clang AND CMAKE_HOST_APPLE` branch, or correct it to match
   `AppleClang`.
4. Remove the `USE_PYLIBSTRIDE` option, which does nothing (superseded by Phase 0b).
5. Drop `-g` from Release builds.
6. **Restore a working local OpenMP build.** Treat this as a correctness prerequisite,
   not a performance option: without it the parallel regions cannot be compiled or
   exercised, and 22 of the 44 C++ test instances are duplicates of the other 22 (F13.2).
   On the current development machine:

   1. use a native shell — `arch -arm64 zsh`, or disable *Open using Rosetta* on the
      terminal application;
   2. install native Homebrew at `/opt/homebrew` (it coexists with the Intel prefix) and
      `brew install libomp`;
   3. replace the inert `LDFLAGS` / `CPPFLAGS` lines in `CMakeLocal.cmake` with
      `set(OpenMP_ROOT /opt/homebrew/opt/libomp)`;
   4. **delete `cmake-build-release/` outright** — `HAVE_CHECKED_OpenMP` is cached, and
      `make clean` does not remove it;
   5. reconfigure **from the native shell**, or CMake will target `x86_64` and produce a
      translated build that links the Intel `libomp` and appears to work;
   6. expect a link failure and fix it: `CMakeCPP.cmake:112` adds `${OpenMP_CXX_FLAGS}` to
      the compile flags but never links the runtime. The correct remedy is the imported
      target, `target_link_libraries(libstride PUBLIC OpenMP::OpenMP_CXX)` (Phase 7b);
   7. verify with `otool -L ~/opt/stride-*/bin/stride | grep omp`.

7. **Make a failed OpenMP detection loud.** Report it prominently at configure time
   instead of silently substituting the `domp` stubs, and drop or narrow
   `-Wno-unknown-pragmas` so dropped pragmas are visible. Document that
   `HAVE_CHECKED_OpenMP` is cached and that recovery requires deleting the build
   directory.

Note that *enabling* threads in production remains a separate decision from repairing the
build: within-run threading changes random-number streams and therefore results (F13.2),
so it needs its own reference set. Once the build works, the gtester's multi-threaded
instantiation becomes meaningful for the first time and may need expected values of its
own, depending on whether the margins in `ScenarioData.cpp` already absorb the difference.

### Phase 1 — Package skeleton (non-destructive)

1. Add `DESCRIPTION` and `NAMESPACE` around the existing files; move nothing.
   Verify `devtools::load_all()` loads all 148 functions and exposes any hidden
   load-order dependency.
2. Move `social_contacts_usa2026.R` out of the library into `main/r/`. Delete the
   exclusion blacklist in `rStride.R`.
3. Declare dependencies in `Imports:`, moving `sf`/`tigris`/`usmap`/`haven`/`VGAM` to
   `Suggests:` so only the USA population factory pays for them.

### Phase 2 — Make failures visible

4. Wrap `system()` with exit-status checking and a clear per-experiment error naming
   the experiment id, config file and exit code. Apply to `rStride.R:385` and
   `rStride_main_abc.R:141`.
5. **Parse the population file by header name** (F12.2). Build a name-to-index map from
   the header row once, then read every field by name:

   ```cpp
   std::map<std::string,size_t> col;
   for (size_t i = 0; i < headers.size(); ++i) col[Trim(headers[i],"\"")] = i;
   const auto age = IntFromString(values[col.at("age")]);
   if (col.count("household_cluster_id")) { ... }
   if (col.count("collectivity_id"))      { ... }   // both may now coexist
   ```

   Column order becomes irrelevant, new columns are safe to add anywhere, the two extra
   pool ids stop being mutually exclusive, and a missing required column throws an error
   naming it rather than silently shifting every subsequent field. This belongs in this
   phase because the current behaviour is a *silent* misread, which is the same class of
   defect as a swallowed exit status. It is also small, and it protects every later phase
   that touches population input.

### Phase 2b — Calibration management and migration to the `mean` rule

`mean` is intended to become the default contact-probability rule. This phase delivers
that migration. It is gated on earlier phases rather than optional, and its design is
specified in §8.

Prerequisites:

- Phase 0 — the rule is configurable, default `min`.
- Phase 2 — `system()` error handling, so a failed run in a multi-thousand-run fitting
  grid cannot be silently excluded from the fit (F4).
- F7 install-list repair, so `rStride_r0_measles.R` reaches `bin/`.
- Phase 3 step 6 — a durable output location, so refits are not orphaned in a
  superseded install root (F1).
- A decision on the quadratic term disabled at `TransmissionAnalyst.R:96` (F10,
  open decision 3).

Steps:

1. Build the calibration infrastructure of §8: separated calibration artefacts,
   hash-based provenance, `promote_calibration()`, and the validation-gate check.
2. Refit the disease files under `mean` into new calibration artefacts. Keys differ by
   rule, so the existing `min` calibrations are never overwritten.
3. Compare the two rules on identical populations — possible only once both
   calibrations coexist.
4. Flip the default in a single dedicated PR that changes only the default and the
   regression references, carrying the step-3 comparison as justification (§7.5).
5. Retain the `min` calibrations in version control as history.

### Phase 3 — Decouple from the install directory

5. Replace hardcoded paths with a single `stride_paths()` object resolved once at
   `run_rStride()` entry.
6. Relocate `sim_output` outside the install root (`STRIDE_OUTPUT_DIR`, default
   `~/stride_runs/`). Add a `~/opt/stride-current` symlink.
7. Flip experiment scripts to `library(rStride)` and drop the R library copy from
   `main/r/CMakeLists.txt`. This also removes the "edits in `bin/rstride/` are destroyed
   on rebuild" trap.

### Phase 3b — Remove `-ffast-math`

Addresses F13.1. **Not behaviour-preserving** — it changes floating-point results — so it
requires its own pull request and a deliberate reference reset under the rule in §7.5.

Sequenced immediately before Phase 4 so that the single clean reference reset performed
there is taken against a build with deterministic IEEE semantics. Resetting references
first and removing the flag afterwards would invalidate them again.

### Phase 4 — Repair and extend the regression harness

8. Replace `rrv_repo()`'s hardcoded absolute path with a derived repository root, so
   reference promotion works for any checkout and any user (F6).
9. Resolve F8 (min vs average), commit it, then run one clean full `rrv()` so all six
   reference files derive from a single known commit. Record that hash alongside them.
10. Extract the comparison engine into `rstride/RegressionTester.R`, leaving
    `rStride_gtester_covid19.R` as scenario definitions only.
11. Add `rStride_gtester_measles.R` on that engine: USA populations, measles config,
    household clustering, contact adjustment factors. Add it and the three orphans
    from F7 to `main/r/CMakeLists.txt`.

Step 10 must precede step 11, otherwise the measles suite is created by copy-paste and
becomes another twin.

### Phase 5 — Collapse the twins

12. `ImmunityProfileFactory` (smallest, proves the pattern) → `HealthEconomist` →
    `CalendarFactory` (largest, ~1,800 lines across three files). Extract the shared
    body; push differences into a region/disease configuration object; delete the fork.
    One pair per commit, each verified green against the now-meaningful harness.

### Phase 5b — Unify the venue extension with the ordinary contact pools

**Goal: an efficient implementation of the extension, not a smaller feature set.** The
memory footprint of F12.3 arises from the *representation* rather than from the feature
itself, and is recovered in full by step 3 without changing anything the model can
express.

Two properties to keep in view while working here. The extension is exercised by two of
the twenty-two regression scenarios (`covid_subpools`, `covid_airborne`), so it is under
test and changes to it are caught. And airborne transmission is not independently usable:
all four of its gate sites test `m_airborne_transmission && m_subpools_community`, so the
two behave as one feature.

Addresses F12.1, F12.3, F12.4, F12.5 and F12.6. **The extension stays; it stops being a
category.** Per F12.5 the venue attributes are already generalized — `ContactPool` owns
ventilation, air mass, day and non-compliance, and `PoolDurations` is already populated
for School, Workplace and Collectivity. The remaining work is to remove the duplication
that makes the venues look special.

The target representation:

- **the pool owns the day.** `ContactPool::m_day_week` already exists and `Sim.cpp:141-146`
  already filters on it. Classic types become the degenerate case where every pool of the
  type shares a day (or is day-agnostic).
- **the pool owns attendance.** `m_members` already lists who attends; a parallel
  `std::vector<unsigned int>` holds each member's duration and computed contact count,
  sized to actual attendance rather than to 12 types x 7 days for every person alive.
- **the person owns only the gate.** `m_in_pools[type]` keeps exactly the role it has for
  School and Workplace: "I would attend today, but I am symptomatic or isolated."

Which day a person attends then *falls out of* membership plus the pool's day, rather than
being stored again per person. This is also strictly more expressive than the present
format: a membership list can represent attending two pools of one venue type on the same
day, which a single id per type per day cannot.

Sequence after Phase 4 — steps 2 and 3 alter `Person` and the transmission loops, and need
a trustworthy regression suite.

1. **Introduce a trait table** and replace identity tests with property tests:

   ```cpp
   struct PoolTraits {
           bool day_scoped;           // pools of this type carry a specific day
           bool uses_age_profile;     // rate from AgeContactProfile vs per-person
           bool individual_variation; // apply GetIndividualContactFactor
           bool has_venue_physics;    // air mass / duration / airborne transmission
   };
   constexpr IdSubscriptArray<PoolTraits> Traits = { /* one row per type */ };
   ```

   `Infector.cpp:222` becomes `Traits[pType].uses_age_profile`, `Infector.cpp:275` becomes
   `Traits[pType].individual_variation`, and `PoolCharacteristicsSeeder.cpp:94` becomes
   `Traits[typ].has_venue_physics`. This makes the divergence between those sets explicit
   rather than accidental, and reduces adding a venue to one row. Purely
   behaviour-preserving: reference `.rds` files must not change.

2. **Collapse the two day-gating mechanisms** (F12.5). Give every pool a day scope and let
   `Sim` apply one rule, replacing the type-level skip at `Sim.cpp:131-137` and the
   venue-only pool match at `Sim.cpp:141-146`. Behaviour-preserving.

3. **Move attendance data pool-side — the efficiency target.** Duration and per-person
   contact counts move from `Person::m_pool_durations` / `m_pool_contacts` into vectors
   parallel to `ContactPool::m_members`, and `m_pool_ids` reduces to one id per
   non-day-scoped type. `ContactDivider.cpp:48-61,117-120` inverts from iterating persons
   to iterating today's pools, which also matches how `Sim` is already parallelised.

   Cost moves from *per person x per type x per day* to *per actual attendance*:

   | | Now | Target |
   |---|---:|---:|
   | `m_pool_ids` | 336 B (12 types x 7 days) | 44 B (one id per type, 11 types) |
   | `m_pool_durations` | 336 B | 0 — pool-side |
   | `m_pool_contacts` | 336 B | 0 — pool-side |
   | other members | 184 B | 184 B |
   | **`sizeof(Person)`** | **1192 B** | **~230 B** |

   At 474k persons: roughly **565 MB -> ~110 MB**. The pool-side vectors cost one entry per
   attendance record rather than per person-type-day; at roughly four venue-days per person
   per week that is on the order of 15 MB, against ~455 MB recovered.

   The secondary benefit is cache behaviour. The transmission loop walks pool members
   through `Person*`, so shrinking `Person` by roughly a factor of five directly reduces
   cache misses in the hottest loop in the simulator. Measure it with the suite's existing
   per-scenario `run_time` tracking; treat any speedup as a welcome side effect rather than
   the justification.

   Strictly behaviour-preserving: reference `.rds` files must not change.

4. **Correct `NumOfTypes()` to 11** (F12.3). A one-character change worth about 40 MB at
   474k, independent of the rest and safe to land first.

5. **Fix the defects of F12.4** — `IsId()`, and the debug output plus ad-hoc quote
   stripping in `ToId()`.

6. **Resolve F12.6** — the shared random draw in `PoolCharacteristicsSeeder::Seed()`.
   Unlike the rest of this phase this is **not** behaviour-preserving: correcting it
   introduces genuine per-pool variability and will move results. It therefore belongs in
   its own PR with a deliberate reference reset, per §7.5, and requires confirmation from
   the author that per-pool variation was the intent.

**One question for the extension's author before step 3:** can a person attend two
different pools of the same venue type on a single day? The current format cannot express
it, so if the generator assumes that constraint the membership-list design is a superset
and nothing breaks — but it should be confirmed rather than assumed.

### Phase 6 — Configuration as data

13. One config file per study (e.g. `config/studies/measles_gaines_tx.yml`) holding
    population file, contact matrix, immunity, dates and clustering parameters.
    Experiment scripts reduce to a loader plus `inspect_*` calls. This makes a run's
    provenance recoverable, which it currently is not.

### Phase 7 — Internal structure

14. Extract the `foreach` body of `run_rStride()` into a testable
    `run_one_experiment()`.
15. Split `Infector.cpp` (668 lines). Make the contact-probability rule a named,
    switchable strategy rather than a commented-out block (already delivered as a
    configuration option in Phase 0; this step gives it a proper home).
16. Modernise the Infector template dispatch — see Phase 7a.

#### Phase 7a — Infector template dispatch

**Assessment: the pattern is sound; keep its core.** The transmission loop is dispatched
as follows:

| Stage | Location | Frequency |
|---|---|---|
| `InfectorMap` lookup resolved to an `InfectorExec*` | `SimBuilder.cpp:72` | once, at construction |
| Selection between default and tracing infector | `Sim.cpp:78` | once per timestep |
| Indirect call into `Infector<LL,TIC,TO>::Exec` | `Sim::TimeStep()` | once per contact pool |
| Contact evaluation loop | inside `Exec` | per candidate contact pair |

The indirection is therefore amortised over an entire contact pool and is not a
per-contact cost. This layering is correct and should not be changed.

`LOG_POLICY<LL>` (`Infector.cpp:101-201`) is the part that earns the design. For
`EventLogMode::None` its `Contact()` is an empty inline function, so in the innermost
loop it compiles to nothing at all — including the argument setup. A runtime
`if (logMode == ...)` cannot reliably match this, because the compiler would additionally
have to prove the arguments free of side effects before eliding them. This is
policy-based design used correctly and should be retained.

**What no longer earns its keep:**

1. **`TIC` doubles the instantiation count for negligible benefit.** It appears as a
   plain `if (TIC)` at `Infector.cpp` lines 406, 418, 485, 573 and 638 — all within
   transmission handling, which executes only when a transmission actually occurs, not in
   the inner contact loop. Because it is already written as an ordinary `if` rather than
   requiring specialisation, moving it from template parameter to function argument is
   close to mechanical and reduces explicit instantiations from **ten to five**.

2. **The two `Exec` bodies are approximately 73% duplicated.** The general form
   (`Infector.cpp:330-498`) is 169 lines; the optimised form (`504-651`) is 148 lines.
   Whitespace-normalised, only 84 of 317 lines differ. The core transmission algorithm
   therefore exists twice and is kept in step by hand. This is the C++ counterpart of the
   twin-file problem in §2.4, and it carries the same correctness hazard: an edit applied
   to one body and not the other. The `min` -> `mean` change of F8 was unaffected only
   because `GetContactProbability` is a shared free function — placement, not structure.

3. **The code predates the language feature that resolves this.** `CMakeCPP.cmake:44`
   sets `-std=c++17`, but these headers date from 2018-2020 and use the pre-C++17 idiom.
   With `if constexpr` the two bodies collapse into one function:

   ```cpp
   if constexpr (TO) { /* optimised path */ } else { /* general path */ }
   ```

   The discarded branch is not instantiated, so every bit of the compile-time elimination
   is preserved while the duplication disappears.

**On the general question.** Template metaprogramming as a style has moved on: SFINAE,
tag dispatch and explicit specialisation are largely superseded by `if constexpr` and
concepts. Policy-based design for zero-cost customisation, however, is not dated — it
remains the correct way to express this. What is dated here is the expression, not the
idea. This is not a matter of removing cleverness, but of saying the same thing in the
language standard the project already compiles against.

**Steps, in order:**

1. `TIC` is deleted in Phase 0b rather than converted to a runtime argument; see
   `removed_features/track_index_case.md`. Instantiations drop from ten to five and the
   `InfectorMap` key simplifies from a tuple to a plain enum.
2. Merge the two `Exec` bodies behind `if constexpr (TO)`. This is the step that matters,
   because it removes a standing correctness hazard rather than only reducing code size.
3. Leave `LOG_POLICY` as it stands.
4. Do **not** template the contact-probability rule. It is a runtime configuration value,
   so a member variable is correct; templating it would multiply instantiations again and
   `if constexpr` does not apply to a value unknown at compile time.

**Verification.** The R regression suite already records `run_time` per scenario and
reports the largest difference between runs, so step 1 needs a suite run rather than new
benchmarking tooling. Per §7.5, none of these changes may alter the reference `.rds`
files — which is precisely the correct test for a behaviour-preserving template refactor.

#### Phase 7b — CMake modernisation

Addresses F13.4. Migrate from directory-scoped to target-scoped commands
(`target_include_directories`, `target_compile_definitions`, `target_link_libraries` in
place of the `LIBS` variable), adopt the `OpenMP::OpenMP_CXX` imported target so link
flags are handled as well as compile flags, drop `CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS`, and
raise `cmake_minimum_required` to 3.16 to enable `target_precompile_headers()` and
`CMAKE_UNITY_BUILD`. Revisit the interprocedural-optimisation block disabled for a
Travis-era compiler (F13.3).

### Working rules

- **Never mix a behaviour change and a structure change in one commit.** In a simulation
  codebase this is the difference between a bisectable history and an unexplainable result.
- **Every phase ends green.** A stalled phase must leave a working simulator.

---

## 5. Open decisions

1. **`Infector.cpp` min vs average (F8, F10)** — **decided.** `mean` is to become the
   default. It is reached by migration, not cutover: the rule becomes configurable with
   `min` as the default in Phase 0, and the default is flipped in Phase 2b once `mean`
   calibrations exist and have been compared against `min`. No open question remains on
   the destination; the remaining choices are timing and the quadratic term (decision 3).
2. **Branch reconciliation** — resolved into a concrete plan; see §6. The remaining
   judgement calls are which of the stale local branches and collaborator remotes
   (`as/`, `ek/`, `ic/`) are live, and who owns each conflict resolution in §6.3.
3. **Quadratic term in the R0 fit (F10)** — `TransmissionAnalyst.R:96` hardcodes
   `fit_b2 <- 0`. Four committed disease files carry a non-zero `b2` from an earlier
   quadratic fit. Restore the quadratic term, or accept linear-only and refit those four
   deliberately? Must be settled before Phase 2b.
4. **`ClearContactPools()` after the MDP excision** — compiled into `libstride`, sole
   caller `MDP.cpp:454`. Delete alongside the MDP removal, or retain deliberately with a
   recorded reason? Must not be left undecided.
5. **Package installation route** — install `rStride` into the R library at build time,
   or `devtools::load_all()` against the repo during development? Affects Phase 1 step 7
   and how collaborators set up.
6. **Generated data in git** — once `sim_output` is relocated (Phase 3), should the
   generated population CSVs currently in `main/resources/data/` be removed from version
   control in favour of a documented regeneration step?
_(A question about whether R0 fitting should use index-case tracking was raised and
resolved: it should not. The fit is a calibration mapping and must be derived under the
same natural flow of infection under which it is applied, including competition between
secondary cases for remaining susceptibles. See
`removed_features/track_index_case.md` §2.1.)_

---

## 6. Branch reconciliation (prerequisite to the refactoring)

The refactoring in §4 must start from a single consolidated baseline. Two branches
carry live measles/USA work and have diverged for two months.

### 6.1 State as of 2026-09-25

| Branch | Tip | Date | Author |
|---|---|---|---|
| `master` / `origin/master` | `c79a76f` | 2026-08-19 | LieseB-1746743 |
| `measles_usa` / `origin/measles_usa` | `6226413` | 2026-09-09 | lwillem |
| `measles_usa_rm` (**local, stale**) | `151f11e` | 2026-08-25 | Regina Manansala |
| `origin/measles_usa_rm` | `a8f5095` | 2026-09-23 | Regina Manansala |

The local copy of `measles_usa_rm` is **21 commits behind** its remote. All analysis
below uses `origin/measles_usa_rm`. Run `git fetch --all` before acting on any of this.

Divergence, measured from the merge base `c2e209f` (2026-07-23, "Update files to model
Dane Co"):

| Comparison | Left-only | Right-only |
|---|---:|---:|
| `master` ↔ `measles_usa` | 1 | 48 |
| `master` ↔ `origin/measles_usa_rm` | 1 | 51 |
| `measles_usa` ↔ `origin/measles_usa_rm` | 32 | 35 |

### 6.2 Critical: the School presence fix is missing from one branch

`master`'s single exclusive commit `c79a76f` adds one line to
`Person::UpdatePresence()`:

```cpp
m_in_pools[Id::School] = true;
```

Without it, a person who stays home from school — for example after developing
symptoms — is never marked present at school again for the remainder of the
simulation.

| Branch | Fix present |
|---|---|
| `master` | yes |
| `measles_usa` | yes (same one-line change, applied independently) |
| `origin/measles_usa_rm` | **no** |

For a school-driven measles model this affects results, not just behaviour. Any
simulation run from `measles_usa_rm` since the merge base should be treated as
suspect until re-run. **This should be communicated before, and independently of,
the merge work.**

### 6.3 Conflict surface

A dry-run merge (`git merge-tree --write-tree measles_usa origin/measles_usa_rm`,
object-database only) reports **five conflicting files, all R**:

| File | Type | Notes |
|---|---|---|
| `main/r/rstride/social_contacts_usa2026.R` | content | The substantive one. They generalised `get_cnt_data(max_age=)` and moved to CT/Fairfield; this branch reworked it for TX/Gaines with `uniroot` adjustment factors. Both independently wrote small-workplace removal — live on theirs, commented out here at lines 171-174. Needs a modelling decision, not a textual merge. |
| `main/r/rStride_measles_explore.R` | content | Both heavily edited (~136 / ~141 lines). |
| `main/r/rstride/rStride.R` | content | Same block: both changed `recursive=F→T` and edited the exclusion blacklist. Mechanical. |
| `main/r/rStride_r0_measles.R` | content | Small (2 vs 32 lines). |
| `main/r/rstride/USA_PopulationBuilder.R` | **modify/delete** | Structural. This branch renamed it to `factories/PopulationFactory_USA.R` (`6ccf4d6`); theirs edited it in place (+20). Git cannot resolve; their changes must be ported onto the renamed file by hand. |

**There are no C++ conflicts.** Their `ImmunitySeeder` rework (+167 lines, including the
hang fix when `immunity_link_probability` is set) and new `PopSnapshotWriter` (+256
lines) merge cleanly with this branch's `AgeContactProfile` adjustment-factor work (+18).

Five files exist on both branches **byte-identical** — `HealthEconomist_USA.R`,
`factories/ImmunityProfileFactory_USA.R`, `immunity_measles_dummy.xml`,
`disease_measles_usa.xml`, `Misc.R`. These were synchronised outside git (copied rather
than merged). The practice has kept the conflict count low so far, but it is how work
gets silently lost; see §7.6.

### 6.4 Consolidation sequence

```
1. git fetch --all
2. merge master → measles_usa_rm          # the School fix (§6.2); do first, independently
3. merge master → measles_usa             # near no-op, same line already present
4. branch integration/measles-usa from measles_usa
   merge measles_usa_rm into it           # resolve the 5 conflicts of §6.3
5. build; run the C++ gtester; run the R regression suite
6. tag pre-refactor-2026-09; merge to master
7. all refactoring branches start from master
```

Notes on step 4: resolve it **jointly with the author of the other branch** — the
`social_contacts_usa2026.R` conflict requires deciding whose workplace treatment is
correct, which is a modelling judgement. Resolve F8 (`Infector.cpp` min vs average) in
the same session, so it does not enter the baseline unreviewed.

Note on step 5: per F6 the reference `.rds` files are internally inconsistent, so they
cannot cleanly validate this merge. Expect diffs, attribute each one deliberately, then
perform a single reference reset on the tagged baseline. That reset is the point at
which the golden master becomes trustworthy again.

---

## 7. Proposed branching and pull-request model

The current model is a direct cause of the situation in §6: long-lived personal branches
(`measles_usa`, `measles_usa_rm`, `superspreading`, `dev`, plus `as/`, `ek/`, `ic/`
remotes carrying parallel `measles`/`universal`/`GeoClustering` forks), two months of
unmerged parallel work, and file-level synchronisation outside version control.

### 7.1 `master` is the trunk

Everything lands on `master`. It must always build and pass the C++ gtester. Stale
branches and collaborator remotes should be explicitly archived (by tag) or deleted,
not left ambiguous.

### 7.2 Separate study branches from code branches

This is the most important change. `measles_usa` is currently both a kernel-development
branch and the Gaines TX study, which is precisely why it is two months old and cannot
merge.

| Kind | Naming | Lifetime | Merges to master |
|---|---|---|---|
| Code | `feature/…`, `fix/…`, `refactor/…` | days to ~2 weeks | yes, always |
| Study | `study/measles-gaines-tx` | indefinite | no — rebases onto master |

A study branch holds run configurations and results. A code branch holds changes to the
kernel or the workbench. Conflating them is what makes both unmergeable.

### 7.3 Pull requests for everything

Including self-merged ones. With a team of this size the value is not gatekeeping but
the written rationale. The repository's commit subjects are good, but the reasoning
behind model changes — min vs average contact probability, the household clustering
ratio, the cluster-size adjustment factors — is not recorded anywhere. A PR body is the
cheapest durable place for it.

### 7.4 Continuous integration

`.travis.yml` is present but presumed dead. Move to GitHub Actions:

- **Per PR:** build + C++ gtester.
- **Nightly on master:** the full R regression suite (too slow for per-PR).

### 7.5 The rule that makes the refactoring safe

> **A refactoring PR must not modify the regression reference `.rds` files.**

This turns "behaviour-preserving" from a promise into something CI can check, using the
harness that already exists. If a change genuinely requires new references, it belongs in
a separate PR whose sole purpose is that change, with the justification in the body.
This is the §4 working rule — never mix behaviour and structure in one commit —
expressed in enforceable form.

### 7.6 No out-of-band file synchronisation

The five byte-identical files in §6.3 should each have been a pull request. Copying files
between checkouts bypasses history, review and attribution, and eventually overwrites
someone's work without a trace.

---

## 8. Calibration management

Specification for the infrastructure delivered in Phase 2b. It replaces the current
practice of writing a refitted disease file into `sim_output/` and copying it into
`main/resources/data/` by hand.

### 8.1 What is wrong with the current arrangement

The existing practice has one correct property that must be preserved: **a refit never
overwrites a committed disease file**. The refit is written as
`sim_output/<timestamp>_r0/<run_tag>_disease_<name>.xml` and promoted manually.

Three properties are missing:

1. **The output is volatile.** `sim_output` lives inside the install root, which moves on
   every commit (F1). An unpromoted refit is orphaned.
2. **Promotion is unrecorded.** The copy is a manual step with no trace of who promoted
   what, when, or on what evidence.
3. **Provenance names paths, not content.** `disease_measles_usa.xml` records
   `population_file = sim_output/20260710_111118_WI-Dane_.../...csv`, which no longer
   resolves. The fit cannot be reproduced or even verified.

### 8.2 Root cause: two different kinds of thing in one file

`disease_*.xml` conflates:

| Content | Nature | Should be |
|---|---|---|
| `start_symptomatic`, `time_infectious`, `time_symptomatic`, `time_asymptomatic` | natural history of the pathogen | hand-curated, stable, in git |
| `<transmission>` `b0`/`b1`/`b2` | a **derived** artefact | generated, never hand-edited |

The fit is not a property of the disease. It is a property of
**(disease x population x contact matrix x contact rule x fit form)**.

Consequently a disease file carrying a single fit can be used with a population it was
never fitted against, and nothing detects it. `disease_measles_usa.xml` is labelled as
fitted against a WI-Dane population; both WI-Dane and TX-Gaines populations appear in
the branch history. The `min` -> `mean` change widens this hazard but does not create
it — it exists today.

### 8.3 Proposed layout

```
main/resources/data/
  disease/
    measles_usa.xml                                  # natural history only
  calibration/
    measles_usa__dane-wi-c1000__min__linear.xml
    measles_usa__gaines-tx-c1000__mean__linear.xml
    INDEX.csv                                        # generated manifest
```

The calibration key encodes every dimension the fit depends on, including the contact
rule and the fit form, so that restoring the quadratic term (F10) cannot silently change
the meaning of an existing artefact.

### 8.4 Calibration artefact format

```xml
<calibration>
  <transmission>
    <b0>...</b0><b1>...</b1><b2>...</b2>
  </transmission>
  <provenance>
    <population_file_sha256>...</population_file_sha256>
    <age_contact_matrix_file_sha256>...</age_contact_matrix_file_sha256>
    <disease_file_sha256>...</disease_file_sha256>
    <contact_probability_rule>mean</contact_probability_rule>
    <fit_form>linear</fit_form>
    <stride_commit>...</stride_commit>
    <r_squared>...</r_squared>
    <fit_r0_range>0-24</fit_r0_range>
    <num_runs>3645</num_runs>
    <date>...</date>
  </provenance>
</calibration>
```

**Content hashes rather than paths** is the specific remedy for 8.1.3. A path breaks when
a file is moved, re-exported or stranded in a superseded install root; a SHA-256 remains
verifiable indefinitely and identifies the input unambiguously.

### 8.5 Enforcement

`run_rStride()` already applies a validation gate before any experiment launches
(`.rstride$data_files_exist`, `log_levels_exist`, `valid_r0_values`,
`valid_immunity_profiles`, `valid_seed_infected`, `valid_cnt_param`). Add
`.rstride$valid_calibration()` to it:

1. Resolve the calibration artefact for the run.
2. Recompute SHA-256 of the population file, contact matrix and disease file actually
   configured.
3. Compare against the recorded provenance, including the contact rule.
4. **Abort on mismatch.**

This converts a silent wrong-calibration into an impossible one. `SimBuilder.cpp:103` is
the single equivalent hook on the C++ side if the same check is wanted in the kernel.

### 8.6 Promotion

```r
promote_calibration(project_dir, force = FALSE)
```

Responsibilities:

- validate fit quality — R-squared, coverage of the requested R0 range, convergence;
- compute the calibration key and the input hashes;
- write the artefact into the **repository** path, taken from configuration or an
  environment variable, never the install root;
- refuse to overwrite an existing key unless `force = TRUE`, and when forced, record the
  hash of the superseded artefact.

The result of a promotion is a git diff, so promotion becomes a pull request: reviewed,
attributed and permanently recorded. The non-overwrite guarantee is then enforced by
tooling and version control rather than by manual discipline.

### 8.7 Migration to `mean`

See Phase 2b for the sequenced steps. The property that matters: because calibration keys
differ by rule, `min` and `mean` calibrations coexist in version control. There is never
a window in which the model is uncalibrated, the two rules can be compared on identical
populations, and the change of default is a single reviewable diff rather than ten files
quietly mutating.

### 8.8 Lighter first step

If separating the files is too invasive at the point this is picked up, the following
obtains most of the safety for roughly a day's work and is forward-compatible with §8.3:

- keep a single disease file, adding a `<provenance>` block with hashes alongside the
  existing `<label>`;
- implement `promote_calibration()` with the repository target and non-overwrite
  semantics of §8.6;
- add the §8.5 hash check to the validation gate — warning first, aborting once trusted.

Separating natural history from calibration can then follow in Phase 6, when those files
are being touched anyway.

---

## 9. Continuous integration

### 9.1 The two suites are complementary — keep both

| | C++ gtester | rStride regression suite |
|---|---|---|
| Location | `test/cpp/gtester/` | `main/r/rStride_gtester_covid19.R` |
| Scope | kernel only | **whole pipeline**: config generation -> stride -> log parsing -> aggregation |
| Scenarios | 22 (19 covid, 3 influenza) x 2 thread settings | 23 scenarios x 5 seeds = 115 runs |
| Assertion | `num_cases` within a **margin** | **byte-exact** across 6 output streams |
| Cost | fast | **1336 s** of simulation time |

Porting the R suite into C++ would be a mistake. Its value is that it exercises `rStride`
— 13,668 lines with no other tests — and its assertions depend on R's own parsing and
aggregation, which are themselves the code under test. The two suites sit at different
levels with different sensitivities: a fast tolerance-based kernel smoke test, and a slow
exact full-pipeline regression.

### 9.2 Tiering

| Trigger | Job | Target |
|---|---|---|
| Every pull request | build gcc + clang, run C++ gtester | ~8 min |
| Every pull request | `R CMD check` on the rStride package (after Phase 1) | ~3 min cached |
| **Nightly on `master`** | full rStride regression suite | ~25 min |
| Weekly | macOS build + gtester | — |

The regression suite must not gate pull requests: twenty-plus minutes in front of every
change trains people to route around it. Nightly on `master` detects regressions within a
day, which is the right resolution here.

### 9.3 R runs in GitHub Actions without difficulty

`r-lib/actions` is mature. The obstacle is not R but the dependency load:
`smd_load_packages()` pulls twenty packages including `sf`, `tigris` and `usmap`, which
need GDAL/GEOS/PROJ and take roughly ten minutes to build from source. Two mitigations,
both already in this plan: Phase 1 step 3 moves the USA-only packages to `Suggests`, and
Posit Package Manager serves precompiled Linux binaries.

Note that the present suite is **not** testthat — it is a script with hand-rolled
comparison. Phase 1 makes `testthat` available; Phase 4 step 10 extracts the comparison
engine, which is what makes it callable from CI.

### 9.4 Prerequisites, and one genuine design decision

1. **`-ffast-math` must go first (F13.1, Phase 3b).** CI runs Linux/gcc while development
   is macOS/AppleClang. With that flag the two legitimately produce different
   floating-point results, so an exact-comparison job could never pass reliably.

2. **Removing it is necessary but not sufficient — this is the decision.** Bit-exact
   floating-point reproducibility *across platforms* is not attainable even under strict
   IEEE, because `std::exp` differs between glibc and macOS libm and the transmission path
   is dense with `exp()`. Two options:

   - **Pin the regression suite to one platform** (`ubuntu-latest`) and treat it as the
     source of truth. Local macOS runs then cannot be compared byte-exactly against CI
     references.
   - **Add a tolerance to the floating-point comparison.** The engine already computes the
     order of magnitude of differences for FP columns but does not use it for pass/fail,
     which is exact (`colSums(project_output != reference_output)`). Introducing a
     threshold there makes the suite portable.

   The second is preferable — it makes the suite meaningful on any machine and promotes
   the magnitude reporting from diagnostic to mechanism — but it slightly weakens the
   guarantee, so it is a deliberate choice rather than an obvious one.

3. **The suite must return an exit code.** `rStride_gtester_covid19.R` currently only
   prints; CI cannot detect a regression from it. Part of Phase 4 step 10.

4. **Reference consistency (F6).** The six `.rds` files span two dates six weeks apart.
   CI built on them would report unattributable failures.

5. **Pin `simid.rtools`** to a tag or commit rather than installing from GitHub HEAD, or a
   change in that repository can silently alter a regression result.

### 9.5 Drafted workflows

`.github/workflows/ci.yml` and `.github/workflows/regression.yml` are drafted and carry
their prerequisites in header comments. Both use `fetch-depth: 1`, which reduces a
checkout to HEAD blobs and makes F14's history weight irrelevant to CI.

Two notes on the drafts. Neither requires Boost — `util/Ptree.h` is a pugixml-based
replacement for `boost::property_tree`, and trng, pugixml, spdlog and tclap are vendored.
And both install `libomp-dev`, so **OpenMP is available in CI although absent locally**
(F13.2); the gtester's second instantiation runs with
`ConfigInfo::NumberAvailableThreads()`, so CI will genuinely exercise the multi-threaded
path that local builds currently cannot. That is additional coverage, but it also means CI
can fail on something never seen locally — which is a reason to fix local OpenMP
detection, not a reason to disable it in CI.
