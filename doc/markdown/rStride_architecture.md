# rStride and Stride: As-Built Architecture

**Date:** 2026-10-05
**Baseline:** `master` at `pre-refactor-2026-10`
**Scope:** how the system works today — the R workbench (`main/r/`), the C++ kernel
(`main/cpp/`), and the build and test machinery that connects them.

This document is **descriptive**. It records how the system is built, not what should
change about it. Every problem, judgement and proposed remedy lives in
`rStride_refactoring_plan.md`, which references the sections here rather than restating
them. Where a structure described below is known to be a problem, the finding that says so
is cited inline, like this: *(see F12 in the refactoring plan)*.

---


## 1. Current flow (as built)

### 1.1 From repository to a runnable workbench

```
repo/stride_2026/
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
- 16 `if(0==1){ attach(...) }` interactive-debug blocks (two distinct patterns; see F2.1)

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

---

## 3. Contact pools and the venue extension

### 3.1 The two families of contact pool

Four contact-pool types were added alongside the original set: `OtherHouse`, `RestoCafe`,
`OtherPlace` and `Transport`. They are structurally different from the classic pools:

| | Classic pools | New venues |
|---|---|---|
| Pool id per person | one (`m_pool_ids[type][0]`) | one **per day of week** |
| Source | population CSV column | separate `subpools_community_file` |
| Contact rate | `AgeContactProfile`, by age | per-person `CPoolContacts(type)[day]` |
| Duration | not tracked | tracked per day |

*(The cost of expressing that difference by naming the four types at each decision point
is F12 in the refactoring plan.)*


### 3.2 Which predicate applies where

The literal four-type list appears **nine times across five files**:

```
PopBuilder.cpp:155,165,188,240   SimBuilder.cpp:93   Sim.cpp:141   Infector.cpp:222
Person.cpp:106-109,126-134       ContactDivider.cpp:48-61,117-120 (unrolled by hand)
```

Those sites encode **three different predicates that do not coincide**:

| Predicate | Sites | Members |
|---|---|---|
| loaded from subpool file / day-indexed | PopBuilder x4, `SimBuilder.cpp:93`, `Sim.cpp:141`, `Infector.cpp:222` | the 4 venues |
| has individual contact variation | `Infector.cpp:275` | 4 venues **+ Workplace + both Community**, **minus OtherHouse** |
| presence toggled per day | `Person.cpp:106-134` | the 4 venues |
| has physical venue characteristics (airborne) | `PoolCharacteristicsSeeder.cpp:94` | **School, Workplace, Collectivity + the 4 venues** |

`Infector.cpp:275` and `PoolCharacteristicsSeeder.cpp:94` are genuinely different sets.
*(F12.1.)*


### 3.3 The population file format

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

### 3.4 `Person` memory layout

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

### 3.5 What is already generalized, and the three day-gating mechanisms

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
---

## 4. Transmission, R0 and calibration

Each disease configuration carries a fitted regression of R0 on the transmission
probability:

```xml
<transmission>
  <b0>0.820990685751851</b0>
  <b1>18.91204185976</b1>
  <b2>0</b2>
</transmission>
```

`TransmissionProfile.cpp:57-77` inverts `E(R0) = b0 + b1*p + b2*p^2` to derive the
transmission probability from a requested R0; `TransmissionProfile.cpp:62` branches on
`b2 == 0` to choose between the linear and quadratic inversion. All ten disease files in
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

A fit is therefore not a property of the disease alone: it is a property of
**(disease x population x contact matrix x contact rule x fit form)**. Nothing in the
configuration records which population a fit was made against, and nothing checks it.
*(F10, and §8 of the refactoring plan, which specifies calibration artefacts.)*

The fit is produced by `bin/rStride_r0.R`, which runs an R0 sweep and calls
`analyse_transmission_data_for_r0()` (`TransmissionAnalyst.R`). That function regresses
secondary cases on the transmission probability, writes the coefficients into the parsed
disease configuration together with a `<label>` provenance block, and saves the result as
`sim_output/<timestamp>_r0/<run_tag>_disease_<name>.xml` — **it never overwrites the
original**, so promotion into `main/resources/data/` is a manual copy.
`TransmissionAnalyst.R:96` currently hardcodes `fit_b2 <- 0`, so new fits are linear even
where a committed file carries a quadratic one.

### 4.1 The contact-probability rule

`Infector.cpp :: GetContactProbability` combines the two age-specific contact
probabilities of a candidate pair. Since 2026-09-28 the rule is a configuration value:

```
run.contact_probability_rule = Min (default) | Mean
```

read once in `SimBuilder`, stored on `Sim`, and threaded through `InfectorExec` and both
`Infector::Exec` bodies. `Min` is the rule every committed disease-file fit was produced
under. An unknown value throws rather than defaulting silently.

---

## 5. The Infector dispatch

The transmission loop is dispatched in four stages:

| Stage | Location | Frequency |
|---|---|---|
| `InfectorMap` lookup resolved to an `InfectorExec*` | `SimBuilder.cpp:72` | once, at construction |
| Selection between default and tracing infector | `Sim.cpp:78` | once per timestep |
| Indirect call into `Infector<LL,TIC,TO>::Exec` | `Sim::TimeStep()` | once per contact pool |
| Contact evaluation loop | inside `Exec` | per candidate contact pair |

The indirection is amortised over an entire contact pool and is not a per-contact cost.

`LOG_POLICY<LL>` (`Infector.cpp:101-201`) is the part that earns the design. For
`EventLogMode::None` its `Contact()` is an empty inline function, so in the innermost loop
it compiles to nothing at all, including the argument setup.

There are two `Exec` bodies: a general form (`Infector.cpp:330-498`, 169 lines) and an
optimised form (`504-651`, 148 lines), selected by the `TO` template parameter.
Whitespace-normalised, only 84 of 317 lines differ. Ten explicit instantiations exist,
`TIC` doubling the count. *(Phase 7a of the refactoring plan.)*

---

## 6. Build configuration as it stands

The flags reaching the compiler in a Release build:

```
-g -fvisibility=hidden -std=c++17 -Wall -Wextra -pedantic -Weffc++ \
-Wno-unknown-pragmas -std=c++1z -O3 -DNDEBUG -arch arm64 ...
```

`-ffast-math` was removed on 2026-10-05 (F13.1). `-g` in Release, the duplicated `-std=`,
and the never-set `CMAKE_CXX_STANDARD` remain. *(F13.3.)*

### 6.1 OpenMP is not linked on the development machine

| | |
|---|---|
| `~/opt/stride-*/bin/stride` | Mach-O **arm64**; `otool -L` shows only `libSystem`, `libc++` |
| `/usr/local/Cellar/libomp/.../libomp.dylib` | Mach-O **x86_64** |
| `/opt/homebrew` (Apple Silicon prefix) | does not exist |
| `CMakeCache.txt` | `HAVE_FOUND_OpenMP:BOOL=FALSE`, `OpenMP_CXX_FLAGS:STRING=NOTFOUND` |

An arm64 binary cannot link an x86_64 library, so detection fails and the build falls back
to the dummy-OpenMP stubs in `main/resources/lib/domp`. Every `#pragma omp parallel`
compiles to nothing. `rStride.R` hardcodes `num_threads <- 1`, so the runtime impact is
nil; the testing impact is not. *(F13.2.)*

### 6.2 Resource installation

`main/resources/CMakeLists.txt` unpacks `data/*.zip` at configure time and installs the
contents alongside the non-zipped files, the config directory and the regression
references. `main/r/CMakeLists.txt` installs `rStride_*.R` into `bin/` and the `rstride`
library directory wholesale. Both use `CONFIGURE_DEPENDS` globs, so adding a script or an
archive requires no edit. *(This was not always so — F7, F7.1.)*

---

## 7. The two test suites

| | C++ gtester | rStride regression suite |
|---|---|---|
| Location | `test/cpp/gtester/` | `main/r/rStride_gtester_covid19.R` |
| Scope | kernel only | **whole pipeline**: config generation -> stride -> log parsing -> aggregation |
| Scenarios | 22 (19 covid, 3 influenza) x 2 thread settings | 23 scenarios x 5 seeds = 115 runs |
| Assertion | `num_cases` within a **margin** | **byte-exact** across 6 output streams |
| Cost | ~91 s | ~6 min wall-clock, ~360 s of simulation |

The six reference `.rds` files live at `main/resources/rstride_test/` and are installed to
`tests/`. `rrv()` rewrites them from a completed run; `rrv_repo()` writes both the install
copy and the repository copy. They were reset from a single run on 2026-10-05.

The gtester instantiates each scenario twice, once with one thread and once with
`ConfigInfo::NumberAvailableThreads()`. Without OpenMP (§6.1) both are single-threaded, so
22 of the 44 instances are duplicates.

---

## 8. Repository composition


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
