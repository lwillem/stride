# Removed feature: Python / MDP interface

**Status:** specified for removal; not yet removed.
**Reference implementation:** tag `pre-refactor-2026-09` (see §6.4 of
`../rStride_architecture_and_refactoring.md`).
**Removal commit:** _to be filled in once the removal lands._
**Specified:** 2026-09-25

This document specifies a feature deliberately removed from Stride so that it can be
re-implemented later without archaeology. It is a pointer plus semantics, not an archive
— the code itself is preserved by git at the tag above.

---

## 1. Summary

A Python binding layer exposing Stride as a **step-wise, controllable simulator** for
reinforcement-learning agents. Rather than running a simulation to completion from a
configuration file, it lets an external agent advance the simulation one day at a time,
inspect the population state, and apply vaccination actions between steps.

The intended application was bandit / RL selection of vaccination strategies: which age
group to vaccinate with which vaccine type, re-decided every N days.

---

## 2. Architecture

Three layers, all part of the same feature:

```
main/r/rStride_MDP.R          generates a run configuration XML for the workflow
        │                     (does NOT invoke the simulator itself)
        ▼
main/python/pybind/           pybind11 module `pylibstride`
        │                     exposes MDP, AgeGroup, VaccineType to Python
        ▼
main/cpp/mdp/                 MDP / MDPRunner: step-wise simulation control API
        │                     wraps the normal Sim with day-level stepping + vaccination
        ▼
libstride                     the ordinary Stride core (unchanged)
```

The Python agent loop is illustrated by `main/python/pybind/test.py`:

```python
mdp = stride.MDP()
mdp.Create(default_file)
for i in range(numDays // step_size):
    age_group, vaccine_type = choose_action()        # bandit decision
    for _ in range(step_size):
        mdp.Vaccinate(availableVaccines=10000, ageGroup=age_group, vaccineType=vaccine_type)
        infected = mdp.SimulateDay()
mdp.End()
```

---

## 3. Current state: already dormant

**The feature is not built.** `main/python/CMakeLists.txt` contains:

```cmake
# exclude for now
#add_subdirectory(pybind)
```

**The MDP core is not compiled either.** `main/cpp/CMakeLists.txt` defines `STRIDE_SRC`
as an explicit source list, and no `mdp/*.cpp` file appears in it. The MDP sources are
compiled *only* by the pybind target, via `PYLIBSTRIDE_SRC` in
`main/python/pybind/CMakeLists.txt`, which is itself disabled.

Consequently this code is **dead weight in the source tree and absent from every build
artefact**. Removing it changes no binary and no simulation result.

### 3.1 Why it was disabled, and what re-enabling would actually require

The subdirectory was commented out because the build failed when a suitable Python
version was not available. **Re-enabling it is not simply a matter of uncommenting the
line**, and the reasons differ from what one might assume. Verified 2026-09-25:

**What is *not* broken:**

| Check | Result |
|---|---|
| `main/cpp/mdp/MDP.cpp` compiles against the current core | clean, 0 errors |
| `main/cpp/mdp/MDPRunner.cpp` compiles against the current core | clean, 0 errors |
| `STRIDE_WC_REVISION_HASH`, `PROCCOUNT`, `PYTHON_LIBRARY_DIR` defined outside `pybind/` | yes, all three |
| `libstride` target exists for `target_link_libraries` | yes (`main/cpp/CMakeLists.txt:69`) |
| pybind wrappers consistent with `MDP.h` | yes, including `VaccineProperties` subclasses |

The C++ has therefore **not** drifted away from the core API despite not being compiled
by the main build since 2024. That is the part most likely to have rotted, and it has not.

**What *is* broken:**

1. **`test.py` is stale.** It calls `mdp.Create(default_file)` with a single argument,
   but `MDP::Create` (`MDP.h:59-63`) requires `configPath`, `mRNA_properties` and
   `adeno_properties`. The example was never updated when vaccine properties were added.
   A working call needs something of the form:

   ```python
   mrna  = stride.ConstantVaccineProperties("mRNA",  ve_susceptible, ve_infectiousness, ve_severe)
   adeno = stride.ConstantVaccineProperties("adeno", ve_susceptible, ve_infectiousness, ve_severe)
   mdp.Create(default_file, mrna, adeno)
   ```

2. **The pybind11 pin is the real blocker.** `pybind11.cmake` fetches
   `v2.6.0` (October 2020). Consequences on a current toolchain:

   - it predates Python 3.12 support — 3.12 removed `distutils` and changed C-API
     internals that pybind11 accommodated only from roughly 2.11/2.12;
   - its own CMake declares compatibility with CMake < 3.5, which CMake 4.0 rejects
     outright;
   - `FetchContent_Populate(pybind11)` uses the single-argument form, deprecated in
     CMake 3.30 and removed in CMake 4.0;
   - it requires network access at configure time.

   At the time of writing the development machine ran **CMake 3.24.4 and Python 3.8.3**,
   both inside v2.6.0's supported window — which is why the failure appeared only when a
   different interpreter was picked up, rather than unconditionally.

**Therefore a re-implementation should not attempt to revive the pinned build.** The
minimum work is: bump pybind11 to a current release, adapt to whatever binding-API
changes that entails, modernise the `FetchContent` call, and rewrite `test.py` against
the actual `Create` signature. The C++ core underneath can be taken as-is.

---

## 4. Inventory

| Path | Lines | Role |
|---|---:|---|
| `main/cpp/mdp/MDP.cpp` | 473 | step-wise simulation control, vaccination logic |
| `main/cpp/mdp/MDP.h` | 156 | public API (see §5) |
| `main/cpp/mdp/Vaccines.h` | 93 | vaccine type definitions |
| `main/cpp/mdp/AgeGroup.h` | 88 | age-group enumerations |
| `main/cpp/mdp/MDPRunner.h` | 51 | runner |
| `main/cpp/mdp/MDPRunner.cpp` | 37 | runner |
| `main/python/pybind/mdp/MDP.cpp` | 73 | pybind wrapper for MDP |
| `main/python/pybind/mdp/AgeGroup.cpp` | 33 | pybind wrapper |
| `main/python/pybind/mdp/Vaccines.cpp` | 30 | pybind wrapper |
| `main/python/pybind/libstride.cpp` | 19 | module entry point |
| `main/python/pybind/CMakeLists.txt` | 55 | pybind build |
| `main/python/pybind/pybind11.cmake` | 15 | dependency fetch |
| `main/python/pybind/INSTALL.txt` | 33 | setup notes |
| `main/python/pybind/test.py` | 44 | example agent loop |
| `main/python/pybind/py_test/run_default.xml` | — | test configuration |
| `main/python/CMakeLists.txt` | 26 | currently only disables the above |
| `main/r/rStride_MDP.R` | 108 | configuration generator (installed; see §7) |

Approximately **1,200 lines** in total.

---

## 5. Public API surface (`main/cpp/mdp/MDP.h`)

The contract a re-implementation must restore:

**Lifecycle**

- `Create(configPath, seed, outputDir, outputPrefix, childless, uptake)`
- `End()`
- `ClearSimulation()`

**Stepping**

- `SimulateDay() -> unsigned int` — advance one day, return infected count
- `Simulate(numDays) -> unsigned int`
- `SimulateVaccinate(numDays, availableVaccines, ...)`

**Actions**

- `Vaccinate(availableVaccines, AgeGroup, VaccineType)`
- `VaccinateChildless(availableVaccines, ChildlessAgeGroup, VaccineType)`
- `UpdateCntReduction(workplace_distancing, community_distancing, ...)`

**Observation**

- `GetNumberOfDays()`, `GetPopulationSize()`, `GetAtRisk()`
- `GetAgeGroupSizes()`, `GetChildlessAgeGroupSizes()`
- `GetVaccinatedAgeGroups()`, `GetVaccinatedChildlessAgeGroups()`
- `GetTotalInfected()`, `GetTotalHospitalised()`
- `CountInfectedCases()`, `CountExposedCases()`, `CountInfectiousCases()`,
  `CountSymptomaticCases()`, `CountHospitalisedCases()`

The Python module exposes three names: `stride.MDP`, `stride.AllAgeGroups`,
`stride.AllVaccineTypes`.

---

## 6. Residual coupling into the live core

One piece of MDP support exists in code that **is** compiled into `libstride`:

```cpp
/// Added for MDP memory management: clear the contact pools
void ClearContactPools();
```

- declared at `main/cpp/contact/ContactPoolSys.h:68`
- defined at `main/cpp/contact/ContactPoolSys.cpp:41`
- **sole caller:** `main/cpp/mdp/MDP.cpp:454`

After removal this becomes dead code in the shipping binary. Either delete it in the same
commit and record that here, or retain it deliberately and note why. It must not be left
undecided.

---

## 7. `rStride_MDP.R` — note before removing

Unlike the other pieces, this file **is installed** (it is listed in
`main/r/CMakeLists.txt`). It does not invoke the simulator: it builds a configuration,
creates a calendar, and calls `save_config_xml()`. Its output was consumed by the Python
agent loop.

It carries workflow-specific assumptions that will not survive relocation:

- `exp_param_list$immunity_distribution_file <- "../FullPop/immunity_covid_belgium.xml"`
- `exp_param_list$output_prefix <- "config/vsc_1/"` and
  `config_exp$output_prefix <- "runs/vcs_1/"` — VSC cluster paths, and note the
  `vsc`/`vcs` inconsistency between the two
- `num_threads <- 16`

It should be removed together with the rest and its install entry dropped from
`main/r/CMakeLists.txt`.

---

## 8. Rot found at time of specification

Recorded because it indicates how long the feature has been unmaintained, and because a
re-implementation should not reproduce it:

1. **`mdp/Health.cpp` does not exist**, yet is referenced in two places as commented-out
   code: `libstride.cpp` (`//void init_health(py::module &);` and `//init_health(m);`)
   and `pybind/CMakeLists.txt` (`#    "mdp/Health.cpp"`). A health-observation binding was
   evidently planned or deleted.
2. **`test.py` carries unresolved TODOs**, including
   `# TODO: use same seed in python and c++?` — meaning reproducibility between the
   Python driver and the C++ core was never established.
3. **`test.py` does not match the bound API** — see §3.1.
4. **The pybind11 dependency is pinned to a 2020 release** — see §3.1.
5. Pybind sources are dated 2024; `main/python/CMakeLists.txt` carrying the disabling
   comment is dated 2026.

---

## 9. Acceptance test for a re-implementation

The defining property is **equivalence with the ordinary simulator when no action is
taken**:

> For a given configuration and seed, stepping the MDP interface `N` times with no
> vaccination must produce the same epidemic trajectory as running the standard
> `stride -c <config>` binary for `N` days.

This is the check the original never established (§8.2) and is the single most valuable
thing to build first. Supporting checks:

- `GetPopulationSize()` equals the population file record count.
- `sum(GetAgeGroupSizes().values()) == GetPopulationSize()`.
- After `Vaccinate(n, group, type)`, `GetVaccinatedAgeGroups()[group]` increases by
  `min(n, eligible_in_group)`.
- `ClearSimulation()` followed by `Create()` reproduces a fresh run bit-identically —
  the memory-management path that `ClearContactPools()` exists for.

---

## 10. Removal requirements

1. **One self-contained commit**, so a single `git revert` restores the feature.
2. **Delete, do not comment out.** The current `#add_subdirectory(pybind)` is exactly the
   failure mode this specification replaces: disabled-in-place code that decays silently.
   Remove `main/python/` entirely, remove `add_subdirectory(python)` from
   `main/CMakeLists.txt`, remove `main/cpp/mdp/`, and remove `rStride_MDP.R` with its
   install entry.
3. **Decide `ClearContactPools()`** per §6 and record the decision here.
4. **No tombstone is required.** Unlike `track_index_case`, this feature has no
   configuration key that could silently change behaviour — it is a separate build target
   that simply ceases to exist.
5. **Prove behaviour preservation.** The regression reference `.rds` files must be
   byte-identical. Because nothing in `main/cpp/mdp/` or `main/python/` is compiled into
   the current binary, this is guaranteed by construction and should be confirmed rather
   than assumed.
6. **Sequence after the branch consolidation** of §6.4, consistent with the other
   excision.

---

## 11. Recovering the original implementation

```sh
git show pre-refactor-2026-09 --stat -- main/python main/cpp/mdp
git checkout pre-refactor-2026-09 -- main/python main/cpp/mdp main/r/rStride_MDP.R

# or, once the removal commit is recorded above:
git revert <removal-commit>
```

Re-enabling additionally requires restoring `add_subdirectory(python)` in
`main/CMakeLists.txt`, uncommenting `add_subdirectory(pybind)` in
`main/python/CMakeLists.txt`, and re-adding `rStride_MDP.R` to `main/r/CMakeLists.txt`.
