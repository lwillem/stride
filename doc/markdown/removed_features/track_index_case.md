# Removed feature: `track_index_case` (TIC)

**Status:** specified for removal; not yet removed.
**Reference implementation:** tag `pre-refactor-2026-09` (see §6.4 of
`../rStride_architecture_and_refactoring.md`).
**Removal commit:** _to be filled in once the removal lands._
**Specified:** 2026-09-25

This document is the specification of a feature deliberately removed from Stride. It
exists so the feature can be re-implemented later without archaeology. It is a pointer
plus semantics, not an archive — the code itself is preserved by git at the tag above.

---

## 1. Summary

`track_index_case` is a run-configuration flag that suppresses onward transmission from
every newly infected person. Cases are still created, attributed and logged, but they
never become infectious. The effect is that each index case's infections can be counted
without competition or susceptible depletion caused by the secondary cases themselves.

It is a **measurement instrument**, not an epidemiological model feature: it makes the
simulator produce a clean offspring count for the first generation.

---

## 2. Why it existed

TIC isolates the first generation of transmission. By marking every newly infected person
recovered immediately, it removes competition between secondary cases for the remaining
susceptibles, so the offspring count of each index case can be observed in an
undepleted population.

### 2.1 It is **not** needed for R0 calibration — do not reintroduce it for that purpose

This is recorded explicitly because it is a plausible but incorrect future use.

The `<transmission>` `b0`/`b1`/`b2` block is **a calibration mapping, not an estimate of a
theoretical quantity**. `analyse_transmission_data_for_r0()`
(`main/r/rstride/TransmissionAnalyst.R:47-49`) fits `E(R0) = b0 + b1*p` by counting
offspring for every infector at every generation:

```r
tbl_infections <- table(data_transm$infector_id)
```

and `TransmissionProfile.cpp:62-77` inverts that relation so a requested `r0` yields the
corresponding transmission probability.

The fit must therefore be derived under **the same conditions under which it is used** —
a natural flow of infection in which secondary cases are themselves infectious and
compete for the remaining susceptibles. Calibrating under TIC, with chains truncated at
generation one and competition suppressed, and then applying the resulting curve to
ordinary runs would introduce precisely the mismatch the calibration exists to prevent.

This matters more in Stride than in a homogeneous-mixing model. Competition for
susceptibles within households, school classes and workplaces is a first-order mechanism
of a structured population, not a nuisance to be removed: in a four-person household,
secondary cases competing for the last susceptible member is real model behaviour. A
TIC-derived curve would systematically overstate transmission relative to how the model
actually runs.

Depletion is bounded instead by a short horizon — `rStride_r0.R` uses `num_days = 20`.

### 2.2 What TIC would legitimately be for

Measuring the first-generation offspring distribution as a quantity in its own right —
for example, comparing the model's realised secondary-case distribution against
observed data on index-case transmission, or characterising individual-level
heterogeneity and overdispersion without the confounding of chain depletion. Those are
measurement exercises, distinct from parameter calibration.

---

## 3. Semantics

At **every point where a newly infected person's infection is started** in the
transmission loops, when the flag is enabled, immediately terminate that infection:

```cpp
hX.StartInfection(id_index_case, id_infector, rel_inf);
if (track_index_case)
        hX.StopInfection();
LP::Trans(eventLogger, infector, infectee, pType, simDay, id_index_case, ...);
```

Required properties:

1. **The case is still registered.** `StartInfection()` records the index-case id and the
   infector id. The transmission event is still emitted to the event log. The case
   therefore appears in `data_transmission` and in incidence output exactly as it would
   otherwise.
2. **The case never becomes infectious.** `Health::StopInfection()`
   (`main/cpp/health/Health.cpp:74`) sets `m_status = HealthStatus::Recovered` and clears
   `m_hospitalised`. Because `IsInfectious()` and `IsSusceptible()` are both false for
   `Recovered`, the person can neither transmit nor be re-infected.
3. **Ordering matters.** `StopInfection()` begins with
   `AssertThrow(IsInfected(), ...)`, so it must follow `StartInfection()` directly. The
   logging call may come after; it does not depend on health status.

### Health status semantics relied upon

| Status | Value | Relevant predicate |
|---|---:|---|
| `Susceptible` | 0 | `IsSusceptible()` true |
| `Exposed` | 1 | |
| `Infectious` | 2 | `IsInfectious()` true |
| `Symptomatic` | 3 | |
| `InfectiousAndSymptomatic` | 4 | `IsInfectious()` true |
| `Recovered` | 5 | `IsRecovered()` true; not susceptible, not infectious |

---

## 4. Configuration surface

| Location | Content |
|---|---|
| `main/cpp/sim/SimBuilder.cpp:55` | `m_config.get<bool>("run.track_index_case")` — **no default**, so the key was mandatory |
| `main/cpp/util/RunConfigManager.cpp` | lines 106, 137, 172 — three embedded default configurations, all `false` |
| `main/resources/config/run_default.xml:26` | `false` |
| `main/resources/config/run_default_air.xml:30` | `false` |
| `main/r/rstride/rStride.R:95` | `config_default$track_index_case <- 'false'` — hardcoded |

**It was never enabled from the R workbench.** No script in `main/r`, `test` or `doc`
overrode the hardcoded `'false'`. It remained reachable only by invoking
`./bin/stride -c <config>` directly with the key set to `true`.

---

## 5. Implementation touch points

Line numbers are those of tag `pre-refactor-2026-09` and will drift; the semantic
descriptions are authoritative.

### 5.1 The five infection sites (`main/cpp/contact/Infector.cpp`)

| Line | Context |
|---:|---|
| 406 | general `Exec`, pairwise transmission, p1 infectious -> p2 |
| 418 | general `Exec`, pairwise transmission, p2 infectious -> p1 |
| 485 | general `Exec`, airborne / subpool transmission path |
| 573 | optimized `Exec`, pairwise transmission |
| 638 | optimized `Exec`, airborne / subpool transmission path |

Semantically: **every call site of `Health::StartInfection()` inside the two `Exec`
bodies.** Any re-implementation must cover all of them, including both the pairwise and
the airborne paths, in both the general and the optimized template specialisation.

### 5.2 Template dispatch

| File | Role |
|---|---|
| `contact/Infector.h` | `TIC` is the second template parameter of `Infector<LL, TIC, TO>`; ten `extern template` declarations (5 log modes x 2 boolean values) |
| `contact/Infector.cpp` | matching ten explicit instantiations (lines 656-665); two `Exec` definitions templated on `TIC` |
| `contact/InfectorMap.h` | map keyed on `std::tuple<EventLogMode::Id, bool>`; `Add<bool B>()` called for `true` and `false` |
| `sim/Sim.h:87`, `sim/Sim.cpp:42` | `bool m_track_index_case` member |
| `sim/SimBuilder.cpp` | 55 (read), 71, 76, 80, 83 (four tuple constructions for map lookup) |

Note: the tuple site at `SimBuilder.cpp:76` lies inside a branch testing
`event_log_level == "ContactTracing"`, but `EventLogMode::Id` defines no such value
(only `None`, `Incidence`, `Transmissions`, `Participants`, `All`). That branch appears
unreachable and should be verified separately during removal; it is not part of this
feature.

---

## 6. Acceptance test for a re-implementation

This is the property the feature exists to guarantee, and it is mechanically checkable
from ordinary run output:

> With `track_index_case = true`, no `infector_id` appearing in `data_transmission` may
> belong to a person who is not an index case.

Equivalently, in R:

```r
# every infector must itself have been seeded, not infected during the run
stopifnot(all(!data_transm$infector_id %in% data_transm$local_id))
```

Supporting checks:

- With the flag enabled, the maximum generation depth in the transmission tree is 1.
- Total case count equals the number of index cases plus their direct offspring.
- With the flag disabled, all output must be **bit-identical** to the pre-feature build;
  this is what makes the feature safe to add or remove.

---

## 7. Known subtleties

1. **TIC does not control seeding.** It has no interaction with `DiseaseSeeder`;
   `m_track_index_case` was consumed only as an `InfectorMap` key. The number of index
   cases is set independently via `num_infected_seeds`. Using several seeds together with
   TIC yields several independent index cases, each truncated at generation one, which is
   usually the desired configuration.
2. **The source comments disagreed.** `Infector.cpp:395` said "no tertiary infections
   with TIC"; `Infector.cpp:572` said "No secondary infections with TIC". The code was
   identical at both sites and the first comment is the correct one: secondary cases *are*
   created and counted; the tertiary generation is what is prevented.
3. **`StopInfection()` does not reset the disease counter.** The call to
   `ResetDiseaseCounter()` is commented out in `Health::StopInfection()`. Any
   re-implementation should confirm whether that matters for the intended measurement
   before changing it.
4. **Hospitalisation is cleared.** `StopInfection()` sets `m_hospitalised = false`, so
   TIC also suppresses hospitalisation of secondary cases. If a future use needs burden
   estimates alongside offspring counts, this interaction must be revisited.

---

## 8. Simplification gained by removal

- `Infector<LL, TIC, TO>` becomes `Infector<LL, TO>`.
- Explicit template instantiations drop from **ten to five**.
- `InfectorMap`'s key simplifies from `std::tuple<EventLogMode::Id, bool>` to a plain
  `EventLogMode::Id`, and `Add<bool B>()` collapses to a single non-templated fill.
- `Sim::m_track_index_case` and its four lookup sites disappear.
- Five `if (TIC)` branches leave the transmission loops.

This supersedes step 1 of Phase 7a in the refactoring plan: the parameter is deleted
rather than converted to a runtime argument.

---

## 9. Removal requirements

1. **One self-contained commit**, so a single `git revert` restores the feature.
2. **Leave a tombstone.** After removal the key is no longer read, and Boost ptree
   silently ignores unknown keys — so an existing configuration setting
   `track_index_case = true` would run *without* the feature and produce quietly wrong
   results. Prevent this:

   ```cpp
   if (m_config.get_optional<bool>("run.track_index_case").value_or(false)) {
           throw std::runtime_error(
               "track_index_case was removed; see doc/markdown/removed_features/track_index_case.md");
   }
   ```

3. **Prove behaviour preservation.** The regression reference `.rds` files must be
   byte-identical after removal. Because no scenario enables the flag, this is a genuine
   proof rather than a formality, and satisfies the rule in §7.5 of the refactoring plan.
4. **Sequence after the branch consolidation** of §6.4. Removing a template parameter
   from `Infector.h` / `InfectorMap.h` beforehand would add conflict surface to a merge
   that is currently clean on the C++ side, and removing after the tag keeps
   `pre-refactor-2026-09` a faithful snapshot of the last version containing the feature.

---

## 10. Recovering the original implementation

```sh
git show pre-refactor-2026-09:main/cpp/contact/Infector.h
git show pre-refactor-2026-09:main/cpp/contact/Infector.cpp
git show pre-refactor-2026-09:main/cpp/contact/InfectorMap.h
git show pre-refactor-2026-09:main/cpp/sim/SimBuilder.cpp

# or, once the removal commit is recorded above:
git revert <removal-commit>
```
