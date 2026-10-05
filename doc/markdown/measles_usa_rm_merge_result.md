# `measles_usa_rm` — what happened to your changes in the merge

**Date:** 2026-10-05
**Branch:** `integration/measles-usa` (`ac987f8`)
**Merged:** `origin/measles_usa_rm` @ `2b927f4` (36 commits) into `measles_usa` (53 commits), from merge base `c2e209f` (2026-07-23)
**Decisions applied:** `measles_usa_rm_discussion.md`
**Status:** C++ gtester 22/22 in 107 s. R regression suite not yet run.

This is a record of how each part of your work was treated, so nothing has to be
discovered from a diff. **One file of yours was reverted** — `ImmunitySeeder.cpp` — and
that is the main thing to read here (§3).

---

## 1. Taken unchanged

Verified byte-identical to your branch in the merged tree:

| File | Lines | |
|---|---:|---|
| `main/cpp/util/PopSnapshotWriter.{h,cpp}` | +256 | plus the `SimController` call and the install entry |
| `main/r/rstride/MeaslesClustering.R` | +327 | |
| `main/r/rstride/factories/HouseholdClusterFactory_USA.R` | +329 | |
| `main/r/rstride/factories/ImmunityProfileFactory_USA.R` | +86 | |
| `main/resources/data/immunity_measles_dummy.xml` | +116 | |
| `main/r/rStride_contacts.R` | +10 | |

`4234e7d` (copy-constructor fix, `RnMan::Shuffle` on an index vector) is in — B8.

## 2. Taken, with adaptation

| Yours | What was done |
|---|---|
| `max_age` parameterisation of `get_cnt_data()` (`5e9bd33`) | **Adopted** (B4). Kept our `contactdata::` qualifier over the bare `contact_matrix()` to avoid colliding with `socialmixr`. Because the contact-data calls now read `max(pop_usa$age)`, the population load had to move *above* them; the later duplicate generator call was dropped. |
| CT `data_year = 2021` FIPS fix (`a53ced1`) | **Ported by hand** onto `factories/PopulationFactory_USA.R`. We had renamed `USA_PopulationBuilder.R` (C1), so git could not carry it across. The rename also removed the `source("~/Documents/Repositories/stride/...")` absolute path (C3). |
| Age-specific hospital probabilities (`2b927f4`) | **Adopted** into `rStride_measles_explore.R`: `hospital_probability_age = 0.18,0.06,0.12` over `hospital_category_age = 0,5,20`, `hospital_mean_delay_age = 2,2,5`. These replaced our `#TODO: set real value` placeholders. |
| `## THIS NEEDS TO STAY 'FALSE' TO GET HOUSEHOLD AGG` | **Kept** — the reasoning was worth preserving. |
| `MeaslesClustering.R` excluded from the library load | **Kept** in the blacklist. Your entries for `generate-population-community.R` and `USA_PopulationBuilder.R` were dropped because neither file exists any more. |

## 3. Reverted: `ImmunitySeeder`

**`a8f5095` and the pruning rework are NOT in the baseline.** This was not a merge
conflict — the file merged cleanly — it was reverted afterwards, deliberately.

The `influenza_c` gtester scenario (`ScenarioData.cpp:94`, `immunity_rate = 0.9991`)
became effectively non-terminating. Reverting only this file, with everything else in
the merge left in place, isolates it — same tree, same machine, same scenario filter:

| `ImmunitySeeder` | `influenza_c` | full suite |
|---|---|---|
| pre-merge version | **6,073 ms** | **22/22, 107 s** |
| `measles_usa_rm` | **>180,000 ms** (timed out) | 24 min, never completed |

**Why.** The sampling loop is structurally identical to the original except that it now
calls `isExhausted()` up to twice per draw, and each call is a linear scan of the
household. A pool is pruned only when *every* member is vaccinated or out of quota — so a
household holding one eligible member among vaccinated ones is never pruned, yet still
yields a rejected draw most of the time. At 99.91 % the rejection tail dominates, and
every rejected draw now pays two household scans it did not pay before.

So the pruning costs more than it saves in exactly the regime it was written for. That is
a property of the *mechanism*, not of the idea:

- The intent is right. `immunity_clustering_plan.md` §1.4 describes this pathology.
- The sampling distribution is preserved — your validation showed a household
  pair-correlation excess of −0.02 against the original's −0.02, across 8 seeds.
- Both versions terminate on the gtester's configurations. The old one is simply ~30×
  faster there, at 600k people and ~250k households, which is a larger population than
  these runs normally use.

**Suggested fix, for when the immunity code is picked up:** keep a per-pool count of
eligible members and decrement it on vaccination, so `isExhausted()` becomes an O(1) test
instead of a scan. That should beat both the original and the current version, and keeps
`a8f5095`'s intent. `immunity_clustering_plan.md` §4 Phase I step 3 asks for a seeding
benchmark at 50/70/90/95 % — that benchmark is what caught this, and it should be part of
the fix.

Nothing outside the file referenced the new API, so the revert left no dangling callers:
`RandomIndependent`, `run.mass_immunize` and `run.vaccine_hesitancy_rate` have no users
elsewhere in `main/cpp`.

## 4. Not taken, and why

| Yours | Reason |
|---|---|
| CT/Fairfield and Spartanburg SC targeting | The data they need exists on **neither** branch: `data/population_6region/`, `data/contacts_6region/`, `data/immunity_6region/`, `disease_measles_6region.xml`. Those files appear to be local to your machine. The study scripts stay on Gaines TX, which is in the tree. |
| `workplace_ages <- 18:69` | B2 — ours derives the band from the ages actually employed, because people over 65 are assigned to workplaces in the US population but got a zero conditional rate under the fixed band. |
| `num_people_workplace_leq7` exploration | B3 — replaced by the `uniroot` cluster-size adjustment factors, which are now written into the contact-matrix XML and read by `AgeContactProfile`. |
| `export = TRUE` and the `export` argument | C2 — `getFREDdata()` no longer takes it. |
| Removing the five geospatial packages from the load | C4 — you are right that this is where the refactoring is going (Phase 1 moves them to `Suggests:`), but doing it now would break the population generators, which are excluded from the library load and rely on them. Deferred, not rejected. |
| The vaccine-hesitancy / mass-immunisation config block | Kept in `rStride_measles_explore.R` **commented out**, because `vaccine_distribution_file` points at `data/immunity_6region/`. **The C++ side is merged and available** — see §5. |

### Workplaces of size 1 (B1) — SETTLED: they are kept

Both of us wrote this removal independently; it was live on your branch and commented out
on ours. **The decision is to keep one-person workplaces in the population**, so the
removal has been deleted from `social_contacts_usa2026.R` rather than left commented out.

This changes no results. The block was inert on this side — every line of it was a
comment — so deleting it is behaviour-identical, and
`contact_matrix_usa_tx_gaines_c1000.xml` does not need regenerating.

**It does mean runs made on `measles_usa_rm` are not comparable on this point**, because
there the removal was live: a one-person workplace contributes no workplace contacts but
still counts as employment in `age_distr_workplace`, the denominator of the conditional
contact rate, so removing them raised that rate for everyone else and shifted the
workplace size distribution feeding the `uniroot` adjustment factor. Your
`rStride_r0_measles.R` carried `## With single-person workplaces` and
`## Without single-person workplaces` variants; the baseline is now the *with* case.

## 5. One thing to know about `PopSnapshotWriter`

It is in, unchanged. But `SimController.cpp:149` calls it **unconditionally** — there is no
config gate — so it writes `households.csv`, the susceptibles-by-age file and the
population snapshot next to *every* run. That is the remaining ~16 % of the gtester's
92 s → 107 s, and at calibration scale it means three files per experiment across
thousands of runs.

B7 asked whether it should be gated by a config flag and the answer was left open. Worth
settling before the next large grid.

## 6. What still needs you

1. **A1 — runs to repeat.** The School presence fix (`c79a76f`) was missing from
   `measles_usa_rm`, so a person who stayed home from school was never marked present
   again. Runs made on that branch since the merge base are suspect. Which ones need
   re-running?
2. **B5 — the clustering and hesitancy capability.** It merged cleanly and is in the
   baseline, which is exactly the risk the item named: nobody on this side has reviewed
   the modelling. A walkthrough would be good, and `immunity_clustering_plan.md` Phase I
   is waiting on it.
3. **B1** — above.
4. **The `immunity_6region` data.** If those populations and immunity profiles are meant
   to be shared, they need to reach the repository; otherwise the Spartanburg and CT
   configurations cannot be run by anyone else.

## 7. Getting the branch

```sh
git fetch origin
git checkout integration/measles-usa
```

Note: the build unpacks population archives at configure time, so if you have an existing
build directory, reconfigure — `cmake` is now set to re-check for new archives, but an
old build directory predates that.

---

## 8. Regression suite result — the merge changes nothing

Added 2026-10-05, after running the R regression suite on both sides.

**Headline: your merged code does not change any result.** The full suite (23 scenarios
x 5 seeds = 115 runs) was run twice — once on `measles_usa` before the merge, once on
`integration/measles-usa` after it — and the outputs are identical:

- same 115 rows, same 74 columns
- `num_cases` identical in every run
- every numeric output column identical

The only column that differs is `holidays_file`, and only because it embeds each run's own
timestamped output directory (`..._124032_...` vs `..._123317_...`). Same filenames.

**The diffs against the stored references are pre-existing, not caused by the merge.**
Both runs were compared against `tests/*.rds` and returned *identical* lists of changed
scenarios:

| stream | scenarios differing from reference | before merge | after merge |
|---|---|---|---|
| summary | 23 | ✓ | identical set |
| incidence | 22 | ✓ | identical set |
| prevalence | 22 | ✓ | identical set |
| contacts | `covid_logParticipants` | ✓ | identical set |
| participants | `covid_daily` | ✓ | identical set |

The references are stale relative to `measles_usa`, which has had 53 commits since the last
reference reset, and the six files were already split across two dates six weeks apart
(Jul 6 and Aug 25 — plan F6). Magnitudes are small and consistent with drift rather than a
break: `covid_collectivity_mixing` +3.7 %, `covid_collectivity` +2.3 %, three scenarios
−1.5 %, most under 1 %, and `covid_airborne`, `covid_hosp` and `covid_suscept_adapt`
unchanged at 0.0 %.

A single clean `rrv()` on the consolidated baseline is what resolves this, and the result
above is what makes it safe to do: the reset will describe one known commit rather than a
mixture of code states.

### A pre-existing failure, unrelated to the merge

`rStride_gtester_covid19.R` does not reach its comparison section. It halts in the ABC
test, identically on both builds:

```
Error in abc_out[i_ref] <- data_incidence_all[sim_date == sum_stat_obs$date[i_ref],  :
  replacement has length zero
Calls: run_rStride_abc
```

`rStride_main_abc.R:221`. `sum_stat_obs.rds` is absent from the project directory, so the
fallback at line 180 derives the reference period from the simulation's own dates, yet a
row of `sum_stat_obs$date` then matches no `sim_date`. Under investigation separately.

Consequences worth knowing: the comparison above had to be driven directly against the
completed runs rather than through the script, `out_abc` is the one reference stream that
could not be checked, and the suite cannot gate CI until this is fixed (plan section 9.4,
point 3).
