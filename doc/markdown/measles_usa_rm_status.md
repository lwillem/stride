# Status of the refactoring, for the `measles_usa_rm` discussion

**Date:** 2026-10-15 (state of play in `rStride_refactoring_plan.md` §0 is dated 2026-10-07)
**For:** the `measles_usa_rm` branch owner
**Baseline:** `master` at `7f43e0e`; `measles_usa_rm` was fast-forwarded to `master` on 2026-10-05, tag `pre-refactor-2026-10`.

Detail on how each of their changes was treated is in `measles_usa_rm_merge_result.md`;
the decisions behind it are in `measles_usa_rm_discussion.md`.

---

## 1. The consolidation (2026-10-05)

- `origin/measles_usa_rm` (36 commits, tip `2b927f4`) was merged with `measles_usa`. The
  result is `master`. Their commits were not rewritten.
- **Taken unchanged:** `PopSnapshotWriter`, `MeaslesClustering.R`,
  `HouseholdClusterFactory_USA.R`, `ImmunityProfileFactory_USA.R`,
  `immunity_measles_dummy.xml`, the `RnMan::Shuffle` copy-constructor fix.
- **Taken with adaptation:** the `max_age` parameter of `get_cnt_data()`; the Connecticut
  FIPS fix (ported by hand, the population factory was renamed); the age-specific
  hospital probabilities in `rStride_measles_explore.R`.
- **Three changes to their workflow:** `ImmunitySeeder.cpp` was reverted (see §2); the
  population snapshot is now opt-in; workplaces of size 1 are kept rather than removed.
- **Results at risk:** the School-presence fix (`c79a76f`) was missing on their branch. A
  person who stayed home from school was never marked present at school again. Any run
  from their branch since the merge base (2026-07-23) is suspect until re-run.

## 2. Immunity seeding: their pruning was replaced by a fast path

- Their `ImmunitySeeder` rework (`a8f5095`) made the `influenza_c` test (immunity 99.91 %)
  run for more than 3 minutes instead of 6 s: `isExhausted()` is a linear household scan
  per rejected draw.
- **Phase 5c** (merged 2026-10-06, `f57bfd4`) replaced the sampler. The test now takes
  2.1 s (was 16 min); seeding is 25-150× faster (0.16 s or less on Gaines and Dane).
- **The old and new sampler are not equivalent in distribution.** Per-age marginals are
  identical, but the old sampler over-sampled people in small households (Dane, 50 %
  coverage: 77 % immune if living alone, 23 % in households of 8+). The fast path does
  not. This was accepted as a correction. Runs seeded before `f57bfd4` carry the bias.
  User-facing note: `rStride_architecture.md` §4.2.

## 3. Contact probability: `min` today, `mean` intended (F8, F10, Phase 2b)

**What it is.** For a pair of people, each age's reference number of contacts gives an
individual contact probability (`Infector.cpp`, around line 265). The two are combined
into one probability for the pair:

- **`min`** (committed baseline, default): the smaller of the two.
- **`mean`**: the average of the two. This raises realised contacts per person, so the same
  transmission probability gives a higher R0.

**What happened.**
- An uncommitted `mean` edit was found in a working tree (2026-09-26) and discarded.
  `min` is the baseline, so the regression references describe the committed code.
- Why it matters (F10): all ten disease files carry an R0 fit made under `min`.
  `TransmissionProfile` inverts it to get the transmission probability from a requested
  R0. Under `mean` every fit is silently wrong: `r0 = 12` would give a different
  effective R0 with no warning.
- **Option C was taken** (`f314514`, 2026-09-28): the rule is a configurable setting
  (`ContactProbabilityRule { Min, Mean }`, key `contact_probability_rule`, validated in
  `Misc.R`). The default is `min`, so existing calibrations and references stay valid.
- `mean` is **deferred, not cancelled**. Phase 2b is: calibration infrastructure with
  provenance, refit all ten files under `mean`, compare both rules on identical
  populations, then flip the default in one dedicated PR with a reference reset.
- Related: the cluster-size adjustment factors in `social_contacts_usa2026.R` were
  calibrated with averaging in mind (averaging raises realised contacts vs the minimum).

**Why it matters for the branch.**
- **Their measles scripts already opt in to `Mean`:** `rStride_measles_explore.R:72` and
  `rStride_r0_measles.R:65` set `contact_probability_rule <- "Mean"`.
- The committed `disease_measles_usa.xml` fit is, per the plan, from a `min` run. So
  measles runs under `Mean` with that fit give a higher effective R0 than requested.
- The fit's provenance is already broken: its `<label>` points to a population file in a
  superseded install directory that no longer exists (F10).
- **To settle:** is `disease_measles_usa.xml` the `min` fit, or has it been refitted
  under `Mean`? If it is the `min` fit, their `Mean` runs need a refit (Phase 2b) or must
  go back to `min`. *(My reading of the plan and the scripts; to be confirmed with the branch owner.)*

## 4. Other changes since the merge

- **Infrastructure:** OpenMP works (it was silently off); kernel failures are no longer
  silent; stable install root `~/opt/stride`; CI builds on gcc and clang; nightly
  regression run (informational).
- **Venue memory (PR #11):** `sizeof(Person)` 1192 → 224 bytes; identical results in
  single-threaded runs.
- **Population file by header name (PR #18):** the kernel reads `community_weekend` /
  `community_weekday`; `primary_community` / `secondary_community` still work as
  synonyms for now. The USA population generator writes the new names; the separate
  `stride_population` repo still has to follow.

## 5. Planned work that touches their area

- **Phase 5d** (decided, not started): household clusters stop being a separate
  mechanism and become generic leisure venues. Membership, meeting days and symptomatic
  withdrawal move into the population data and per-type settings. Their
  `HouseholdClusterFactory_USA.R` and the `household_clustering` calendar category would
  be reworked or replaced.
- **F12.7:** the venue generator makes 17 pools spanning two days: a modelling defect
  and a data race. Fixing it changes results and needs a reference reset.

## 6. Questions for the branch owner

1. Do they need household clusters as they work today (per-pair probability, community
   reduction by cluster size), or can a venue pool with a membership list do the job?
2. Which contact rule do their measles runs use, and which rule produced
   `disease_measles_usa.xml`? Do they need a `mean` refit?
3. Have they re-run anything seeded with the old immunity sampler, or run with the missing
   School-presence fix?
