# Immunity seeding: age-split and household clustering

**Date:** 2026-09-26
**Context:** measles/USA work with Regina Manansala
**Code inspected:** `measles_usa` @ `97f9b39` and `origin/measles_usa_rm` @ `a8f5095`
**Related:** `measles_usa_rm_discussion.md` (B5, B8), `rStride_architecture_and_refactoring.md`

## Aim

1. **Age ≥ 40 — unclustered immunity.** Treated as acquired when measles still
   circulated; no reason for it to correlate within households.
2. **Age < 40 — clustered immunity.** Vaccine-derived, so correlated within households
   because vaccination sentiment and behaviour are household attributes.
3. **Compare clustering mechanisms.** First candidate: a person is more likely to be
   vaccinated when another household member is, and conversely more likely to be
   unvaccinated when another is unvaccinated.
4. **Make high uptake fast.** At a 90 % target the current sampler spends most of its
   draws on people who are already immune.

---

## 1. What the code does today

### 1.1 Call path

`SimBuilder` → `ImmunitySeeder::Seed()` runs **two independent passes**, both over
**Household** pools:

| Pass | Profile key | Rate source | Clustering knob |
|---|---|---|---|
| `"immunity"` | `run.immunity_profile` | `run.immunity_distribution_file` | `run.immunity_link_probability` |
| `"vaccine"` | `run.vaccine_profile` | `run.vaccine_distribution_file` | `run.vaccine_link_probability` |

(`vaccine_profile = "Teachers"` switches the second pass to School pools.)

Profiles: `None`, `Random`, `AgeDependent`, `Teachers`, `Cocoon` — validated in
`Misc.R:473`, for *both* passes.

### 1.2 The `AgeDependent` algorithm (`ImmunitySeeder::Random`)

```
count unvaccinated per age                      -> populationBrackets[age]
quota[age] = floor(count[age] * rate[age]);  numImmune = sum(quota)

while numImmune > 0:
    pick a household uniformly at random          (WITH replacement)
    shuffle all its member indices
    for each member in that random order:
        if not vaccinated and quota[age] > 0: vaccinate; quota[age]--; numImmune--
        draw U01; if U01 < (1 - link_probability): break   # leave, pick a new household
```

So **`*_link_probability` is the existing clustering knob**: the probability of
*continuing* within the same household after dealing with one member.

- `link = 0` → exactly one member is considered per household visit ⇒ effectively
  independent sampling, reached by rejection.
- `link → 1` → the walk continues through the whole household ⇒ strong within-household
  clustering.

### 1.3 Four properties worth naming

1. **Clustering is a by-product of the sampling walk, not a parameter.** `link` has no
   interpretation as a correlation, odds ratio or ICC. Its realised effect is confounded
   with household size and with how much quota is left when a household is drawn. You
   cannot currently state "the odds of vaccination given a vaccinated household member
   are X" — only "the walk continued with probability `link`".
2. **Both passes write the same thing.** Each sets `ConstantVaccine{"immunity",1,1,1}`
   and the eligibility test is `!IsVaccinated()`. After seeding, natural and
   vaccine-derived immunity are **indistinguishable** — per person and in the output.
3. **`AgeDependent` does not log.** `log_immunity` is `true` only on the
   `Random`/`Cocoon` path, so the `[VACC]` event lines are absent for exactly the profile
   this work uses.
4. **`floor()` per age systematically undershoots.** Summed over ~95 age classes the
   realised total falls below the target by up to ~95 people; at small county scale that
   is visible.

### 1.4 Why 90 % is slow — precisely

The loop is **rejection sampling with replacement**, and each rejection is expensive:

- Every iteration allocates a `vector<unsigned int>`, runs `iota`, and **shuffles the
  whole household** — then, at `link = 0`, considers exactly **one** person.
- The probability that a uniformly drawn household yields an acceptable person falls as
  quota fills. Filling the last few slots of one age class requires on the order of
  `N_households / (remaining candidates)` draws **each** — the coupon-collector tail.
- At 90 % uptake most drawn people fail `!IsVaccinated()`; at 95 % almost all do.
- It can fail to terminate: Regina's `a8f5095` fixes a **hang** when
  `link_probability = 0` and few household members remain to vaccinate.

Cost is therefore roughly *O(H · Σ_age 1/p_accept(age))* with a per-draw constant of
*O(household size)*, against an achievable *O(N)*.

---

## 2. What `measles_usa_rm` already solves

Regina's branch (`+167` in `ImmunitySeeder.cpp`) has already built much of this. **Adopt
it rather than rewriting it.**

| Addition | What it does |
|---|---|
| `RandomIndependent()` | Buckets unvaccinated people by age, shuffles, takes the first `quota`. **O(N), exact, no rejection.** Ignores household structure entirely. |
| Household **pruning** in `Random()` | Drops a household from the draw pool once no member is both unvaccinated and in an age class with open quota. Removes the wasted-draw pathology **without changing the sampling distribution**. |
| `run.mass_immunize` | Switch between the original workflow (default, unchanged) and "mass immunise adults → vaccinate children iteratively". |
| `run.vaccine_hesitancy_rate` | Excludes a **fixed count** of households (`round(rate·N)`, sampled without replacement) from the vaccine pass entirely. |
| `MeaslesClustering.R` (+327) | Analysis harness: `analyze_immunity_clustering()`, Cramér's V, household and community distributions, household-size effect. |
| `PopSnapshotWriter` (+256) | Writes population-wide immunity/susceptibility at *t = 0* — the measurement instrument for all of this. |

Her code comments record a validation result worth keeping: a variant that resolved each
household fully in one pass **over-clustered** (household pair-correlation excess +0.05
vs. the original's −0.02, consistent across 8 seeds), whereas pruning reproduced the
original (−0.02). That is the right acceptance test, and it already exists.

---

## 3. Gap analysis

| Requirement | Today | `measles_usa_rm` | Gap |
|---|---|---|---|
| Unclustered immunity for ≥ 40 | rejection sampling | `RandomIndependent` | **adopt** |
| Clustered immunity for < 40 | `link_probability` | same + pruning | mechanism is "continue walking", not "conditional on a vaccinated member" |
| Different mechanisms comparable | — | — | **new: a named, switchable mechanism** |
| Clustering as an interpretable parameter | no | no | **new: ICC / odds ratio, not `link`** |
| Fast at 90 % | no | yes | **adopt** |
| Distinguish natural vs vaccine immunity | no | no | decide whether it matters |

---

## 4. Proposed approach

### Phase I — adopt and verify (no new science)

1. Take `RandomIndependent()`, the pruning fix and `a8f5095` from `measles_usa_rm`
   (items B5/B8 of `measles_usa_rm_discussion.md`).
2. Confirm that with `mass_immunize = false` the regression references are **unchanged** —
   the default path must stay byte-identical.
3. Benchmark seeding time at 50 / 70 / 90 / 95 % on Gaines TX and Dane WI. This is the
   baseline the rest is measured against.
4. Adopt `PopSnapshotWriter` + `MeaslesClustering.R` as the standard instrument.

### Phase II — express the age split with **no C++ change**

The two existing passes already give exactly the split:

| Pass | Ages | File | `link_probability` |
|---|---|---|---|
| `immunity` (natural) | 40 … maxAge, 0 below 40 | `immunity_measles_usa_natural.xml` | **0** |
| `vaccine` (uptake) | 0 … 39, 0 above | `immunity_measles_usa_vaccine.xml` | **> 0** |

Both accept `AgeDependent` (`Misc.R:473` validates `vaccine_profile` against the same
list), and each reads its own distribution file and its own link probability. Zero the
complementary age range in each file and the passes do not interfere — the second pass
counts only `!IsVaccinated()`, so it operates on exactly the people the first left alone.

To do:
- Add `vaccine_distribution_file` and `vaccine_link_probability` to `config_default` in
  `rStride.R` (only `vaccine_profile` and `vaccine_rate` are there today, both defaulting
  to `None`/`0`).
- Extend `ImmunityProfileFactory_USA.R` to emit the two complementary files from one
  age-specific source, so the split is declared once.
- Enable logging on the `AgeDependent` path (§1.3.3) or rely on `PopSnapshotWriter`.

**This yields a working age-split model immediately**, before any sampler work, and it is
the baseline every later mechanism is compared against.

### Phase III — clustering as a named, switchable mechanism

Introduce `run.immunity_clustering` (and the `vaccine_` counterpart), defaulting to the
current behaviour:

| Value | Mechanism | Parameter |
|---|---|---|
| `none` | independent within age | — |
| `link` | **current** walk-continuation | `*_link_probability` |
| `household_propensity` | θ_h ~ Beta(α,β), mean = age target; each member vaccinated w.p. θ_h | ICC, or Beta dispersion |
| `hesitant_class` | households are *hesitant* (uptake p_lo) or *compliant* (p_hi) | hesitant fraction, p_lo |
| `conditional` | P(vaccinate ∣ k vaccinated in household) = logistic(β₀ + β₁·1{k>0}) | odds ratio β₁ |

Notes on the choice:

- `conditional` is the **literal statement of the stated hypothesis** — "more likely when
  another household member is". It is sequential, so the result depends on visiting
  order, and β₀ must be calibrated per age to recover the marginal target (`uniroot`, the
  same device already used for the contact adjustment factors).
- `household_propensity` and `hesitant_class` are the same family — a continuous vs. a
  two-point mixing distribution — and produce the same *kind* of positive within-household
  correlation **without** order dependence, which makes them exactly reproducible and
  trivially parallel. `hesitant_class` is a direct generalisation of Regina's
  `vaccine_hesitancy_rate` from all-or-nothing to two rates, and maps onto the CDC
  hesitancy data already in use.

**Recommendation:** implement `household_propensity` first (one interpretable parameter,
order-independent, exact-quota-able), keep `hesitant_class` as the
data-linked variant, and add `conditional` to test whether the sequential formulation
gives a materially different clustering shape. Keep `link` available so old runs remain
reproducible.

### Phase IV — one sampler that is both clustered and fast

All the mechanisms above reduce to the same computational problem: **sample an exact
per-age quota without replacement, with unequal, household-derived inclusion weights.**
That admits a direct algorithm with no rejection at all:

```
for each age stratum a:
    candidates  = unvaccinated persons of age a          # O(N) bucket pass
    quota       = largest_remainder(count * rate[a])     # not floor(), see 1.3.4
    for each candidate i:
        w_i  = weight from the household state/propensity
        key_i = U_i^(1/w_i)                              # Efraimidis-Spirakis
    take the `quota` candidates with the largest keys
```

- **O(N log N)** worst case, single pass, no redraws — the 90 % and 95 % cases cost the
  same as 10 %.
- Exact quota by construction; `largest_remainder` removes the `floor()` undershoot.
- `w_i ≡ 1` reproduces `RandomIndependent` exactly.
- For `conditional`, weights depend on the evolving household state, so keep the
  sequential walk but over a **pruned, shuffled candidate list** rather than
  with-replacement household draws.

### Phase V — measure, then calibrate

Report for every run, from `PopSnapshotWriter` + `MeaslesClustering.R`:

1. realised per-age marginal vs. target (must match within tolerance);
2. within-household **ICC** / pair-correlation excess over binomial;
3. distribution of the number of immune per household vs. the binomial expectation;
4. Cramér's V, and the household-size effect.

Acceptance: marginals match, and the clustering statistic moves **monotonically** with
the mechanism's parameter. Then calibrate that parameter to an external target (CDC
hesitancy, NIS-Child clustering) by `uniroot`, exactly as the contact adjustment factors
are fitted.

---

## 5. Questions to settle first

1. **Is 40 the intended boundary?** In 2026 a 40-year-old was born in 1986 — after
   licensure (1963) and inside the one-dose era, but before two-dose became routine
   (1989 ⇒ age ≈ 37). US convention treats those **born before 1957** (age ≈ 69+) as
   presumptively immune from circulation. So age 40 reads as *"pre-two-dose era"* rather
   than *"natural immunity"*. Worth confirming which is meant — it changes both the label
   and the rate profile. A graded transition may fit better than a sharp cut.
2. **Does the natural/vaccine distinction need to survive seeding?** Today both passes are
   indistinguishable afterwards (§1.3.2). If waning, dose count or boosting is ever
   wanted, the immunity source must be recorded per person.
3. **Household only, or wider?** Vaccine sentiment also clusters by school and
   neighbourhood. The pool machinery would support School or Community clustering with
   the same mechanism — is that in scope?
4. **Exact quota, or stochastic marginal?** A per-household coin flip gives a random
   realised total; exact quotas are reproducible but constrain the tail. The phase IV
   sampler can do either.
5. **How is clustering strength to be chosen** — fitted to data, or swept as a
   sensitivity parameter?

---

## 6. Suggested first step

Phase I + Phase II together: adopt Regina's sampler work, then express the ≥40 / <40
split using the two existing passes and two complementary distribution files. That gives
a running, fast, age-split model **without touching the C++ sampler**, and produces the
baseline against which every clustering mechanism in Phase III is then judged.
