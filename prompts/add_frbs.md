# Add CHIME Repeaters and Hosts to the FRB Repo

## Goals

Bring the FRB repo up to date with the localized CHIME repeaters and their host galaxies, so that `zdm/papers/Mo_Repeaters/CHIME_FRB_hosts.csv` (built by `zdm/papers/Mo_Repeaters/py/build_host_table.py`) can be regenerated from the repo alone, without the literature values hard-coded in that script.

All decisions behind this doc are recorded in `zdm/papers/Mo_Repeaters/claude_prompts/frb_host_table.md` (Q&A Q1–Q40, and the "Missing Repeaters" table). The literature values there were verified against the paper PDFs in `zdm/papers/Mo_Repeaters/papers/`.

## Context

- `frb/data/FRBs/FRBs_base.csv`: one row per FRB, with `repeater`, `telescope`, `z`, `refs` and `P(O|x)`.
- `frb/data/FRBs/FRB*.json`: per-FRB JSON files.
- `frb/data/Galaxies/public_hosts.csv`: the input list of hosts, read by `frb/builds/build_hosts.py` (`main(frbs, ...)` → `run(host_row)`). `run()` queries the public surveys (DECaL, Pan-STARRS, WISE, ...), merges any literature tables, corrects for Galactic extinction (`photom.correct_photom_table`, with EBV from `nebular.get_ebv`), and writes `frb/data/Galaxies/<FRB>/FRB<FRB>_host.json`.
- `frb/data/Galaxies/Literature/all_refs.csv`: an ordered list of literature tables (`Table, Format, Reference, DOI`). Order matters: the last measurement of a quantity wins. Each table has `Name` (`HG<FRB>`), `ra`, `dec` and quantity columns, with `_err` or `_loerr`/`_uperr` (see `gordon2023_derived.csv`).
- `frb/galaxies/defs.py`: the valid keys. Upper limits use `_err = 999`; "no measurement" uses `-999`.
- Papers (PDFs): `zdm/papers/Mo_Repeaters/papers/`.

Reference labels to use:
- `Michilli2023` (arXiv:2212.11941)
- `Ibik2024a` (arXiv:2304.02638)
- `Ibik2024b` (arXiv:2409.11533)
- `Moroianu2025` (arXiv:2509.05174)
- `Ravi2023` (arXiv:2211.09049)
- `Hewitt2024` (arXiv:2410.17044)
- `Bhardwaj2025` (arXiv:2506.11915)
- `Eftekhari2024` (arXiv:2410.23336)
- `Shah2024` (arXiv:2410.23374)
- `Leung2025` (KKO catalog, ApJS 280, 6; arXiv:2502.11217)
- `Leung2025b` (ApJL 991, L25; arXiv:2507.16816)

## Prompts

1. **FRB20181030A host.**
   (a) Run `python zdm/papers/Mo_Repeaters/py/fix_frb20181030A_host.py` in the `astro` env. Then `git rm frb/data/Galaxies/FRB20181030A_host.json` and `git add frb/data/Galaxies/20181030A/` (the user does the git steps).
   (b) Check that `FRBHost.by_frb(FRB.by_name('FRB20181030A'))` loads it.
   (c) `Literature/bhardwaj2021_derived_FRB20181030A.csv` actually holds the FRB20201124A values (HG20201124A). Rename it or remove it, and check `all_refs.csv`.
   Log your work below.

2. **Lower-limit convention.** In `frb/galaxies/defs.py` (next to line 30), document `_err = -998` as a **lower limit** for derived quantities (agreed in Q37). Search the repo for anything that treats negative errors specially (e.g. `chk_fill` in `build_hosts.py`, and the table builders in `frb/galaxies/utils.py`), and make sure -998 survives a build. Log your work below.

3. **Rename the two repeaters listed under later bursts (Q32).** FRB20231204A is the same source as **FRB20190303A**, and FRB20231128A is the same source as **FRB20191106C** (KKO §6.1, §6.6). For each:
   - rename the `FRBs_base.csv` row, keeping the localization of the later burst (it is the better one) and noting the burst name in `refs`;
   - rename `Galaxies/<old>/` to `Galaxies/<new>/` and the `FRB` field in the host JSON;
   - update `public_hosts.csv`;
   - check whether an `FRB*.json` exists under either name.
   Ask the user before deleting anything. Log your work below.

4. **Base-table fixes.**
   - FRB20231201A: `z` = 0.119 should be 0.1119 (host JSON and KKO Table 2; Leung+2025b Table 1 repeats the typo).
   - The 12 CHIME rows with a `z` but P(O|x) < 0.9 have no published redshift in the KKO paper: 20230410A, 20230616A, 20230702A, 20230828A, 20230918A, 20230923A, 20230924A, 20231006B, 20231102A, 20231223D, 20231224A, 20240210C. Ask the user where these z values came from, and whether to blank them.
   - FRB20181119A: the `FRB20181119A.json` z = 0.26064 has no literature support (no host is proposed; Ibik+2024b App. A.4). Remove it.
   - Consider cross-matching the CHIME rows against `CHIME_catalog-2021-1-27.json` (`repeater_of`) to flag any other repeaters.
   Log your work below.

5. **Add the missing CHIME repeaters to `FRBs_base.csv`.** Use the localizations in the papers, set `repeater=TRUE`, and record the papers in `refs`. Read coordinates, DM, z and the localization ellipse from the PDFs, not from memory.

   | FRB | Loc. telescope | z | Source |
   |---|---|---|---|
   | FRB20180814A | CHIME | 0.06835 | Michilli2023 |
   | FRB20190110C | CHIME | 0.12244 | Ibik2024a |
   | FRB20200223B | CHIME | 0.06024 | Ibik2024a |
   | FRB20190417A | EVN | 0.12817 | Moroianu2025; Ibik2024b |
   | FRB20220912A | DSA (+EVN) | 0.0771 | Ravi2023; Hewitt2023 (arXiv:2312.14490) |
   | FRB20240114A | MeerKAT / EVN | 0.130287 | Bhardwaj2025; Tian2024 (arXiv:2408.10988) |
   | FRB20240209A | CHIME (+KKO) | 0.1384 | Shah2024; Eftekhari2024 |
   | FRB20190208A | EVN | — | Hewitt2024 |
   | FRB20181119A | CHIME | — | Michilli2023 |

   Also create or update the `FRB*.json` files with `frb/builds/build_frbs.py`. Log your work below.

6. **Build host JSONs for the new repeaters.** Add rows to `public_hosts.csv` (FRB, host `Coord`, `P_Ox`, `z`, `Projects=CHIME`, `References`). Host coordinates:
   - 20180814A: PanSTARRS-DR1 J042256.01+733940.7
   - 20190110C: 16h37m16.43s +41d26m36.30s (Ibik2024a Table 3)
   - 20200223B: 00h33m04.68s +28d49m52.60s (Ibik2024a Table 3)
   - 20190417A: take the host centroid from Moroianu2025. The FRB is at 19h39m05.8919s +59d19m36.99s.
   - 20220912A: PSO J347.2702+48.7066
   - 20240114A: DESI J212739.84+041945.8 (Bhardwaj2025 Table 3)
   - 20240209A: 289.85036 +86.06090 (Eftekhari2024 Table 2)
   - 20190208A: host O4 at the EVN position 18h54m11.27s +46d55m21.67s (Hewitt2024); r ≈ 27, so there is no survey photometry.

   Then run `build_hosts.main([...])`. Check each JSON's survey photometry against the papers. The papers quote these uncorrected catalog values, so they should agree before the extinction correction:
   - 20190110C: DECaL r = 18.009
   - 20200223B: DECaL r = 16.080
   - 20220912A: PS1 r = 19.65
   Log your work below.

7. **Literature derived-quantity tables.** Create one table per paper in `Galaxies/Literature/` (e.g. `ibik2024a_derived.csv`, `ibik2024a_nebular.csv`), add each to `all_refs.csv`, and rebuild the hosts. Values (verified; see `frb_host_table.md`):

   | Host | Mstar [Msun] (method) | SFR [Msun/yr] |
   |---|---|---|
   | 20180814A | 10^10.78 +0.12/−0.18 dex (Prospector; Michilli2023 Tab 3) | `SFR_SED` < 10^-0.5 (upper limit: `_err=999`) |
   | 20190110C | 2.5 +0.10/−0.17 e10 (Prospector; Ibik2024a Tab 3) | `SFR_nebular` 0.1575(6); `SFR_SED` 0.54 ± 0.04 |
   | 20200223B | 5.6 +1.14/−0.93 e10 (Prospector) | `SFR_SED` 0.59 ± 0.04 (no Hα: AGN) |
   | 20191106C | — (see prompt 8) | `SFR_nebular` 1.53 (SDSS fiber; Ibik2024a) |
   | 20190303A | 10^10.75(3) (SDSS-collaboration value; Michilli2023 Tab 3) | log SFR 0.84(4) (SDSS-collaboration value). There is no `defs` key for an SDSS-pipeline SFR; ask the user |
   | 20190417A | 10^7.88 +0.12/−0.14 (Prospector; Moroianu2025) | `SFR_nebular` 0.19 ± 0.01 |
   | 20220912A | 10^10.0 ± 0.1 (Prospector; Ravi2023 Tab 2) | `SFR_nebular` 0.1 lower limit (`_err=-998`) |
   | 20240114A | 10^8.55 +0.12/−0.14 (Prospector, Kroupa; Bhardwaj2025 Tab 3) | `SFR_nebular` 0.061 +0.004/−0.003 |
   | 20240209A | 10^11.34 ± 0.01 (Prospector, Kroupa; Eftekhari2024 Tab 2) | `SFR_SED` < 0.36 (upper limit) |

   Notes:
   - `Mstar` is documented in `defs.py` as "CIGALE (or Prospector if Gordon2023)". Update that comment to say the method is set by `Mstar_ref`, and list the Prospector refs.
   - Leung2025b (Table 1) CIGALE / NED-LVS masses for the KKO hosts lacking `Mstar`: 20230926A 10.49, 20231011A 9.59, 20231123A 9.42, 20191106C 9.47, 20231201A 9.47, 20231229A 9.87, 20231230A 10.04 (CIGALE); 20230222B 10.19 and 20231223C 10.40 (NED-LVS). All are Chabrier, with no per-object errors. These need a way to record the method (CIGALE vs NED-LVS). Propose one to the user.
   Log your work below.

8. **Redo stellar masses with Prospector (Q26, Q38).**
   - FRB20191106C: Chang+2015 (via Ibik2024a) gives 4.5e10, but Leung2025b CIGALE gives 10^9.47, a 1.2 dex conflict. Fit it with Prospector (Gordon2023 setup).
   - Do the same for the other KKO hosts that have only CIGALE or NED-LVS masses, and for FRB20181030A, whose Prospector fit (Bhardwaj2021b) used a different setup.
   - Where a host already has a repo value, the repo value stays until the new fit is reviewed (Q26).
   Log your work below.

9. **Regenerate the CHIME table.** Once prompts 1–7 are done:
   - remove the corresponding entries from `LIT_HOSTS`, `LIT_SUPPLEMENT`, `RENAME` and `FIX_181030A` in `zdm/papers/Mo_Repeaters/py/build_host_table.py`;
   - rerun it;
   - diff the new `CHIME_FRB_hosts.csv` against the previous one. Every change should be explainable (e.g. the magnitude source).
   Log your work below.

## Q&A

## Logs

### Prompt 1: FRB20181030A host (2026-10-04, Claude Opus 5.5)

**(a)** Ran `python py/fix_frb20181030A_host.py` in the `astro` env, from `zdm/papers/Mo_Repeaters/`. It wrote `frb/data/Galaxies/20181030A/FRB20181030A_host.json`.
- Photometry renamed to `defs.valid_filters` (DECaL_grz, 2MASS_jhk, WISE_W1–W4, GALEX_FUV/NUV), each with `_ref = Bhardwaj2021b`. The Herschel entries were dropped.
- The old file had `DESI_g_ref` twice and no `DESI_r_ref`. This does not matter, because the script sets every `_ref` itself.
- `EBV` changed from 0.42 (the host Balmer-decrement value) to the Galactic 0.0955 (`nebular.get_ebv`).
- Derived values: `Mstar` 5.8e9 (+1.6e9/−2.0e9), `AV_nebular` 1.3 ± 0.2, `SFR_photom` 0.36 ± 0.08, and `SFR_nebular` 0.033 with `_err = -998` (lower limit). All have `_ref = Bhardwaj2021b`. The script dropped `Z_spec` and `Mstar_spec`.
- The redshift block is unchanged (z = 0.0039).
- Left for the user: `git rm frb/data/Galaxies/FRB20181030A_host.json` and `git add frb/data/Galaxies/20181030A/`. The old file still has the uncommitted trailing-comma fix, so use `git rm -f`.

**(b)** `FRBHost.by_frb(FRB.by_name('FRB20181030A'))` loads the new file. `by_frb` reads `Galaxies/<name>/`, so the old top-level JSON is never used. The result is `<FRBHost: 10:34:23.10 +73:45:53.92, FRB=FRB20181030A z=0.0039>`, and `derived`, `photom` and `redshift` are as above.
- Prompt 2 still needs to check that -998 survives any later `build_hosts` run or table builder.

**(c)** `Literature/bhardwaj2021_derived_FRB20181030A.csv` (commit bdec78a8, 2023) is byte-identical to `fong2021_derived.csv` (`cmp`), which holds the HG20201124A row. Nothing referenced it: it is not in `all_refs.csv`, and no code uses it. I **removed** it (plain `rm`); the user should `git rm` it.
- The real FRB20181030A values (Bhardwaj+2021b, ApJL 919 L24) are now in the host JSON itself (a). The Q26-style alternative would be a proper `bhardwaj2021b_derived.csv` with `Name = HG20181030A`, added to `all_refs.csv`. This would only matter if 20181030A is ever rebuilt with `build_hosts`, and it is not yet in `public_hosts.csv`.
- `all_refs.csv` is unchanged. Its only Bhardwaj entry is `bhardwaj2021_derived.csv` (Bhardwaj2021, DOI 10.3847/2041-8213/abeaa6, the FRB20200120E/M81 paper).
- Side note: `bhardwaj2021_derived.csv` (ra 148.888, dec 69.065 = FRB20200120E) has **no `Name` column**, unlike the other literature tables. It is probably not matched to any host by the literature merge. This was not changed.

### Prompt 2: Lower-limit convention (2026-10-04, Claude Opus 5.5)

**defs.py.** Documented `_err = -998.0` as a lower limit (derived quantities only) next to the 999/-999 lines, at about line 31. I also added a short block at the head of `# Derived quantities` that lists all three codes (999 upper limit, -998 lower limit, -999 no error/measurement).

**Code that treats special or negative errors.** I searched for 999, 998, `err < 0`, `chk_fill` and `fill_value`:
- `build_hosts.chk_fill` uses `np.isclose(value, -999.)`. -998 is not "close" (the tolerance is about 0.01), so -998 errors are not treated as fill values. They get a `_ref` and are merged. This function is only used for the photometry tables.
- The `build_hosts.run` derived-literature loop copies the value, error and `_ref` verbatim. It has no sign checks, so -998 passes through.
  - Caveat: this loop also has no `-999` skip. A -999 in a derived table is copied as a value.
- `build_hosts.run` nebular-literature loop: verbatim copy.
  - `nebular.calc_lum` (L178) sets `Lum_err = -999` for any line-flux error ≤ 0. This applies only to `neb_lines` errors, not to derived quantities.
- `FRBGalaxy.vet_one` checks keys only, not values.
- `galaxies/utils.build_table_of_hosts` copies dict values verbatim.
- `photom.correct_photom_table`, `parse_photom` and `surveys/catalog_utils` (L326, L384) only look at photometric magnitudes and errors.
- `figures/utils.log_me` (used by `figures/galaxies.py` for Mstar and SFR) treats `err < 0` as "no error bar". A -998 point plots without an arrow and does not crash.
  - Pre-existing issue, not changed: `err = 999` would give `log10(val - 999)` = NaN error bars.
- Historical note: `build_fg.py:88` and `old_build_hosts.py` used `SFR_nebular_err = -999` and `[NII]_err = -999` ("Upper limit") inconsistently. Existing JSONs with derived `_err = -999` therefore mean "no error", not a limit.

**Bug fixed: `build_hosts.read_lit_table`.** It had `idx = int(np.where(match)[0])`, which raises `TypeError` under NumPy 2.5 (the `astro` env). Every `run()` that matched a literature table crashed. The line is now `int(np.where(match)[0][0])`. This is the only instance of the pattern in `frb/`. Prompts 6 and 7 need this fix.

**Round-trip test.** I wrote a temporary `Literature/zz_tmp_lowerlim_derived.csv` (since deleted) with an HG20181030A row: `SFR_nebular_err = -998`, `SFR_SED_err = 999`, and `Mstar` with `_loerr`/`_uperr`. A temporary `all_refs` pointed at it. I then ran `build_hosts.run(row, lit_refs=..., skip_surveys=True, out_path=<scratch>)`.
- The output JSON has `SFR_nebular_err = -998.0`, `SFR_SED_err = 999.0`, and the `Mstar` lo/up errors intact.
- `FRBHost.from_json` keeps both values.
- `utils.build_table_of_hosts()` (repo JSONs) gives HG20181030A `SFR_nebular_err = -998.0`.
- Conclusion: -998 survives a build.

**Tests.** `test_frbhosts`, `test_galaxies` and `test_photom` pass. `test_build::test_host_build` fails on survey access before the literature step: first on PS1 metadata, then on "NEDLVS environment variable not set". This is environmental.

**Warnings for later prompts:**
- The `astro` env has a **stale non-editable install** of `frb` (site-packages, dated Sep 11, with no `Galaxies/20181030A/`). A script run outside the repo root picks it up, for example `python py/build_host_table.py` from `zdm/` (prompt 9), or `fix_frb20181030A_host.py` in prompt 1, which only used it for `get_ebv`. Run `pip install -e .` in the `astro` env, or set `PYTHONPATH`.
- `build_hosts` survey queries need `$NEDLVS` set (prompt 6).
- Correction to the Prompt 1 side note: `read_lit_table` matches literature rows by **ra/dec (< 1″)**, not by `Name`. So `bhardwaj2021_derived.csv`, which has no `Name` column, still matches its host.

### Prompt 3: Rename the two repeaters to their first-burst names (2026-10-04, Claude Opus 5.5)

| Old (later burst) | New (source) | Check |
|---|---|---|
| FRB20231204A | FRB20190303A | CHIME Cat 1 (`chimefrbcat.csv`, `CHIME_catalog-2021-1-27.json`) lists FRB20190303A at RA 208.03, Dec 48.24, DM 221.67. This is consistent with the 2023 burst (207.999, 48.116, DM 222). |
| FRB20231128A | FRB20191106C | Not in Cat 1. Identification per KKO §6.6. |

**`FRBs_base.csv`.** Only the `Name` and `refs` fields changed. The localization, DM, fluence, z and P(O|x) of the later (KKO) burst were kept. `refs` changed from `Leung+2025` to `"Leung+2025,burst:FRB20231204A"` (and `burst:FRB20231128A`). `build_frbs` splits `refs` on commas, so the burst name becomes its own `refs` entry (`burst:` prefix) in the FRB JSON.
- The rows stay in their old (2023) positions in the file, to keep the diff small.
- The file has **CRLF** line endings. My first edit with Python text mode converted the whole file to LF. I redid it in binary mode from the HEAD version, so the diff is two lines. The same applies to `public_hosts.csv`.

**Host directories.** `git mv Galaxies/20231204A → Galaxies/20190303A` and `Galaxies/20231128A → Galaxies/20191106C`, with the JSONs renamed to `FRB20190303A_host.json` and `FRB20191106C_host.json`. The `"FRB"` field was updated in each. Nothing else in the host JSONs referred to the old names.

**`public_hosts.csv`.** The two rows were renamed (20231204A → 20190303A, 20231128A → 20191106C). Coordinates and references are unchanged.

**FRB JSONs.** Both `FRB20231204A.json` and `FRB20231128A.json` existed; nothing existed under the new names.
- Built `FRB20190303A.json` and `FRB20191106C.json` with `build_frbs.main([...])`. They are identical to the old files apart from `FRB`, `refs` and the last digit of DMISM.
- `frb.list_of_frbs()` and `galaxies.utils.list_of_hosts()` glob `FRB*.json`, so the old files would have listed each source twice. With the user's OK, I removed them with `git rm` (git shows them as renames).
- Checks:
  - `FRB.by_name` and `FRBHost.by_frb` load both new names.
  - `list_of_hosts()` returns 95 hosts with no duplicates, including FRB20190303A and FRB20191106C.

**Other references.** The only other repo files that mention the old or new names are `chimefrbcat.csv` and `CHIME_catalog-2021-1-27.json`. These are catalog data and were left alone. `FRB_pulses.csv` and the literature tables have no entries for either name.
- Downstream: `zdm/.../build_host_table.py` `RENAME` handles this mapping. It can be dropped in prompt 9.

### Prompt 4: Base-table fixes (2026-10-04, Claude Opus 5.5)

All `FRBs_base.csv` edits were made in binary mode and keep the CRLF line endings. Each change touches only the `z` field.

**FRB20231201A.** `z` changed from 0.119 to **0.1119** in `FRBs_base.csv`. I rebuilt `FRB20231201A.json` with `build_frbs`. The host JSON and `public_hosts.csv` already had 0.1119.
- The row still has `repeater = FALSE`. KKO/Leung+2025b do not mark it as a repeater, so I left it.

**The 12 insecure z (P(O|x) < 0.9).** The values arrived in commit 0c665727 ("kko", profxj, 2025-02-27), with `refs = Leung+2025` and no source file or note. KKO Table 3 gives no z for any of them. Their origin is unknown; the user did not say.
- Per the user's choice, I **blanked** z in `FRBs_base.csv` for 20230410A, 20230616A, 20230702A, 20230828A, 20230918A, 20230923A, 20230924A, 20231006B, 20231102A, 20231223D, 20231224A and 20240210C.
- I rebuilt their `FRB*.json` with `build_frbs.main([...])`. The JSON diffs drop `"z"`; the only other change is DMISM in the 15th decimal place (NE2001 rounding).
- `P(O|x)` was kept, so the candidate associations are still recorded.
- None of the 12 has a host JSON or a `public_hosts.csv` row.

**FRB20181119A.** I removed `"z": 0.26064` from `FRB20181119A.json` by editing the text, since the source is not yet in `FRBs_base.csv`. `FRB.by_name('FRB20181119A').z` is now `None`.
- The JSON still has `refs = ["Astroflash"]` and a placeholder ellipse (a = b = 0.01). Prompt 5 replaces both when the row is added and the JSON is rebuilt.

**Cross-match with CHIME Cat 1 (`CHIME_catalog-2021-1-27.json`, `repeater_of`).** Cat 1 has 18 repeating sources. I matched all 174 base rows with coordinates against the mean burst position of each source (< 1°, |ΔDM|/DM < 5%), and checked names as well. The only matches are rows already flagged `repeater = TRUE`:
- FRB20180916B
- FRB20190303A (4.4′ / 7.5′)
- FRB20181030A, which Cat 1 lists under the source name **FRB20181030B** (Cat 1's burst FRB20181030A has `repeater_of = FRB20181030B`)
- FRB20121102A, which Cat 1 lists as source **FRB20181119D**

**No new repeaters to flag.** The Cat 1 repeating sources missing from the base table are:
- 20180814A, 20181119A, 20190208A and 20190417A, which prompt 5 adds;
- 20180908B, 20181017A, 20181128A, 20190116B, 20190117A, 20190209A, 20190212A, 20190213A, 20190222A and 20190604A. These have only Cat 1 (arcmin) positions and no host, so they are out of scope.

**Other issues noticed (not changed):**
- **Name mismatch for the ASKAP FRB.** The `FRBs_base.csv` row for the ASKAP FRB (RA 326.105) is named **FRB20180924A**. Everything else in the repo uses **FRB20180924B**: the `FRB20180924B.json`, `Galaxies/20180924B/`, `public_hosts`, the literature/PATH tables and `test_build`. FRB20180924A in TNS and Cat 1 is a CHIME one-off at RA 35.46. `build_frbs.main(['all'])` would therefore write a stray `FRB20180924A.json`.
- **Missing coordinates.** Five DSA rows have no ra/dec: FRB20220121B, FRB20220424E, FRB20220801A, FRB20220926A and FRB20221002A (Sherman23). `build_frbs` would fail on them.
