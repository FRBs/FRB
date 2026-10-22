# Rebuild the Host JSONs Not Touched by GitHub #240

## Goals

Rebuild the **48 host JSONs** that GitHub PR #240 (`repeater_updates`) did not touch, with main's current code, so that all 103 hosts in `frb/data/Galaxies/public_hosts.csv` come from the same survey code, extinction law and flag conventions. The work happens on a new branch, `rebuild_hosts`, created from `main` after #240 is merged, and ends in a PR with a per-host summary.

Background and every decision behind this file: `prompts/add_frbs.md`, in particular prompts PR-5 to PR-13, the Q&A sections "PR review by Sunil" and "Finishing up" (FU-Q1 to FU-Q18), and the first attempt at this rebuild (Logs → Run 1 below; its JSONs were discarded).

## Context

- **Hosts:** all of `public_hosts.csv` (103) minus `PR_HOSTS` in `prompts/rebuild_pr_hosts.py` (55). That is 48 hosts; the "49" in add_frbs.md counted the stray `frb/data/Galaxies/20190614D/` directory, which has no `public_hosts.csv` row:
  ```
  20121102A,20180916B,20180924B,20181112A,20190102C,20190520B,20190523A,20190608B,20190611B,20190714A,20191001A,20191228A,20200120E,20200430A,20200906A,20201123A,20201124A,20210117A,20210320C,20210410D,20211127I,20211212A,20220105A,20220501C,20220610A,20220725A,20220918A,20221106A,20230526A,20230708A,20230718A,20230902A,20231226A,20240201A,20240210A,20240304A,20240310A,20230203A,20230222A,20230311A,20230703A,20230730A,20231005A,20231017A,20231025B,20231206A,20231230D,20240117B
  ```
- **Machine and environment:** `profx`, conda env `astro`. Last seen: Python 3.13.2, astropy 7.0.1, numpy 2.2.5, astroquery 0.4.10, dust_extinction 1.7.
  - `frb` is installed editable from this repo. `~/Projects/FRBs` is a symlink to `/mnt/tank/Astronomy/FRBs`, which is fine.
  - `$FRB_GDB` = `/mnt/tank/Astronomy/Galaxy_DB/` (`~/.bashrc` line 30), and `$NEDLVS` = `/mnt/tank/Astronomy/FRBs/NEDLVS_20210922_v2.fits`.
  - `$SPECDB` points to `/u/xavier/...`, which does not exist here; spectra are not needed (`build_ppxf` is off).
- **Build:** `frb/builds/build_hosts.py` (`run(host_row)`) queries the public surveys, merges the literature tables listed in `frb/data/Galaxies/Literature/all_refs.csv` (literature comes last and wins), corrects for Galactic extinction (G23, `dust_extinction>=1.7`), and reads the pPXF, CIGALE and Galfit products from `$FRB_GDB/<Project>/<Reference>/`. A missing pPXF or CIGALE file is skipped silently; a missing Galfit file raises.
- **Rules in the code since GitHub #240** (prompts PR-6, PR-11, PR-12):
  - **Multiple survey matches** within 1″: the **nearest** source is taken. The `public_hosts.csv` column `Survey_match` overrides this per survey (`DECaL:2` = 2nd nearest, `DECaL:skip` = drop the survey). Of the 48 hosts, only **20211127I = `DECaL:2`** is set; the read-only pass (FU-Q16) found no other host that needs one.
  - **Upper limits:** a survey error of `-99` with a real magnitude is written as `_err = 999`; `_flux_err` stays `-99`. No survey `-99` errors should remain.
  - **No measurement:** `-999`. A magnitude at or below `PHOTOM_MISSING_MAX = -90` is dropped. The DELVE `98`/`99` sentinels and DECaL `inf` are masked.
  - **Removed columns:** the survey `separation` column is dropped before the merge, and `Bad_photom` bands are removed.
- **Driver:** `prompts/rebuild_pr_hosts.py`.
  - `--hosts <list>` overrides the default (the 55 #240 hosts). `--retry N` reruns a host that lost bands. `--report <file>` names the report (default `rebuild_report.md`). `--dry-run` checks the setup.
  - **`--keep-morphology <hosts>`** builds those hosts without Galfit and copies the committed `morphology` block into the new JSON (FU-Q10).
  - For every host the driver compares the new JSON with `git HEAD` and flags lost bands, lost `neb_lines`/`derived`/`morphology` keys and changed values. It commits nothing.
- **Known special cases among the 48:**
  - **Galfit files missing (FU-Q4/Q10):** 20180916B (WFC3_F110W), 20190523A (LRIS_r), 20190608B (WFC3_F160W), 20200120E (SDSS_g). Build these with `--keep-morphology`. The first two were probably lost when Galaxy_DB moved.
  - **20190520B:** needs `frb/data/analysis/CIGALE/SOAR_U.dat` (renamed from `SOAR_u.dat` in #240; SDSS u′, FU-Q3).
  - **SDSS outages in Run 1:** 20231005A, 20231017A, 20240117B.
  - **Multiple matches (FU-Q16):** 20211127I DECaL (`DECaL:2`), 20240210A DECaL, 20231206A DECaL, 20231230D DELVE/NSC, 20240310A VISTA/GALEX, 20200906A VISTA, 20210320C VISTA, 20201123A VISTA, 20211212A GALEX, 20230703A GALEX.
  - **Literature 2MASS now read (FU-Q2):** 20201124A and 20211212A (`gordon2023_photom.csv`, lowercase since #240).
- **Helper scripts** (in the old session scratchpad; recreate if needed): `multimatch_pass.py` (the read-only multiple-match pass) and `flag999.py` (the scripted `-99` → `999` edit). Neither should be needed here.
- **Git:** Xavier runs all git write commands (branch, add, commit, push). Claude may run read-only git and `gh` commands, and opens the PR with `gh` only when asked.

## Prompts

1. **Setup check (on `rebuild_hosts`, before building).** Use Opus 5.5. Log your work.
   - Check that GitHub #240 is merged (`gh pr view 240 --json state`), that the current branch is `rebuild_hosts`, and that it contains `origin/main`. Stop if there are local changes you did not expect (the untracked `prompts/rebuild_hosts.md` is expected).
   - From a directory outside the repo, `python -c "import frb; print(frb.__file__)"` must resolve, through `realpath`, to this repo.
   - Check `$FRB_GDB` and `$NEDLVS`, and `dust_extinction` ≥ 1.7.
   - Check that `public_hosts.csv` has the `Survey_match` column with 20211127I = `DECaL:2`, and that `SOAR_U.dat` exists.
   - For the 48 hosts, list every host whose `Galfit_filter` file (`_galfit.fits` or `_galight.json`) is missing. Expect exactly the 4 above; any other host means Galaxy_DB changed, so stop and ask.
   - `python prompts/rebuild_pr_hosts.py --hosts <48> --dry-run` must list 48 hosts with no warnings.
   - Run `pytest frb/tests/test_frbhosts.py frb/tests/test_galaxies.py frb/tests/test_photom.py frb/tests/test_surveys_offline.py frb/tests/test_build.py`.
   - Record the Python, astropy, numpy, astroquery and dust_extinction versions. Do not build anything yet.

2. **Rebuild the 48 hosts.** Use Opus 5.5. Log your work.
   - Run from the repo root, in the background, with the console output kept in a file (expect 2–3 hours):
     ```
     python prompts/rebuild_pr_hosts.py --hosts <48> \
         --keep-morphology 20180916B,20190523A,20190608B,20200120E \
         --retry 2 --report rebuild_report_<date>.md
     ```
   - Post a progress update every 30 minutes (hosts done, retries, failures, survey errors).
   - Rerun any host that failed on a survey outage (SDSS 503 or timeout, MAST, Data Lab, Euclid) alone with `--retry 3`, a few times if needed. Say which survey keeps failing.

3. **Drive** You have suggested that I make links for the CIGALE files.  Instead, I wish to update these on the Google Drive.  Can you do so using `rclone`?  Or create a script for me to do so.  Discuss in Q&A if need be.  Use Opus 5.5. Log your work.  I would then rsync from the Drive to the local drive.

4. **Review every changed JSON.** Use Opus 5.5. Log your work.
   - Compare every key in every changed JSON with `git HEAD`, section by section (`photom`, `derived`, `neb_lines`, `morphology`, `offsets`, `redshift`), as in prompts PR-6 and PR-8 of add_frbs.md.
   - **Expected changes** (from Run 1 and the FU decisions):
     - DELVE sentinels dropped: 20190611B z, 20200430A i, 20210410D z, 20211127I i, 20220105A i, 20220918A g, 20230718A g/r/i.
     - DECaL `inf` dropped: 20210117A r, 20211127I r.
     - Gordon2023 `-999` placeholders dropped (`agn_tau*`, `f_AGN*`), and NaN `UV_global`/`UVff` dropped. `f_AGN` from our CIGALE fits where they exist: 20121102A, 20180924B, 20181112A, 20190102C, 20191001A.
     - Extinction shifts of 0.001–0.04 mag. Up to about 0.09 mag for 20121102A and 20201124A, last built in 2023-07 before G23. `AV_nebular` moves by 1–6%.
     - Offsets follow updated FRB positions, e.g. 20200430A (+1.5″), 20191001A, 20191228A, 20220501C, 20240201A.
     - Survey upper limits as `999`; no `-99` survey errors.
     - Multiple matches:
       - 20211127I DECaL from the 2nd source (g ≈ 15.2, as committed);
       - 20240210A, 20231230D and 20240310A keep their committed survey bands (nearest source);
       - 20231206A **gains** real DECaL g/r/z (about 16.9/16.1/15.4; it had `inf`);
       - GALEX from the nearest source.
     - 2MASS: 20201124A and 20211212A carry the `gordon2023_photom.csv` 2MASS (j/h/k ≈ 17.40/17.02/16.91 and 15.90/15.90/15.39, raw) with `_ref` = Gordon2023.
     - New bands: GALEX, 2MASS, DELVE, VISTA, as in Run 1.
     - The 4 `--keep-morphology` hosts: `morphology` identical to the committed one.
     - 20190520B builds (`SOAR_U`).
   - **Not acceptable:** any other lost photometry band, a lost `neb_lines`, `derived` or `morphology` key, a magnitude shift above about 0.1 mag that the extinction law does not explain, a host at the wrong position, or a failed host.
     - For a flaky survey, rerun the host.
     - For a code or data problem, describe it in the Q&A section with a recommendation.
     - If a host cannot be fixed, ask me to restore its old JSON (`git checkout -- <file>`).
   - Spot-check three hosts by hand and explain every changed key: 20121102A (pre-G23 build), 20211127I (`DECaL:2`, sentinels) and 20231206A (DECaL gained).

5. **Verify, finish the PR text, and hand over.** Use Opus 5.5. Log your work.
   - `FRBHost.by_frb(FRB.by_name(...))` must load all 103 hosts.
   - `python prompts/check_repo_csvs.py` must report only the known items (20190614D, the 20200120E/M81 offset, prochaska2019 row 0).
   - Run `pytest` on the five test files of prompt 1.
   - Fill in the PR draft below: one row per changed host (bands lost, bands gained, other changes), and the hosts left at the old JSON with the reason.
   - List the files to commit (only host JSONs, this file and the dated report) and the git commands for me to run. After I push, open the PR with `gh pr create` only if I ask, and post nothing else.

6. **Follow-ups (optional, after the PR is open).** Use Opus 5.5. Log your work.
   - `frb/data/Galaxies/20190614D/`: a host directory with no `public_hosts.csv` row. Find out what it is and recommend: delete, or add a row.
   - `frb/data/Galaxies/PATH/tmp.csv` (untracked PATH scratch table): confirm it can be deleted.
   - The `-998` lower limit in the figures: `figures/utils.log_me` hides the error bar for `err < 0`, and `figures/galaxies.py` prefers `SFR_nebular` even when it is a lower limit (20181030A, 20220912A). Propose a fix.
   - The CHIME table (`zdm/papers/Mo_Repeaters/py/build_host_table.py`) is not on `profx`. Remind me to rerun it on the laptop and diff `CHIME_FRB_hosts.csv`.

## Q&A

**Q1. GitHub #240 is not merged yet (prompt 1).**
`rebuild_hosts` was created from `repeater_updates` (3fde91da), not from `main`. It contains `origin/main`, but it is 24 commits ahead of it (all of #240), and #240 is still open: approved, with Sunil's re-review requested. Building now would use exactly #240's code, which is what main will have unless the re-review changes it. The rebuilt JSONs stay uncommitted, so they carry over to a new branch.
- **Recommended:** wait for #240 to be merged, then recreate the branch from main. The staged `prompts/rebuild_hosts.md` comes along:
  ```
  git checkout main && git pull
  git checkout -B rebuild_hosts
  ```
  If #240 is squash-merged, recreating the branch avoids carrying its 24 commits into the new PR.
- Alternative: run prompt 2 now on this branch (it takes 2–3 hours), and recreate the branch from main before committing. Only safe if the re-review does not change `build_hosts`, the survey code or `public_hosts.csv`; if it does, rerun prompt 2.

>A.  The PR was merged in.  We are all set 

**Q2. Galaxy_DB products stored under the letterless host names (prompts 2 and 3). Revised 2026-10-21.**
`build_hosts.search_for_file` looks for `HG<FRB>_...` with the TNS letter (e.g. `HG20190523A_CIGALE.fits`). Some Galaxy_DB files still use the pre-letter names (`HG20190523_...`), so the build has not seen them since the hosts were renamed:

| Host | Letterless files (in its `Projects/References` dirs) | Effect |
|---|---|---|
| 20180301A (#240, in main) | `Realfast/Bhandari2021`: CIGALE (3 files), `DEIMOS_ppxf` | no `neb_lines` (Mstar from literature) |
| 20190523A | `Realfast/Bhandari2021`: CIGALE (4 files); `CRAFT/Heintz2020`: `LRIS_r_galfit` | **the rebuild lost all its CIGALE values** (Mstar 4.6e11, SFR_photom 4.19, EBV_photom 0.66, Z_photom, age_mass, f_AGN, M_r, u−r, Lnu_r); morphology was kept by `--keep-morphology` |
| 20191228A | `Realfast/Bhandari2021`: CIGALE (3 files), `LRISr_ppxf` | no `neb_lines`, no Mstar |
| 20201123A | `MeerTRAP/rajwade2022`: `Goodman_ppxf` | no `neb_lines` |
| 20190611B, 20190714A | `Realfast/Bhandari2021`: CIGALE | not copied: these hosts already get CIGALE values from other files, and a lettered copy could change which fit the build picks |

**Done in prompt 3:** `prompts/drive_letter_copies.sh` made server-side **copies** under the lettered names on `GDrive:Astronomy/Research/FRB/Galaxies/Galaxy_DB` (18 files, listed in the script). **Originals, at Xavier's request:** 7 removed (to the Drive trash). 11 could not be removed (403 `insufficientFilePermissions`): the CIGALE files of 20180301, 20190523 and 20191228 and `HG20180301_GMOS-S_r_galfit.fits` in `Realfast/Bhandari2021`, owned by Shivani Bhandari. They are harmless, since the build ignores letterless names. See decision (a).

**Your steps:**
1. Sync these to the local Galaxy_DB, with your rsync or with:
   ```
   for d in Realfast/Bhandari2021 CRAFT/Heintz2020 MeerTRAP/rajwade2022; do
     rclone copy "GDrive:Astronomy/Research/FRB/Galaxies/Galaxy_DB/$d" "/mnt/tank/Astronomy/Galaxy_DB/$d" --include "HG20*A_*" -v
   done
   ```
2. Tell me, and I rebuild 20180301A, 20190523A (now **without** `--keep-morphology`: its Galfit file will be found), 20191228A and 20201123A, then review them.

**Decisions:**
- (a) The 11 remaining letterless originals are owned by Shivani Bhandari, and only she can delete them (or transfer ownership to you). **Recommended:** leave them; the build never reads letterless names. Alternatively, ask her to remove them.
- (b) For 20180301A, 20191228A and 20201123A the pPXF files were **never read** before, so their `neb_lines` and nebular-derived values (e.g. `SFR_nebular`, `AV_nebular`) will be new. Accept them in this PR (**recommended**: they are your fits), or skip the pPXF for now (blank `Spectrum` in `public_hosts.csv` for those hosts)?
- (c) 20180301A is a #240 host already in main. Rebuild it in this PR (**recommended**), or leave it for later?

>A. (b) yes, use them although I suspect they might be junk; (c) Rebuild it

**Q3. 20200120E (M81): survey catalog photometry of a 7th-magnitude galaxy (prompt 4 review).**
The rebuild lost DECaL g/r/z and gained SDSS u–z.
- **DECaL:** within 1″ of `Coord` the only Tractor entry is type `DUP` (all `-99`). The galaxy's DR10 model is a Sérsic 1.7″ away (g 7.92, r 7.02, z 6.95), outside the 1″ search. The committed DECaL (g 7.54, r 7.62, z 7.70, from the 2022-10 build) came from an older catalog and has an implausibly blue g−z = −0.16 for M81's bulge.
- **SDSS** (new): the pipeline photometry of the nucleus (u 11.14, g 9.35, r 8.42, i 7.98, raw), which badly underestimates a ~20′ galaxy.
- **Pan-STARRS** (kept, unchanged apart from extinction): g 9.43, also a nuclear measurement.

None of these is a total magnitude. Nothing in the repo uses M81's integrated photometry for SED fitting (`CIGALE_photom` is blank).
- **Recommended:** set `Survey_match` for 20200120E to `SDSS:skip,DECaL:skip,Pan-STARRS:skip,WISE:skip` (or list those bands in `Bad_photom`), so that no misleading catalog magnitudes are kept. If you want M81 photometry, add a literature table with total magnitudes (e.g. RC3 or the S4G/2MASS LGA) instead. Rebuild the host afterwards.
- Alternative: accept the rebuilt JSON (SDSS and Pan-STARRS nuclear photometry, no DECaL), noting in the PR that these are not totals.

>A.  Use your Recommendation

**Q4. The never-read pPXF fits: are the nebular values usable? (Q2 b follow-up)**
After the Q2 rebuild, three hosts have `neb_lines` for the first time:

| Host | Hα/Hβ (Hβ S/N) | `AV_nebular` | `SFR_nebular` | `SFR_photom` (CIGALE) | Assessment |
|---|---|---|---|---|---|
| 20201123A | 3.74 (25) | 0.75 | 0.24 | — | plausible; [NII]/Hα = 0.82 suggests some AGN/LINER emission |
| 20180301A | **2.18 ± 0.17** (13) | **−0.78** | 1.83 | 3.11 | Balmer decrement 4σ below Case B (2.86), probably Hβ absorption or flux calibration; the negative AV *lowers* the SFR |
| 20191228A | 13 (**0.9**, not detected) | **4.31** | **0.69** | 0.20 | junk: AV from a non-detected Hβ inflates the SFR |

`public_hosts.csv` already has the `Bad_EM_lines` column (used for 20181112A = `Halpha`). Lines listed there are dropped before `AV_nebular` and `SFR_nebular` are computed, so without Hβ the SFR comes from Hα with no dust correction.
- **Recommended:** `Bad_EM_lines` = `Hbeta` for **20191228A** and **20180301A**, and rebuild both (about 5 min). Their `SFR_nebular` then comes from Hα without a dust correction, which is conservative. 20201123A stays as is.
- Alternative: keep all three as built (your Q2 b answer).
- Alternative: drop the pPXF for 20191228A entirely (blank `Spectrum`).

>A.  Use your recommendation

## PR Checklist

- [x] GitHub #240 merged (2026-10-20); `rebuild_hosts` has main's tree (fast-forward to `origin/main` before committing)
- [x] Setup check (prompt 1)
- [x] 48 hosts rebuilt, failures rerun (prompt 2). 48/48 ok; letterless Galaxy_DB files found (Q2)
- [x] Lettered copies on Google Drive (prompt 3); synced; 4 Q2 hosts rebuilt; Q3/Q4 applied
- [x] Every changed JSON reviewed; spot checks (prompt 4, run before prompt 3 was added). Open: Q2 hosts, Q3 (M81)
- [x] Verification; PR text; commit list (prompt 5)
- [ ] PR opened

## PR draft (for `rebuild_hosts`)

Title: Rebuild the host JSONs not touched by #240

Rebuilds the 48 host galaxy JSONs that #240 did not touch, plus 20180301A, with the current code (main after #237 and #240) and `dust_extinction` 1.7. All 103 hosts now share the same survey code, extinction law (G23) and flag conventions (`999` upper limits, `-999` no measurement), and the nearest-source rule with per-host `Survey_match` overrides.

**Summary**
- 49 host JSONs rebuilt; none left at the old JSON. All 103 hosts load, `check_repo_csvs.py` shows only known items, and the tests pass.
- **Extinction:** most magnitudes shift by 0.001–0.04 mag (G23 in `dust_extinction` 1.7 vs 1.5). Hosts last built before 2024-10-21 (Fitzpatrick & Massa 2007) shift by up to −0.12 mag (20180916B, EBV 0.87). Every shift above 0.02 mag was reproduced to 0.001 mag from the old and new corrections.
- **Placeholders removed:** DELVE `98`/`99` sentinels (7 hosts), DECaL `inf` (3), Gordon2023 `-999` and NaN derived values; `f_AGN` from our CIGALE fits where they exist.
- **Multiple survey matches:** the nearest source is used. 20211127I uses `DECaL:2`; 20200120E (M81) skips SDSS/DECaL/Pan-STARRS/WISE.
- **Literature 2MASS** (Gordon2023, lowercase since #240) is now read: 20201124A, 20211212A.
- **Galaxy_DB:** products stored under the pre-TNS-letter names (`HG20190523_...`) were copied to lettered names on Google Drive (`prompts/drive_letter_copies.sh`). 20190523A gets its CIGALE and Galfit back; 20180301A, 20191228A and 20201123A get CIGALE and/or pPXF values for the first time. Hβ is excluded (`Bad_EM_lines`) for 20180301A (Hα/Hβ = 2.18, below Case B) and 20191228A (Hβ not detected).
- **`public_hosts.csv`:** `Survey_match` for 20200120E; `Bad_EM_lines` = `Hbeta` for 20180301A and 20191228A.
- **Offsets** follow updated FRB positions (e.g. 20200430A 0.58″ → 1.10″).

| Host | Bands lost | Bands gained | Other |
|---|---|---|---|
| 20121102A | — | — | max Δmag -0.086 (GMOS_N_i); f_AGN from CIGALE; -999/NaN placeholders dropped; SFR_nebular 0.132→0.13 |
| 20180301A | — | — | CIGALE values added (SFR_photom 3.11, EBV_photom, Z_photom, f_AGN); neb_lines added (pPXF, Hβ excluded via `Bad_EM_lines`); SFR_nebular 3.24 (Hα, no dust correction). A #240 host, rebuilt for its lettered Galaxy_DB files |
| 20180916B | — | — | max Δmag -0.120 (Pan-STARRS_z); f_AGN from CIGALE; -999/NaN placeholders dropped; morphology kept (Galfit file missing) |
| 20180924B | — | — | max Δmag -0.002 (WFC3_F300X); f_AGN from CIGALE; -999/NaN placeholders dropped |
| 20181112A | — | — | f_AGN from CIGALE; -999/NaN placeholders dropped |
| 20190102C | — | — | max Δmag -0.017 (WFC3_F300X); f_AGN from CIGALE; -999/NaN placeholders dropped |
| 20190520B | — | — | max Δmag -0.031 (SOAR_i); -999/NaN placeholders dropped |
| 20190523A | — | — | DECaL g/r/z +0.02 to +0.05 mag (DECaL catalog update since the 2022-07 build; EBV 0.016); CIGALE and Galfit now read from the lettered files (values identical to before) |
| 20190608B | — | — | max Δmag -0.004 (WFC3_F300X); morphology kept (Galfit file missing) |
| 20190611B | DELVE_z (sentinel) | — | max Δmag -0.011 (VLT_FORS2_g); -999/NaN placeholders dropped |
| 20190714A | — | — | max Δmag -0.005 (WFC3_F300X) |
| 20191001A | — | — | max Δmag -0.002 (WFC3_F300X); f_AGN from CIGALE; -999/NaN placeholders dropped; offset 0.95″→0.69″ |
| 20191228A | — | — | CIGALE values added (Mstar 5.4e9, SFR_photom 0.20); neb_lines added (pPXF, Hβ excluded via `Bad_EM_lines`); SFR_nebular 0.028 (Hα, no dust correction); offset 1.44″→1.06″ |
| 20200120E (M81) | DECaL g/r/z, Pan-STARRS g–y, WISE W1–W4 | — | `Survey_match` skips SDSS/DECaL/Pan-STARRS/WISE: catalog photometry of the nucleus, not totals; morphology kept (Galfit file missing) |
| 20200430A | DELVE_i (sentinel) | — | max Δmag -0.002 (VLT_FORS2_g); offset 0.58″→1.10″ |
| 20200906A | — | — | max Δmag -0.002 (VLT_FORS2_g) |
| 20201123A | — | DELVE_i, DELVE_r, VISTA_J, VISTA_Ks | max Δmag -0.018 (WISE_W3); neb_lines added (pPXF) |
| 20201124A | 2MASS_H (uppercase key), 2MASS_J (uppercase key), 2MASS_K (uppercase key) | 2MASS_h, 2MASS_j, 2MASS_k | max Δmag -0.087 (Pan-STARRS_z) |
| 20210117A | DECaL_r (inf) | — | max Δmag -0.003 (VLT_FORS2_g); -999/NaN placeholders dropped |
| 20210320C | — | — | max Δmag -0.004 (VLT_FORS2_g); -999/NaN placeholders dropped |
| 20210410D | DELVE_z (sentinel) | DECaL_g, DECaL_r, DECaL_z | max Δmag -0.020 (DELVE_i); -999/NaN placeholders dropped |
| 20211127I | DECaL_r (inf), DELVE_i (sentinel) | GALEX_FUV, GALEX_NUV | DECaL from the 2nd source (`Survey_match` `DECaL:2`; the nearest is a knot) |
| 20211212A | 2MASS_H (uppercase key), 2MASS_J (uppercase key), 2MASS_K (uppercase key) | 2MASS_h, 2MASS_j, 2MASS_k, GALEX_FUV, GALEX_NUV | max Δmag -0.002 (VLT_FORS2_g); -999/NaN placeholders dropped |
| 20220105A | DELVE_i (sentinel) | — | max Δmag -0.001 (DELVE_g) |
| 20220501C | — | — | offset 0.17″→0.30″ |
| 20220610A | — | — | — |
| 20220725A | — | GALEX_FUV, GALEX_NUV | — |
| 20220918A | DELVE_g (sentinel) | — | — |
| 20221106A | — | — | — |
| 20230203A | — | 2MASS_h, 2MASS_j, 2MASS_k | — |
| 20230222A | — | 2MASS_h, 2MASS_j, 2MASS_k | max Δmag -0.003 (Pan-STARRS_g) |
| 20230311A | — | — | max Δmag -0.007 (DECaL_g) |
| 20230526A | — | — | max Δmag -0.001 (DES_g) |
| 20230703A | — | 2MASS_h, 2MASS_j, 2MASS_k, GALEX_NUV | — |
| 20230708A | — | — | max Δmag -0.003 (DES_g) |
| 20230718A | DELVE_g (sentinel), DELVE_i (sentinel), DELVE_r (sentinel) | — | max Δmag -0.042 (DELVE_z) |
| 20230730A | — | — | max Δmag -0.021 (Pan-STARRS_r) |
| 20230902A | — | — | — |
| 20231005A | — | GALEX_NUV | — |
| 20231017A | — | — | max Δmag -0.005 (Pan-STARRS_g) |
| 20231025B | — | — | max Δmag -0.001 (DECaL_g) |
| 20231206A | — | DECaL g/r/z (real values replace `inf`) | max Δmag −0.004 (Pan-STARRS_g) |
| 20231226A | — | — | max Δmag -0.001 (DELVE_g) |
| 20231230D | — | — | — |
| 20240117B | — | — | max Δmag -0.004 (DES_g) |
| 20240201A | — | — | offset 6.26″→6.08″ |
| 20240210A | — | — | — |
| 20240304A | — | — | max Δmag -0.003 (DELVE_g) |
| 20240310A | — | GALEX_NUV | — |

Not rebuilt: none. Morphology kept from the previous JSON (Galfit file unavailable): 20180916B, 20190608B, 20200120E.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

## Logs

### Run 1 (2026-10-15, Claude Opus 5.5). Superseded: the JSONs were discarded (FU-Q15)

This first attempt ran on `repeater_updates` before the FU-Q fixes. It found the problems that led to the `Survey_match` rule, the 2MASS column rename and the `SOAR_U.dat` rename. The resolution of each is at the end of this entry.

The run was on `profx` in the `astro` env (Python 3.13.2, astropy 7.0.1, numpy 2.2.5, astroquery 0.4.10, dust_extinction 1.7), with `FRB_GDB=/mnt/tank/Astronomy/Galaxy_DB/` set explicitly. Command: `python prompts/rebuild_pr_hosts.py --hosts <44 hosts> --retry 2`.

**Not run (4 hosts; the build would raise).** These have a `Galfit_filter`, but no `_galfit.fits` or `_galight.json` is in `$FRB_GDB`. All 4 committed JSONs have full morphology.

| Host | Galfit_filter | Last built | In the old in-repo `Galaxy_DB` on 2026-10-13 |
|---|---|---|---|
| 20180916B | WFC3_F110W | 2024-01 | found |
| 20190608B | WFC3_F160W | 2024-11 | found |
| 20190523A | LRIS_r | 2024-01 | missing |
| 20200120E | SDSS_g | 2022-10 | missing |

The first two were probably lost when Galaxy_DB moved to `/mnt/tank/Astronomy/Galaxy_DB/` (`F4/mannings2024` is also gone). `/mnt/tank` has no ZFS snapshots.

**Failed (4 hosts):**
- 20190520B: `SOAR_U.dat is not in the Repo`. The curve was `SOAR_u.dat`; the 2023 build probably ran on a case-insensitive file system.
- 20231005A, 20231017A, 20240117B: SDSS SkyServer 503 or connect timeout on all 3 attempts.

**Built (40 hosts).** Every host was checked key by key against `HEAD`. The expected changes are listed in prompt 3 above. The problems found:
1. `build_hosts.run` kept `srvy_tbl[1:]` (line 287, 2024-11-16) when a survey returned more than one source within 1″, dropping the nearest source, the host. This hit 20231230D (DELVE from a neighbour; r/i lost), 20240210A (DECaL lost; the 2nd entry is all `-99`) and 20240310A (VISTA lost).
2. 20211212A lost its 2MASS: the `gordon2023_photom.csv` columns were uppercase `2MASS_J/H/K` (not valid filter names), and the survey query was skipped by `is_inside` (only 1 2MASS object within 1′). `sharma2024_photom.csv` had the same problem for 7 #240 hosts.

**Resolved (add_frbs.md, FU-Q1 to FU-Q18):**
- (1) The nearest source is now the default, with the `Survey_match` overrides. A read-only pass over all 103 hosts found 17 hosts with multiple matches; overrides are needed only for 20211127I (`DECaL:2`) and 20220418A (`Pan-STARRS:skip`, a #240 host).
- (2) The 2MASS columns are lowercase in both tables.
- The `SOAR_U.dat` rename fixes 20190520B.
- The 4 Galfit hosts use `--keep-morphology`.
- The SDSS hosts are rerun.
- All of this is in #240 (commits 6dd70066, fc942a09, 3fde91da).

### Prompt 1: Setup check (2026-10-20, Claude Opus 5.5)

All checks pass **except the first**: GitHub #240 is still open (state OPEN, not merged), so `rebuild_hosts` was branched from `repeater_updates` at 3fde91da rather than from `main` (Q1). Nothing was built.
- **Branch:** `rebuild_hosts`. It contains `origin/main` and is 24 commits ahead of it (the #240 commits). The only local change is the staged `prompts/rebuild_hosts.md` (expected).
- **Editable install:** from `/tmp`, `import frb` resolves, through `realpath`, to `/mnt/tank/Astronomy/FRBs/FRB/frb/__init__.py`.
- **Environment:** `$FRB_GDB` = `/mnt/tank/Astronomy/Galaxy_DB/` (exists); `$NEDLVS` = `/mnt/tank/Astronomy/FRBs/NEDLVS_20210922_v2.fits` (exists, 961 MB).
- **Versions:** Python 3.13.2, astropy 7.0.1, numpy 2.2.5, astroquery 0.4.10, dust_extinction 1.7.
- **`public_hosts.csv`:** the `Survey_match` column is present, with 20211127I = `DECaL:2` (and 20220418A = `Pan-STARRS:skip`, a #240 host). `frb/data/analysis/CIGALE/SOAR_U.dat` exists.
- **Galfit:** among the 48 hosts, exactly the 4 expected have no `_galfit.fits` or `_galight.json`: 20180916B (WFC3_F110W), 20190523A (LRIS_r), 20190608B (WFC3_F160W), 20200120E (SDSS_g).
- **Dry run:** `prompts/rebuild_pr_hosts.py --hosts <48> --dry-run` lists 48 hosts with no warnings. NED-LVS is kept, since `$NEDLVS` is set.
- **Tests:** `test_frbhosts`, `test_galaxies`, `test_photom`, `test_surveys_offline`, `test_build`: 49 passed.

### Prompt 2: Rebuild the 48 hosts (2026-10-21, Claude Opus 5.5)

- **State:** GitHub #240 merged on 2026-10-20 (merge commit f50b77a8). `rebuild_hosts` is at 3fde91da, which is in `origin/main`; `git diff HEAD origin/main` is empty, so the build used main's code exactly. The branch is not at main's merge commit yet: `git merge origin/main` (a fast-forward) before committing.
- **Command:** `python prompts/rebuild_pr_hosts.py --hosts <48> --keep-morphology 20180916B,20190523A,20190608B,20200120E --retry 2 --report rebuild_report_2026-10-21.md` (console log in the session scratchpad). About 2.5 h, 70 attempts in all.
- **Result:** 48 of 48 ok, **Failed: none**. The three hosts that hit SDSS outages in Run 1 (20231005A, 20231017A, 20240117B) built on the first try, so no reruns were needed. No survey failed repeatedly.
- **Multiple matches:** handled as planned. DECaL took the nearest source twice (20240210A, 20231206A) and the 2nd source for 20211127I (`DECaL:2`, three tries); DELVE, NSC, VISTA (×4) and GALEX (×5) took the nearest. `--keep-morphology` copied the committed morphology for all 4 hosts.
- **Flagged for prompt 3 (20 hosts):** mostly the expected categories (DELVE sentinels, DECaL `inf`, Gordon2023 `-999` placeholders, uppercase 2MASS replaced by lowercase). Hosts that retried (3 tries) did so because those lost bands never come back. New beyond Run 1:
  - **20190523A lost its CIGALE values:** Q2, letterless Galaxy_DB files.
  - **20200120E (M81):** DECaL g/r/z (7.5–7.7 mag) lost, SDSS u–z gained. To review in prompt 3.
  - **20180916B:** lost NaN `SSFR`, `UV_SB`, `UV_global` and the `-999` `agn_tau`. Expected (NaN and placeholder values dropped).
  - **20231206A:** DECaL `inf` replaced by real values (g 16.70 corrected), as expected.

### Prompt 4: Review of the 48 rebuilt JSONs (2026-10-21, Claude Opus 5.5). Run while it was still numbered 3, before the Drive prompt was inserted

Every key of all 48 JSONs was compared with `HEAD`, section by section, with thresholds of 0.1 mag for photometry and 5% for derived values. The Q2 symlinks were not made yet, so the four Q2 hosts (20180301A, 20190523A, 20191228A, 20201123A) are **pending a rebuild**; 20180301A is a #240 host and was not part of this run.

**As expected (no action):**
- **DELVE sentinels dropped:** 20190611B z, 20200430A i, 20210410D z, 20211127I i, 20220105A i, 20220918A g, 20230718A g/r/i. **DECaL `inf` dropped:** 20210117A r, 20211127I r.
- **Gordon2023 `-999` placeholders and NaN `UV_global`/`UVff`/`SSFR`/`UV_SB` dropped.** `f_AGN` is now from our CIGALE fits (20121102A 0.019 ± 0.033, 20180924B, 20181112A, 20190102C, 20191001A).
- **Bands gained:**
  - GALEX: 20211127I, 20211212A, 20220725A, 20230703A, 20231005A, 20240310A.
  - 2MASS: 20230203A, 20230222A, 20230703A; and 20201124A and 20211212A, from Gordon2023 now that the columns are lowercase (the uppercase keys are gone).
  - DELVE/VISTA: 20201123A. DECaL: 20210410D, and 20231206A (16.70/15.94/15.34, replacing `inf`; it agrees with Pan-STARRS g 16.68).
- **Multiple matches as configured:** 20211127I DECaL from the 2nd source (g 15.22, as committed); 20240210A, 20231230D and 20240310A keep their survey bands.
- **No survey `-99` errors left;** EBV is unchanged everywhere; no host moved.
- **Offsets follow updated FRB positions:** 20191001A 0.95″ → 0.69″, 20191228A 1.44″ → 1.06″, 20200430A 0.58″ → 1.10″, 20220501C 0.17″ → 0.30″, 20240201A 6.26″ → 6.08″, and small changes elsewhere.
- **Extinction:**
  - Most hosts shift by 0.001–0.04 mag.
  - Shifts above 0.1 mag occur only in **20180916B** (Pan-STARRS i/z/y, −0.10 to −0.12; EBV 0.87).
  - 20121102A (EBV 0.68) and 20201124A (EBV 0.63) shift by up to −0.09 mag.
  - These three JSONs were last built in 2023-07, under the **Fitzpatrick & Massa (2007)** law. G23 replaced FM07 in 6cdcd8cf (**2024-10-21**). Earlier summaries in this session gave 84ec2baf (2026-02), but that commit only updated the `dust_extinction` calls.
  - Recomputing the old correction (`extinction.fm07` over the full filter curve, as the pre-6cdcd8cf code did) reproduces **every** observed shift to 0.001 mag: Pan-STARRS g–y, GMOS-N, MMIRS J, WFC3 F160W, WISE W1/W3.
  - `AV_nebular` and `SFR_nebular` move slightly (20121102A −0.150 → −0.159; 20201124A 1.28 → 1.36), because the line fluxes get the same correction.
- **`--keep-morphology`:** the 4 hosts have morphology identical to the committed one.
- **20190520B** builds with `SOAR_U` (23.23 → 23.22) and loses only Gordon2023 placeholders.

**Not acceptable / open:**
- **20190523A:** lost all its CIGALE values (Q2). Rebuild once the symlinks exist, or restore.
- **20200120E (M81):** catalog photometry of a 7th-magnitude galaxy (Q3).

**Spot checks:**
- **20121102A:** 10 photometric bands shift by +0.006 (g) to −0.087 (i) mag, each exactly the G23−FM07 difference. `AV_nebular` and `SFR_nebular` are as above. 10 placeholder `derived` keys were dropped, and `f_AGN` (0.019 ± 0.033) comes from CIGALE.
- **20211127I:** DECaL `inf` (r) and DELVE i = 98 were dropped. GALEX FUV/NUV (17.74/17.21) are new. DECaL g/z come from the 2nd source (`DECaL:2`), matching the committed values. All other bands shift by ≤ 0.004 mag, and the offset by +0.02″.
- **20231206A:** DECaL g/r/z are now real (16.70/15.94/15.34 ± 0.001, replacing `inf`); the other bands shift by ≤ 0.004 mag.

### Prompt 3: Lettered copies on Google Drive (2026-10-21, Claude Opus 5.5)

Instead of local symlinks (Q2), the lettered files were created on Google Drive, which has no symlinks, so as server-side copies.
- **Drive layout:** the rclone remote `GDrive:` holds Galaxy_DB at `GDrive:Astronomy/Research/FRB/Galaxies/Galaxy_DB`, with the same layout as `/mnt/tank/Astronomy/Galaxy_DB`. Drive had the same letterless files as the local copy, and none of the lettered names.
- **New script `prompts/drive_letter_copies.sh`:** `rclone copyto` for 18 files. 20180301A: CIGALE ×3, DEIMOS pPXF ×2, GMOS-S Galfit (`Realfast/Bhandari2021`). 20190523A: CIGALE ×4 (`Realfast/Bhandari2021`) and LRIS_r Galfit (`CRAFT/Heintz2020`). 20191228A: CIGALE ×3, LRISr pPXF ×2. 20201123A: Goodman pPXF ×2 (`MeerTRAP/rajwade2022`). `--dry-run` shows what would be copied. 20190611B and 20190714A are deliberately left out (Q2).
- **Ran it:** dry run first (18 copies, no errors), then the real run (18 "Copied (server-side copy)"). Verified with `rclone lsl`: all 18 lettered copies exist next to their originals with identical sizes. Nothing on Drive was renamed or deleted.
- **Q2 revised** for the Drive route, with the sync command for Xavier and decisions (a) to (c). Q3 (M81) is unchanged and awaits an answer.
- **Renumbering:** the earlier review log is now "Prompt 4" (it ran while still numbered 3), and the checklist is updated.
- **Next:** Xavier syncs Drive → local and answers Q2/Q3; then I rebuild 20180301A, 20190523A (without `--keep-morphology`), 20191228A and 20201123A, plus 20200120E if Q3 changes its `Survey_match`, and review them.

### Prompt 3, follow-up: removing the letterless originals on Drive (2026-10-21, Claude Opus 5.5)

At Xavier's request ("remove the originals that don't have letters"):
- **New mode** `bash prompts/drive_letter_copies.sh --remove-originals`: for each of the 18 files, it deletes the original only if the lettered copy exists with the same, non-zero size. The copy step runs as a dry run in this mode. Before running it, I checked that nothing in the repo refers to the old names, and that the `GDrive` remote uses rclone's default `--drive-use-trash`, so deleted files go to the Drive trash (recoverable for about 30 days).
- **Removed (7, owned by xavier@ucolick.org):** `Realfast/Bhandari2021/HG20180301_DEIMOS_ppxf.{ecsv,fits}` and `HG20191228_LRISr_ppxf.{ecsv,fits}`, `CRAFT/Heintz2020/HG20190523_LRIS_r_galfit.fits`, `MeerTRAP/rajwade2022/HG20201123_Goodman_ppxf.{ecsv,fits}`.
- **Not removed (11, 403 `insufficientFilePermissions`, owned by shivanibhandari58@gmail.com):** `Realfast/Bhandari2021/HG20180301_CIGALE{,_SFH,_model}.fits`, `HG20180301_GMOS-S_r_galfit.fits`, `HG20190523_CIGALE{,_SFH,_model}.fits`, `HG20190523_CIGALE_photo.dat` and `HG20191228_CIGALE{,_SFH,_model}.fits`. Their lettered copies (owned by Xavier) are in place.
- **Not touched:** the letterless `HG20190611_*` and `HG20190714_*` CIGALE files (no lettered copies; Q2).
- **Drive check (`rclone lsf`):** all 18 lettered files are present. Q2 (a) is revised accordingly.

### Q2/Q3 follow-up: rebuild of 5 hosts (2026-10-21, Claude Opus 5.5)

After Xavier synced Drive → local and answered Q2 ((b) use the new pPXF values, (c) rebuild 20180301A) and Q3 (recommendation):
- **Q3:** `public_hosts.csv` `Survey_match` for 20200120E = `"SDSS:skip,DECaL:skip,Pan-STARRS:skip,WISE:skip"` (quoted, since it contains commas; CRLF kept; one line changed). Rebuilt with `--keep-morphology` (Galfit SDSS_g still missing). The four surveys were skipped and the morphology kept; `photom` now holds only `EBV`, because no other survey or literature table has M81 photometry.
- **Q2:** rebuilt 20180301A, 20190523A (without `--keep-morphology`), 20191228A and 20201123A: 4 of 4 ok, none flagged. The build now finds the lettered files:
  - CIGALE for 20180301A, 20190523A and 20191228A;
  - pPXF for 20180301A (DEIMOS), 20191228A (LRISr) and 20201123A (Goodman);
  - Galfit for 20190523A (LRIS_r) and 20191228A.
- **Results:**
  - **20190523A** is fully restored: all CIGALE values identical to the committed ones (Mstar 4.65e11, SFR_photom 4.19), and the morphology identical apart from float rounding (now from the real Galfit file).
  - **20180301A** gains CIGALE-derived values (SFR_photom 3.11, EBV_photom 0.23, Z_photom, f_AGN, M_r, u−r) and `neb_lines`.
  - **20191228A** gains CIGALE (Mstar 5.4e9, SFR_photom 0.20) and `neb_lines`.
  - **20201123A** gains `neb_lines`.
  - Two of the three pPXF fits give unphysical Balmer decrements: Q4.

### Q4 and prompt 5: Final rebuilds, verification and hand-over (2026-10-22, Claude Opus 5.5)

- **Q4 (recommendation):** `public_hosts.csv` `Bad_EM_lines` = `Hbeta` for 20180301A and 20191228A. The edit was csv-safe: each row was re-serialized byte-identically before the change. Both hosts rebuilt (ok, none flagged); `AV_nebular` is now absent, and `SFR_nebular` comes from Hα without a dust correction: 20180301A 3.24 (CIGALE SFR_photom 3.11), 20191228A 0.028 (SFR_photom 0.20).
- **Verification:**
  - `FRBHost.by_frb` loads all 103 hosts.
  - `check_repo_csvs.py` reports only the known items (20190614D, the 20200120E/M81 1170″ offset, prochaska2019 row 0).
  - `pytest` on the five test files: 49 passed.
- **PR table:** generated from the diffs of all 49 changed JSONs (48 + 20180301A). The magnitude shifts above 0.02 mag were each checked quantitatively:
  - FM07 → G23 for the hosts built before 2024-10-21: 20121102A, 20180916B, 20201124A, 20190520B.
  - `dust_extinction` 1.5 → 1.7 for the later builds: 20230718A DECam z −0.0417, 20230730A Pan-STARRS r −0.0209 and i −0.0157, all matching to 0.0001 mag.
  - The exception is 20190523A, whose DECaL g/r/z are 0.02–0.05 mag *fainter* at EBV 0.016: the DECaL catalog changed since its 2022-07 build.
- **Branch:** Xavier committed `prompts/rebuild_hosts.md` and `rebuild_report_2026-10-21.md` as c7fe5e41 ("wip") and pushed `rebuild_hosts`. The branch is based on 3fde91da, which is in main (main's merge commit f50b77a8 adds no file changes), so no merge is needed; the PR will show only this branch's commits.
- **To commit:** the 49 host JSONs, `frb/data/Galaxies/public_hosts.csv`, `prompts/rebuild_hosts.md`, `prompts/drive_letter_copies.sh`, and `rebuild_report_2026-10-21.md` (with the follow-up rebuilds appended).
