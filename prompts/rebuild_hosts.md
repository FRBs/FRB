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

3. **Review every changed JSON.** Use Opus 5.5. Log your work.
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

4. **Verify, finish the PR text, and hand over.** Use Opus 5.5. Log your work.
   - `FRBHost.by_frb(FRB.by_name(...))` must load all 103 hosts.
   - `python prompts/check_repo_csvs.py` must report only the known items (20190614D, the 20200120E/M81 offset, prochaska2019 row 0).
   - Run `pytest` on the five test files of prompt 1.
   - Fill in the PR draft below: one row per changed host (bands lost, bands gained, other changes), and the hosts left at the old JSON with the reason.
   - List the files to commit (only host JSONs, this file and the dated report) and the git commands for me to run. After I push, open the PR with `gh pr create` only if I ask, and post nothing else.

5. **Follow-ups (optional, after the PR is open).** Use Opus 5.5. Log your work.
   - `frb/data/Galaxies/20190614D/`: a host directory with no `public_hosts.csv` row. Find out what it is and recommend: delete, or add a row.
   - `frb/data/Galaxies/PATH/tmp.csv` (untracked PATH scratch table): confirm it can be deleted.
   - The `-998` lower limit in the figures: `figures/utils.log_me` hides the error bar for `err < 0`, and `figures/galaxies.py` prefers `SFR_nebular` even when it is a lower limit (20181030A, 20220912A). Propose a fix.
   - The CHIME table (`zdm/papers/Mo_Repeaters/py/build_host_table.py`) is not on `profx`. Remind me to rerun it on the laptop and diff `CHIME_FRB_hosts.csv`.

## Q&A

(Questions raised while running the prompts above go here, with a recommendation and a `>A.` line for the answer.)

## PR Checklist

- [ ] GitHub #240 merged; `rebuild_hosts` created from `main`
- [ ] Setup check (prompt 1)
- [ ] 48 hosts rebuilt, failures rerun (prompt 2)
- [ ] Every changed JSON reviewed; spot checks (prompt 3)
- [ ] Verification; PR text; commit list (prompt 4)
- [ ] PR opened

## PR draft (for `rebuild_hosts`)

Title: Rebuild the host JSONs not touched by #240

Rebuilds the 48 host galaxy JSONs that #240 did not touch, with the current survey code (main after #237 and #240) and `dust_extinction` 1.7. All 103 hosts now share the same survey code, extinction law (G23) and flag conventions (`999` upper limits, `-999` no measurement), plus the nearest-source rule with per-host `Survey_match` overrides. *(Table to be filled in by prompt 4.)*

| Host | Bands lost | Bands gained | Other |
|---|---|---|---|
| | | | |

Not rebuilt: *(to be filled in)*. Morphology kept from the previous JSON (Galfit file unavailable): 20180916B, 20190523A, 20190608B, 20200120E.

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
