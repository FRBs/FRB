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

8. **Stellar-mass conflicts and non-Prospector masses (Q26, Q38).** We cannot run Prospector (and will not), so there are no new SED fits. Instead:
   - Keep the literature masses as recorded. The method is set by `Mstar_ref` and its suffix (Q9).
   - FRB20191106C: Leung2025b CIGALE (10^9.47) stays, per the Q38 tiers. Record the Chang+2015 value (4.5 ± 1.2 e10, via Ibik2024a) so that the CHIME table can cite it.
   - Sanity-check every non-Prospector mass (CIGALE, NED-LVS, SDSS) with a colour-based M/L estimate from the repo photometry, calibrated on the Gordon2023 Prospector masses. Flag outliers in the Q&A. This is a check only; the estimates are not stored as masses.
   - FRB20181030A: keep the Bhardwaj2021b Prospector value (no refit).
   - Look in `papers/` for any published Prospector mass for these hosts. If there is one, add it per Q26.
   Log your work below.

9. **Regenerate the CHIME table.** Once prompts 1–7 are done:
   - remove the corresponding entries from `LIT_HOSTS`, `LIT_SUPPLEMENT`, `RENAME` and `FIX_181030A` in `zdm/papers/Mo_Repeaters/py/build_host_table.py`;
   - rerun it;
   - diff the new `CHIME_FRB_hosts.csv` against the previous one. Every change should be explainable (e.g. the magnitude source).
   Log your work below.

10. **TODO updates.**  I have completed the TODO items and added the answers below it.  Please update the files as need be.  Use Opus 5.5.  Log your work.

11. **Database** I am about to open a PR for this work.  But before doing so, I need you to compare the `FRBs_base.csv` against the copy on my Google Drive at:
`GDrive:Astronomy/Research/FRB/public_frbs`. Describe any differences in the files in the "PR Checklist/FRBs" section below.  
Similarly, compare the `public_hosts.csv` against the copy on my Google Drive at:
`GDrive:Astronomy/Research/FRB/Galaxies/Galaxy_DB/Public_Hosts.xlsx`. Describe any differences in the files in the "PR Checklist/Hosts" section below.  
Use Opus 5.5. Log your work.   

12. **Sync the `public_frbs` Google Sheet with `FRBs_base.csv`.** Use Opus 5.5. Log your work below.
   - **Before running, choose the FRB20190711A `ee_b` (PR Checklist → FRBs, B6):** `1.281` (user, 2026-10-07) (1.28 as in the repo, or 1.281 as on Drive). If it is 1.281, also change `FRBs_base.csv` (CRLF line endings) and rebuild `FRB20190711A.json` with `build_frbs`.
   - **Target:** Sheet `public_frbs`, file ID `1nNwhYZWOnTcLq6Uv0KJebxMet4NzAnUKW7SFZ6n3GoY` (owner xavier@ucolick.org; shared as Editor with jxp@ucsc.edu). Load the `google-workspace` skill before the first edit. Confirm with `get_file_permissions` that jxp@ucsc.edu is a writer; stop if not.
   - **Back up first:** export the current Sheet with `rclone copy "GDrive:Astronomy/Research/FRB/public_frbs.xlsx" <tmpdir>` and keep the copy until the sync is verified. Drive version history is the second fallback.
   - **Edits** (see PR Checklist → FRBs for the full lists). Keep the existing rows in place and append the new rows at the end, in the same order as `FRBs_base.csv`:
     1. Column 19: set the header to `P(O|x)` and fill it from `FRBs_base.csv` (84 CHIME/KKO rows; blank elsewhere).
     2. Rename FRB20231204A → FRB20190303A and FRB20231128A → FRB20191106C, and set their `refs` to the CSV values (`Leung+2025,burst:FRB2023...`).
     3. FRB20231201A: `z` = 0.1119.
     4. Clear `z` for the 12 rows listed in A4.
     5. Append the 9 new repeater rows (A2) with every column copied from `FRBs_base.csv`.
     6. FRB20190711A `ee_b`: as chosen above.
     - Do **not** rewrite the ra/dec values that differ only by rounding (B7). Leave the `repeater` column in the Sheet's existing TRUE/FALSE style.
   - **Verify:** re-export with rclone and rerun the Prompt 11 comparison (in the `ocean14` env, keyed on `Name`). Expect the same 188 names in the same order and no value differences except the rounding-only ra/dec (B7). Report anything else, then delete the temporary files.
   - Do not change `public_hosts.csv` or `Public_Hosts.xlsx` in this prompt.

13. **Workstation setup and Sheet sync (finish prompt 12).** Use Opus 5.5. Log your work below.
   - **State check.** On branch `repeater_updates`, `git pull`, then confirm that the laptop's last changes are present:
     - `FRBs_base.csv` FRB20190711A `ee_b` = 1.281;
     - `FRB20190711A.json` `eellipse.b` = 1.281;
     - this file contains prompts 12–15 and Q14.

     If any of these are missing, stop: the laptop changes were not pushed.
   - **Tool check.** Use ToolSearch (and ask the user to run `/mcp` if needed) to confirm that a **Google Sheets** editor connector is loaded (`get_spreadsheet`, `get_values`, `update_values` or similar), signed in as an account with write access to Sheet `1nNwhYZWOnTcLq6Uv0KJebxMet4NzAnUKW7SFZ6n3GoY` (jxp@ucsc.edu is a writer). Load the `google-workspace` skill and read its `references/sheets.md` before the first edit.
   - **If the Sheets connector is available:** run prompt 12 as written. The `ee_b` choice is already made (1.281), and the repo already has it.
   - **If it is not:** stop and ask the user to choose between Q14 (b) and (c). Do not fall back on your own.
     - If (b) is chosen, first check that rclone works on this machine (`rclone listremotes`, `rclone lsf "GDrive:Astronomy/Research/FRB/"`).
     - Then test on a **copy**: copy the Sheet with Drive `copy_file`, upload `FRBs_base.csv` over the copy with `rclone copyto --drive-import-formats csv`, and confirm that the copy's file ID is unchanged and its contents match.
     - Only then do the same to the real Sheet, after a backup export.
   - **Verify** as in prompt 12: rerun the Prompt 11 comparison (any env with pandas + openpyxl; `ocean14` on the laptop). Expect no differences except rounding-only ra/dec (PR Checklist → FRBs, B7).
   - Then update PR Checklist → FRBs to say the Sheet is in sync.

14. **Sync `Public_Hosts.xlsx` with `public_hosts.csv`.** Use Opus 5.5. Log your work below.
   - **Target:** `GDrive:Astronomy/Research/FRB/Galaxies/Galaxy_DB/Public_Hosts.xlsx`. This is an uploaded .xlsx file, not a native Google Sheet, so rclone can replace it in place.
   - **Back up** the current file with `rclone copy` to a temp dir.
   - Apply PR Checklist → Hosts, A1–A3, by editing a local copy with openpyxl. Keep the sheet name (`Sheet1`), column order and existing cell values:
     - rename 20231204A → 20190303A and 20231128A → 20191106C;
     - append the 8 new CHIME rows with every column copied from `public_hosts.csv`;
     - add the `Bad_photom` column header after `Comments`, with `DECaL_g,DECaL_r,DECaL_z` for 20231230A.

     Do not overwrite the full-precision P_Ox/z values that differ only by rounding (B4).
   - Upload with `rclone copyto <local> "GDrive:.../Public_Hosts.xlsx"`. Re-download and rerun the Prompt 11 hosts comparison. Expect no differences except B4. Delete the temp files.
   - Then update PR Checklist → Hosts to say the file is in sync.

15. **Pre-PR checks.** Use Opus 5.5. Log your work below.
   - Run the tests that touch this work:

     ```
     pytest frb/tests/test_frbhosts.py frb/tests/test_galaxies.py frb/tests/test_photom.py frb/tests/test_frb.py
     pytest "frb/tests/test_frbsurveys.py::test_sdss" "frb/tests/test_frbsurveys.py::test_galex"
     ```

     Plus `test_build.py` if `$FRB_GDB` and `$NEDLVS` are set on the workstation. Report every failure, and say for each whether it predates this branch (compare with `main`).
   - `git diff --stat main...repeater_updates`: list every file changed and check that each change is covered by a prompt log above. Flag anything unexplained.
   - Draft (do not open) a PR description for `repeater_updates` → `main`, summarizing:
     - **data:** 9 new repeaters, 2 renames, base-table fixes, new hosts and literature tables;
     - **code:** NumPy 2.x fixes, `read_lit_table`, `galex`/`catalog_utils`/`ppxf`/`sdss`/`survey_utils` fixes, `build_hosts` NaN skip and `Bad_photom`, `defs` conventions;
     - the PR Checklist status.

     Put the draft in a new "PR Draft" section of this file for the user to review.

## PR Checklist

### FRBs

Compared `frb/data/FRBs/FRBs_base.csv` (this branch, 188 rows) with `GDrive:Astronomy/Research/FRB/public_frbs` (Google Sheet, modified 2026-05-30, 179 rows; exported with rclone). Same 26 columns, and the common rows are in the same order.

**A. Changes made in this PR (Drive does not have them yet):**
1. **Renamed rows (prompt 3):**
   - FRB20231204A → **FRB20190303A**, `refs` = `"Leung+2025,burst:FRB20231204A"`;
   - FRB20231128A → **FRB20191106C**, `refs` = `"Leung+2025,burst:FRB20231128A"`.
2. **9 new repeater rows (prompt 5):** FRB20180814A, 20190110C, 20200223B, 20190417A, 20220912A, 20240114A, 20240209A, 20190208A, 20181119A. All have `repeater = TRUE` and `telescope = CHIME`, with positions, ellipses, DMs and z from the papers. They are appended at the end of the file.
3. **FRB20231201A** `z`: 0.119 → **0.1119** (prompt 4).
4. **`z` blanked** (prompt 4) for the 12 rows with P(O|x) < 0.9: 20230410A, 20230616A, 20230702A, 20230828A, 20230918A, 20230923A, 20230924A, 20231006B, 20231102A, 20231223D, 20231224A, 20240210C.

**B. Pre-existing differences (already on `main` before this branch):**
5. **`P(O|x)` column.** In the Drive sheet, column 19 has **no header** and is **empty** (0 of 179 values). The repo has `P(O|x)` for the 84 CHIME/KKO rows (added in f28709c8, 2025-03-24). The Drive sheet needs the header and the values.
6. **FRB20190711A `ee_b`:** Drive 1.281, repo 1.28 (changed in d5fbc4b8 "frbs", 2026-02-11, which also updated this row's ra/dec; Drive has the new ra/dec).
7. **Rounding only:** ra/dec in 41 rows (and `ee_a` of 20181112A) differ by ≤ 3e-8 deg (≈ 0.1 mas). The Drive sheet stores full double precision, and the CSV has 7–9 decimals. Not meaningful.

**To sync the Drive sheet** with this PR: apply A1–A4, add the `P(O|x)` header and values (B5), and decide which `ee_b` is right for FRB20190711A (B6).

### Hosts

Compared `frb/data/Galaxies/public_hosts.csv` (this branch, 103 rows) with `GDrive:Astronomy/Research/FRB/Galaxies/Galaxy_DB/Public_Hosts.xlsx` (uploaded xlsx, modified 2025-08-09, 95 rows, one sheet). The 13 shared columns match by name, and the common rows are in the same order.

**A. Changes made in this PR (Drive does not have them yet):**
1. **Renamed rows (prompt 3):** 20231204A → **20190303A**, 20231128A → **20191106C**. Coordinates and other fields are unchanged.
2. **8 new rows (prompt 6)**, all `Projects = CHIME`:

   | FRB | P_Ox | z | References |
   |---|---|---|---|
   | 20180814A | 0.34 (prompt 10) | 0.06835 | michilli2023 |
   | 20190110C | 0.779 | 0.12244 | ibik2024a |
   | 20200223B | 0.899 | 0.06024 | ibik2024a |
   | 20190417A | 1.0 | 0.12817 | moroianu2025 |
   | 20220912A | 0.95 | 0.0771 | ravi2023 |
   | 20240114A | 0.997 | 0.130287 | bhardwaj2025 |
   | 20240209A | 0.99 | 0.1384 | eftekhari2024 |
   | 20190208A | 0.9995 | — | hewitt2024 |
3. **New column `Bad_photom`** (prompt 9, Q12), read by `build_hosts.run`. It is empty except for 20231230A: `"DECaL_g,DECaL_r,DECaL_z"`.

**B. Pre-existing differences (already on `main`):**
4. **Rounding only:** in the repo, P_Ox for 20181112A, 20190611B and 20191001A, and z for 20190102C (0.29117 vs 0.291168), are rounded to 5 digits. Drive has full precision. Not meaningful.

**To sync the Drive file** with this PR: apply A1–A3.

## TODO

- [x] **20180814A localization ellipse (Q1).** Waiting on Dr. Michilli. Table 1 of Michilli+2023 gives σRA = 18″ and σDec = 20″, but Fig. 2 shows about 57″ × 20″. `FRBs_base.csv` has Table 1 for now. Once confirmed, update `ee_a`/`ee_b`/`ee_theta`, rerun `build_frbs` and `build_hosts` for 20180814A, and check `offsets`.
- [x] **20180814A P_Ox (Q2).** Blank in `public_hosts.csv`, because Michilli+2023 did not run PATH on this host. Set it once Dr. Michilli replies.
- [x] **20191106C stellar mass (Q11).** Waiting on the Leung+2025b authors. Leung Table 1 gives log M* = 9.47 (CIGALE), the same as 20231201A in the next row. The colour-M/L estimate (10.3–10.5) and Chang+2015 (10.65, via Ibik2024a) suggest about 10.6. The repo keeps `Leung2025b_CIGALE` 9.47 for now. If the authors confirm a typo, update `leung2025b_cigale_derived.csv` (or add `chang2015_derived.csv`) and rebuild 20191106C.

### TODO anaswers
- Use σRA = 18″ and σDec = 20″ for 20180814A
- P(O|x) = 0.34 for the primary candidate.  Please use that and update files as need be
- Here is the answer from Calvin: "FRB 20191106C = FRB 20231128A has M* = 10.43 +- 0.13 as measured by CIGALE. The 20231201A value is correct" 

## Q&A

**Q1. 20180814A localization ellipse (from Prompt 5).**
Michilli+2023 Table 1 gives σRA = 18″ and σDec = 20″ ("approximately circular"). Their Fig. 2 1σ ellipse is about **57″ (RA) × 20″ (Dec)**. The host is 51.9″ from the FRB: ≈2.5σ with Table 1, but inside the Fig. 2 2σ ellipse, as the text says. `FRBs_base.csv` currently has Table 1 (20″ × 18″, θ = 0). The host JSON offsets use it too (`ang_best_err` = 17.9″). Keep Table 1, or switch to a = 57″, b = 20″, θ = 90° from the figure (and rebuild the FRB and host JSONs)?

>A.  I will ask Dr. Michilli for the correct value.  Hold on for now and come back to it later.  Record this as a TODO item in the TODO section above

**Q2. P_Ox for 20180814A.**
Michilli+2023 did not run PATH on 20180814A. The association rests on the host being the only galaxy in the region with z < z_max (they ran PATH only for 20190303A). I left `P_Ox` blank in `public_hosts.csv`, which gives NaN in `build_table_of_hosts`. Leave it blank, or set a value? If `build_host_table.py` cuts on P(O|x) ≥ 0.9, a blank may drop this host.

>A. Same as above.  Add to TODO

**Q3. P_Ox for 20190110C and 20200223B: which prior on P(U)?**
Ibik+2024a Table 1 gives PATH with P(U) = 0.0 and with P(U) = 0.1: 0.918 / 0.779 (20190110C) and 0.994 / 0.899 (20200223B). I entered the P(U) = 0 values (0.918, 0.994). With P(U) = 0.1, both fall below 0.9. Which do you want?

>A. Use P(U) = 0.1 which means both fall below 0.9

**Q4. P_Ox for 20220912A.**
Ravi+2023 quote only "a 5% false-association probability" (PATH, standard priors, P(U) = 0.5). I entered `P_Ox = 0.95`. Hewitt+2023 (EVN) gives no PATH value. OK?

>A. Ok

**Q5. Hosts with no or inadequate survey photometry.**
- 20190417A: no survey match. Moroianu fn. 39 says the host is not in the DECaLS catalog, so the JSON has no photometry and no `EBV`. Moroianu gives GMOS g, r, i, z = 23.45 ± 0.15, 22.42 ± 0.06, 22.42 ± 0.07, 23.32 ± 0.12 (AB, uncorrected, calibrated to PS1). Ibik+2024b Tab 4 gives r = 21.47.
- 20190208A: r ≈ 27, so no survey photometry, as expected. Hewitt+2024 gives GTC r = 27.32 ± 0.16 (text) vs 27.17 (Tab 2).
- 20240209A: the PS1 catalog r is 17.34 uncorrected (17.07 corrected). Eftekhari gives GMOS r = 16.79 ± 0.02 (14″ aperture, extinction-corrected). The host is a large elliptical, so the PS1 catalog magnitude probably misses flux.

Should prompt 7 add `*_photom.csv` literature tables (e.g. `moroianu2025_photom.csv` with `GMOS_N_*`, `hewitt2024_photom.csv` with the GTC r, `eftekhari2024_photom.csv`)? Which value should the CHIME table use for 20240209A?

>A.  Yes, add those.  And use the 14'' aperture flux for 20240209A

**Q6. 20190208A has no redshift: NaN in the JSON.**
`build_hosts.run` calls `Host.set_z(nan)`, so the JSON has `"z": NaN` (and `z_FRB`, `z_spec`) and `physical`/`physical_err` = NaN. The JSON loads fine, and other repo JSONs also contain NaN tokens, but NaN is not strict JSON. Leave it, or strip the NaN entries after the build (an empty `redshift` dict and no `physical`)?

>A. Leave it

**Q7. DELVE `98 ± 99` placeholders (repo-wide).**
`DELVE_i = 98.0`, `DELVE_i_err = 99` (after the extinction correction) appears in the new 20200223B and 20240114A JSONs. It also appears in 20200430A, 20211127I, 20220105A, 20230222B and 20231223C (7 hosts in all, more for other bands). `surveys/delve.py` (L81) masks only NaN and 99.99, not DELVE's 99 non-detection value. Should I change `delve.py` to set mag = 99 / err = 99 to -999 and rebuild the affected hosts? This may be a separate PR.

>A.  Let's ignore DELVE for now

**Q8. NED-LVS.**
`$NEDLVS` is not set, and there is no LVS file on this machine. `survey_utils.in_which_survey` asserts on it, so I ran the build with NEDLVS removed from `optical_surveys` in my driver script (no repo change). NED-LVS supplies only z, ebv and Mstar (no `defs` photometry bands), so the photometry is unaffected. If you want the hosts built exactly as on your main machine, set `$NEDLVS` and rerun. Also, `$FRB_GDB` here points to a nonexistent Linux path (`/u/xavier/...`), so no CIGALE, pPXF or Galfit products were found. None exist for these hosts anyway.

>A.  Let's not worry about NEDLVS.  If we thought it were important, I could run on my workstation.

**Q9. How to record the stellar-mass method (Leung+2025b CIGALE / NED-LVS masses; Michilli's SDSS mass for 20190303A).**
Since `Mstar_ref` sets the method, a paper that mixes methods breaks that rule. Michilli2023 uses Prospector for 20180814A but quotes an SDSS-collaboration value for 20190303A, and Leung2025b mixes CIGALE and NED-LVS. Options:
- (a) **Recommended:** add a string key `Mstar_method` (Prospector / CIGALE / NED-LVS / SDSS) to `defs.valid_derived`, set it in each literature table, and fall back to the `Mstar_ref` rule when it is absent.
- (b) Separate keys (`Mstar_CIGALE`, `Mstar_NEDLVS`, ...), following the `Mstar_spec` precedent.
- (c) Encode the method in the ref label (e.g. `Leung2025b_CIGALE`).

The Leung masses have no per-object errors, so they would get `Mstar_err = -999`. I have **not** added the Leung masses or the 20190303A row yet. Per Q26, a host that already has an `Mstar` keeps it.

>A. (c)

**Q10. SFR for 20190303A (log SFR = 0.84 ± 0.04, SDSS-collaboration value, method not stated).**
There is no `defs` key for it: `SFR_SED` is Prospector and `SFR_photom` is CIGALE. Options:
- (a) Add `SFR_SDSS`.
- (b) Use the Q9 approach with a generic key plus `_method`.
- (c) Skip it.

>A. Use SFR_photom.  We will figure it out from the reference

**Q11. FRB20191106C stellar mass: keep Leung (tier rule) or switch to Chang+2015?**
The colour-M/L check (Prompt 8 log) puts this host at log M* ≈ 10.3–10.8 in every survey (DECaL, SDSS, PS1): M_r = −21.2, g−r ≈ 0.7–0.8.
- Chang+2015 (via Ibik2024a): **10.65 ± 0.12**, consistent with that.
- Leung2025b CIGALE: **9.47**, about 1 dex low.
- Leung Table 1 also gives exactly 9.47 for 20231201A (the next row). For 20231201A that value *is* consistent with the colour-M/L estimate (9.7), so 20191106C's 9.47 may be a transcription error in Leung.

Options:
- (a) **Recommended:** adopt Chang+2015, i.e. a `chang2015_derived.csv` with ref `Chang2015` (Ibik2024a does not state the method, and Chang+2015 is not in `papers/`), listed after the Leung table so that it wins. This overrides the Q38 tier rule for this one host because of the evidence above.
- (b) Keep Leung (9.47), as agreed in Q38, and cite Chang in the CHIME-table `Refs` (prompt 9).
- (c) Ask the Leung+2025b authors.

>A. (c).  Add this as a TODO that I will follow up on once I get a response from them

**Q12. 20231230A: DECaL r is a shredded fragment.**
The host JSON has `DECaL_r = 19.89` but `Pan-STARRS_r = 15.00` and `DELVE_r = 14.56`. At z = 0.0298, DECaL gives M_r = −15.8; PS1 gives −20.7, which is consistent with Leung's log M* = 10.04. The DECaL Tractor catalog evidently matched a fragment of this bright, nearby galaxy (DECaL g, r, z are all affected). If `build_host_table.py` prefers DECaL magnitudes, the CHIME table will be about 5 mag too faint for this host. Options:
- (a) **Recommended:** drop the DECaL bands for 20231230A, e.g. via a per-host exclusion such as `Bad_EM_lines`. That column is only for emission lines, so this needs a small `build_hosts` addition such as a `Bad_photom` column. Alternatively, delete them by hand after the build.
- (b) Have the CHIME table prefer PS1 for this host.
- (c) Leave it.

Should I also screen all hosts for DECaL vs PS1 disagreements > 1 mag?

>A. (a), and yes do that screening

**Q13. Leftover files in zdm.**
`py/lit_host_photom_cache.json` (tracked in git) is no longer read by `build_host_table.py`. I also left the previous output as `CHIME_FRB_hosts.prev.csv` (untracked) for comparison. OK to `git rm` the cache and delete the `.prev.csv`?

>A. Yes, that is ok

**Q14. Prompt 12 is blocked: no Google Sheets editor connector in this session.**
jxp@ucsc.edu is now a writer on `public_frbs` (checked with `get_file_permissions`). But this session has only the Google **Drive** connector (search, read, create, rename, share), not the **Google Sheets** connector that edits cells in place (`update_values` etc.). Options:
- (a) **Recommended:** turn on (or add) the Google Sheets connector for this chat, then rerun prompt 12. The edits land in the existing Sheet; the link, sharing and formatting are kept.
- (b) Replace the Sheet's contents by uploading `FRBs_base.csv` with rclone (`--drive-import-formats csv`). This first needs a test on a copy to make sure rclone updates the existing file rather than creating a second one. It replaces everything, including formatting, and the old version stays in version history.
- (c) You paste `FRBs_base.csv` into the Sheet by hand. The columns are in the same order; also add the `P(O|x)` header in column 19.

>A. (a) (user, in chat, 2026-10-07). Waiting for the Google Sheets connector to be turned on in the chat; as of this answer it is still not available.

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

### Prompt 5: Add the missing CHIME repeaters to `FRBs_base.csv` (2026-10-04, Claude Opus 5.5)

I appended nine rows to `FRBs_base.csv` (now 188 rows, no duplicate names), in binary mode with CRLF line endings. The file was not sorted, so the rows go at the end.
- All rows: `repeater = TRUE` and `telescope = CHIME`. In this table `telescope` means the discovery or CHIME-selection instrument, as for FRB20180916B, which has an EVN position but is labelled CHIME. `build_host_table.py` selects on `telescope == CHIME`. The localizing instrument is recorded through the refs.
- `ee_a`/`ee_b` are 1σ semi-axes in arcsec, and `ee_theta` is the major axis in degrees E of N (`build_frbs` sets cl = 68).
- Fluence, RM and P(O|x) were left blank.
- Every value below was read from the PDF (text and tables). Figures were rendered at 300 dpi where needed.

| FRB | Position (source) | 1σ ellipse a × b, θ | DM [pc cm⁻³] | z | refs |
|---|---|---|---|---|---|
| 20180814A | 4h22m44s +73°39′52″ (M23 Tab 1) | 20″ × 18″, 0° (σDec = 20, σRA = 18) | 189.4 ± 0.4 (M23 Tab 1) | 0.06835 (M23 Tab 2) | Michilli2023 |
| 20190110C | 249.33, 41.445 deg (CHIME2023 Tab 1) | 21.9″ × 19.7″, 90° (see note) | 221.92 ± 0.01 (CHIME2023 Tab 1) | 0.12244 (I24a Tab 3) | CHIME2023, Ibik2024a |
| 20200223B | 8.265, 28.831 deg (CHIME2023 Tab 1) | 17.5″ × 15.3″, 90° | 202.268 ± 0.007 (CHIME2023) | 0.06024 (I24a Tab 3) | CHIME2023, Ibik2024a |
| 20190417A | 19h39m05.8919s +59°19′36.828″ (Moroianu §3, EVN) | 5.2 × 4.9 mas, 0° | 1378.9 ± 1.4 (Moroianu, mean of B2 and B5) | 0.12817 (Moroianu) | Moroianu2025, Ibik2024b |
| 20220912A | 23h09m04.8988s +48°42′23.9078″ (Hewitt+23 §3.3, EVN final) | 5 × 5 mas | 219.46 ± 0.04 (CHIME; Hewitt+23 §1, Ravi) | 0.0771 (Ravi Tab 2) | Ravi2023, Hewitt2023 |
| 20240114A | 21h27m39.835s +04°19′45.668″ (Bh25 §3.1.2, EVN) | 93 × 28 mas, **153°** (see note) | 527.723 ± 0.042 (Bh25 §3.1.3) | 0.130287 (Bh25) | Bhardwaj2025, Tian2024 |
| 20240209A | 19h19m33s +86°03′52″ (Shah Tab 1, CHIME–KKO) | 2.12″ × 1.08″, 9.54° | 176.49 ± 0.01 (Shah Tab 2, B1) | 0.1384 (Eft Tab 2) | Shah2024, Eftekhari2024 |
| 20190208A | 18h54m11.27s +46°55′21.67″ (Hewitt+24, EVN) | 260 × 260 mas | 580.03 ± 0.14 (Hewitt+24 Tab 1, B2 DM_S/N) | — | Hewitt2024 |
| 20181119A | 12h41m52s +65°07′02″ (M23 Tab 1) | 29″ × 25″, 0° | 364.05 ± 0.09 (M23 Tab 1) | — | Michilli2023 |

**Notes and decisions:**
- **New reference: 20190110C / 20200223B.** Ibik2024a gives no FRB positions, only DMs (Tab 1) and host positions. Its positions come from CHIME/FRB Collaboration 2023 (ApJ 947, 83; arXiv:2301.08762), which was not in `papers/`. I downloaded it from arXiv and saved it as `papers/CHIME_2023_arXiv2301.08762.pdf`.
  - Its Table 1 gives α, δ in degrees, with baseband uncertainties at **90% confidence**: 249.33(1), 41.445(9) and 8.265(8), 28.831(7).
  - I converted these to 1σ assuming per-coordinate Gaussian errors (÷1.645) and treated the RA error as on-sky degrees. This is the more conservative choice; the paper does not say.
  - The DMs used are the CHIME2023 inverse-variance averages (221.92 and 202.268). Ibik Tab 1 quotes 221.6(1.6) and 201.8(4).
  - `CHIME2023` is a new ref label (not in the list above).
- **Ellipse angle: 20240114A.** Bh25 says "rotated by 27°" without a convention. The combined ellipse in Fig. 1 (RA increasing to the left) leans from N toward **W**, so PA = −27° → θ = 153° E of N.
- **Wrong FRB position in prompt 6: 20190417A.** The adopted FRB position is the fitted EVN one, +59°19′36.**828**″. Prompt 6's "19h39m05.8919s +59d19m36.99s" mixes the fitted RA with the Dec of Moroianu's footnote 37 (19h39m05.892s +59°19′36.99″), which appears to be a correlation phase centre. Use 36.828″ in prompt 6.
- **Coarse position: 20240209A.** Shah Tab 1 rounds the centre to 1 s of RA (≈1.05″ at δ = 86°) and 1″ in Dec, which is comparable to the ellipse. The FRB–host offset is 15.7″ (≈38 kpc at z = 0.138), consistent with "outskirts".
- **Paper inconsistency: 20180814A.** **Michilli+2023 contradicts itself.** Table 1 gives σRA = 18″ and σDec = 20″ ("approximately circular"), which is what I entered. But the 1σ (dashed) ellipse in Fig. 2 measures about **57″ (RA) × 20″ (Dec)** against the 30″ scale bar. Its centre agrees with Table 1 (≈4h22m43.4s +73°39′54″).
  - The proposed host (PS1 J042256.01+733940.7) is 51.9″ from the FRB (ΔRA = +50.7″, ΔDec = −11.3″). That is ≈2.5σ with the Table 1 errors, but inside the 2σ ellipse of Fig. 2, as the text says.
  - 20181119A (same table) may have the same problem; it has no figure to check against.
  - **For the user:** keep Table 1 (as entered), or use a ≈ 57″, b = 20″, θ = 90° from the figure?
- **Other offsets (sanity check).** 20190110C: 31.7″ (≈1.4σ; limited by the 0.01° RA precision). 20200223B: 14.2″. 20220912A: 0.5″ (PSO name truncation). 20240114A: 0.15″. 20190208A: 0″ (O4 at the EVN position).
- **Existing JSONs overwritten.** `FRB20190208A.json` (579.4, 3″ ± 0.4″ sys, refs "Astroflash") and `FRB20181119A.json` (190.4799, 65.1119, 0.01″, "Astroflash") were replaced by `build_frbs` output from the new rows. The old sources are unknown. The other seven JSONs are new.

**JSONs.** `build_frbs.main([...9 names...])` wrote all nine. `FRB.by_name` loads each with the ellipse, refs, `repeater = True` and an NE2001 DMISM.

### Prompt 6: Build host JSONs for the new repeaters (2026-10-04, Claude Opus 5.5)

**`public_hosts.csv`.** I appended 8 rows in binary mode with CRLF line endings, each with `Projects = CHIME` and one lowercase reference. (`run()` asserts that `Projects` and `References` have the same number of entries; the ref is only used to look up `$FRB_GDB/<Project>/<ref>/` products.)

| FRB | Host Coord | Source | P_Ox | z | References |
|---|---|---|---|---|---|
| 20180814A | 04h22m56.01s +73d39m40.7s | PS1 J042256.01+733940.7 (M23 Tab 2) | — (Q2) | 0.06835 | michilli2023 |
| 20190110C | 16h37m16.43s +41d26m36.30s | I24a §2.5 | 0.918 (Q3) | 0.12244 | ibik2024a |
| 20200223B | 00h33m04.68s +28d49m52.60s | I24a §2.5 | 0.994 (Q3) | 0.06024 | ibik2024a |
| 20190417A | 19h39m05.82s +59d19m36.7s | **Ibik+2024b Tab 4** (Gemini; Moroianu gives no host centroid) | 1.0 (Moroianu, P_PATH) | 0.12817 | moroianu2025 |
| 20220912A | 23h09m04.848s +48d42m23.760s | PSO J347.2702+48.7066 (Ravi) | 0.95 (Q4) | 0.0771 | ravi2023 |
| 20240114A | 21h27m39.84s +04d19m45.8s | DESI J212739.84+041945.8 (Bh25 Tab 3) | 0.997 (Tian+24 PATH) | 0.130287 | bhardwaj2025 |
| 20240209A | 19h19m24.086s +86d03m39.240s | 289.85036 +86.06090 (Eft Tab 2) | 0.99 (Shah) | 0.1384 | eftekhari2024 |
| 20190208A | 18h54m11.27s +46d55m21.67s | O4 at the EVN position (Hewitt+24) | 0.9995 (Hewitt+24) | — | hewitt2024 |

20190417A: the FRB–host offset is 0.56″ with the Ibik+2024b centroid, which matches Ibik's 0.56 ± 0.06″. (The FRB position is the EVN fit, Dec 36.828″; see Prompt 5.)

**Code fixes (NumPy 2.5 / missing import).** The first build failed for all 8 hosts:
- `frb/surveys/galex.py`: `u` was used (`self.radius.to(u.deg)`) but never imported. Added `from astropy import units as u`.
- `frb/surveys/catalog_utils.py:103` (`match_ids`) and `frb/galaxies/ppxf.py:226`: `np.in1d` was removed in NumPy 2.x. Replaced with `np.isin`, which gives the same result for 1-D inputs.
- `test_frbsurveys.py` with and without the fixes:
  - `test_galex` now passes;
  - `test_in_which_survey` now gets past the `in1d` crash, then fails on its assertion (needs `$NEDLVS`);
  - 7 other survey tests (euclid, nsc, hsc, first, panstarrs, tully, search_all) fail identically either way. They are pre-existing network or version problems.
- `test_frbhosts`, `test_galaxies` and `test_photom` pass.

**Build.** I ran `build_hosts.main([n])` for each of the 8 FRBs through a driver script that removes NEDLVS from `survey_utils.optical_surveys` (Q8). All 8 wrote `Galaxies/<FRB>/FRB<FRB>_host.json`.
- `FRBHost.by_frb` loads all 8.
- `list_of_hosts()` and `build_table_of_hosts()` now give 103 hosts (95 + 8).

**Survey photometry vs the papers.** "Uncorrected" means the JSON value plus A_λ, recomputed with `photom.extinction_correction` and the JSON `EBV`.

| FRB | EBV | Bands | Check |
|---|---|---|---|
| 20190110C | 0.0067 | DECaL grz, PS1, SDSS, GALEX, WISE (19) | DECaL r: corrected 17.992 → uncorrected **18.009** = paper 18.009 ✓ |
| 20200223B | 0.0415 | DECaL, DELVE, PS1, SDSS, GALEX (19) | DECaL r: 15.974 → **16.080** = paper 16.080 ✓ |
| 20220912A | 0.203 | PS1, WISE (9) | PS1 r: 19.113 → **19.656** vs paper 19.65 ✓ |
| 20180814A | 0.1176 | PS1, 2MASS, WISE, GALEX NUV (13) | PS1 r uncorrected 17.03 vs M23 Kron rK 17.15. The 0.12 mag difference is plausibly the PS1 magnitude type. |
| 20240114A | 0.0582 | DECaL, DELVE, PS1, SDSS (16) | SDSS r corrected **21.94** = Bh25 Tab 8 "extinction corrected" 21.94 ± 0.11 ✓. DECaL r corrected 21.69. |
| 20240209A | 0.1011 | PS1, WISE (9) | PS1 r corrected 17.07 vs Eft GMOS 16.79 (14″ aperture, corrected). Catalog vs large aperture on a big elliptical (Q5). |
| 20190417A | — | none | Host not in the DECaLS/PS1 catalogs (Moroianu fn. 39). No photometry and no EBV (Q5). |
| 20190208A | — | none | Expected (r ≈ 27). z is NaN (Q6). |

**Other notes:**
- `DELVE_i = 98 ± 99` (non-detection placeholder) appears in 20200223B and 20240114A. This is a pre-existing `delve.py` issue (Q7).
- Offsets: 20180814A `ang_best` 51.9″ ± 17.9″ (Q1); 20240209A 15.7″ = 39.7 ± 4.9 kpc, consistent with Shah's 40 ± 5 kpc ✓.
- All 8 have an empty `derived` block; prompt 7 adds the literature tables.

### Prompt 7: Literature derived-quantity tables (2026-10-04, Claude Opus 5.5). (completed 2026-10-05)

**Answers applied (Q1–Q8):**
- Q1 and Q2 recorded as TODO items.
- Q3: `P_Ox` set to 0.779 (20190110C) and 0.899 (20200223B), the P(U) = 0.1 values.
- Q5: photometry tables added (below).
- Q4, Q6, Q7 and Q8: no action needed.

**Done:**
- `build_hosts.run` derived-literature loop: it now **skips** a quantity whose value is NaN or -999. Before, a multi-host table wrote NaN into any host that lacked that quantity.
- New filter `GTC_OSIRIS_r`: added to `defs.valid_filters`, with the transmission curve from SVO (GTC/OSIRIS.sdss_r) in `data/analysis/CIGALE/GTC_OSIRIS_r.dat`. A_r ≈ GMOS_N_r.
- `defs.py`: the `Mstar` comment now says the method is set by `Mstar_ref` and lists the Prospector refs.
- New tables in `Galaxies/Literature/`, appended to `all_refs.csv` in this order (CRLF kept). All values were rechecked against the PDFs. dex values are converted to linear with lo/up errors. Rows use the host-JSON ra/dec.
  - `michilli2023_derived.csv`: 20180814A. Mstar 10^10.78 (+0.12/−0.18); SFR_SED 0.316, `_err = 999`.
  - `ibik2024a_derived.csv`:
    - 20190110C: Mstar 2.5e10 (+0.10/−0.17 e10), SFR_SED 0.54 ± 0.04, SFR_nebular 0.1575 ± 0.0006.
    - 20200223B: Mstar 5.6e10 (+1.14/−0.93 e10), SFR_SED 0.59 ± 0.04.
    - 20191106C: SFR_nebular 1.53, `_err = -999` (SDSS fiber, no error given).
  - `moroianu2025_photom.csv`: 20190417A GMOS_N g, r, i, z = 23.45 ± 0.15, 22.42 ± 0.06, 22.42 ± 0.07, 23.32 ± 0.12 (uncorrected). The text's "23.45 ± 15" is taken as 0.15.
  - `moroianu2025_derived.csv`: 20190417A Mstar 10^7.88 (+0.12/−0.14), SFR_nebular 0.19 ± 0.01.
  - `ravi2023_derived.csv`: 20220912A Mstar 10^10.0 ± 0.1, SFR_nebular 0.1, `_err = -998` (lower limit).
  - `hewitt2024_photom.csv`: 20190208A GTC_OSIRIS_r 27.32 ± 0.16 (§3.3 and abstract; Table 2 gives 27.17).
  - `bhardwaj2025_derived.csv`: 20240114A Mstar 10^8.55 (+0.12/−0.14), SFR_nebular 0.061 (+0.004/−0.003).
  - `eftekhari2024_photom.csv`: 20240209A GMOS_N_r, 14″ aperture (Q5). The paper's value is extinction-corrected (16.79), and `build_hosts` corrects again, so the table stores the **uncorrected** value 16.79 + A_r(0.2655) = 17.0555. The build therefore returns 16.79.
  - `eftekhari2024_derived.csv`: 20240209A Mstar 10^11.34 ± 0.01, SFR_SED 0.36, `_err = 999`.
- Not added yet: the Leung2025b masses and the 20190303A Mstar/SFR (Q9, Q10).
- No `ibik2024a_nebular.csv`. Adding the Hα/Hβ fluxes would make the build compute an `AV_nebular` from Ha/Hb (≈7.5), which nobody asked for.

**Q9 (c) and Q10 applied (2026-10-05).** The method is now encoded as a suffix on the ref label. The `defs.py` `Mstar` comment documents this ("A method suffix on the ref overrides this list").
- `leung2025b_cigale_derived.csv` (`Leung2025b_CIGALE`): 20230926A 10.49, 20231011A 9.59, 20231123A 9.42, 20191106C 9.47, 20231201A 9.47, 20231229A 9.87, 20231230A 10.04.
  - Leung Table 1 marks 20231229A and 20231230A "bc": both estimates exist, and the printed value is the higher-priority CIGALE (b) one.
  - Errors are ±0.16 dex, the paper's stated bound for its SED fits ("≤ 0.16 dex"; no per-object values).
- `leung2025b_nedlvs_derived.csv` (`Leung2025b_NEDLVS`): 20230222B 10.19, 20231223C 10.40, ±0.3 dex (stated in the paper). All values are Chabrier.
- `michilli2023_sdss_derived.csv` (`Michilli2023_SDSS`): 20190303A. The host is SDSS J135159.87+480714.2, the third column of M23 Table 3, matching the repo host position. Mstar 10^10.75 ± 0.03 dex; `SFR_photom` = 10^(0.84 ± 0.04) = 6.92 (+0.67/−0.61) (Q10).
- Added to `all_refs.csv`, after the earlier tables.
- None of these hosts had an `Mstar` before, so Q26 does not apply.

**Rebuild.** I rebuilt 18 hosts with `build_hosts.main([n])`, with NEDLVS dropped from `optical_surveys` in the driver. These are the 9 from above plus 20190303A, 20230926A, 20231011A, 20231123A, 20231201A, 20231229A, 20231230A, 20230222B and 20231223C. All 18 succeeded, after two more code fixes:
- **`frb/surveys/sdss.py`**: the spectroscopic query requested `instrument` in `spec_fields`, which now makes SkyServer return an HTML error page. That raised `InconsistentTableError` and stopped 20190303A every time (reproduced: the same query without `instrument` works). `instrument` is not used anywhere, so I removed it.
- **`frb/surveys/survey_utils.is_inside`**: it now also catches `InconsistentTableError` (an unreadable survey response), warns, and treats the survey as not covering the position, as it already does for DALServiceError, ReadTimeout and HTTPError. For 20231011A, SkyServer returned a 403 page for the 1′ footprint query at the FRB position. The git JSON for 20231011A has no SDSS photometry anyway.

**Checks:**
- `derived` blocks are as tabulated, including 20220912A `SFR_nebular_err = -998`, the 999 limits (20180814A, 20240209A `SFR_SED`) and 20191106C `SFR_nebular_err = -999`. 20190208A has no derived block.
- 20240209A `GMOS_N_r` = **16.790** after the build's extinction correction ✓ (stored uncorrected value 17.0555).
- 20190417A: GMOS-N g, r, i, z and `EBV` added. 20190208A: `GTC_OSIRIS_r` 27.32 and `EBV` added.
- New hosts from prompt 6: all other content is identical to git.
- **The 10 previously built hosts** (20191106C, 20190303A and the KKO hosts) also changed outside `derived`:
  - **Magnitude shifts:** −0.0005 to −0.013 mag, bluer bands more. These scale with the unchanged `EBV` (e.g. 20231123A, EBV 0.245: PS1 g −0.013, r −0.010). That is the current G23 extinction law in `photom.extinction_correction`; the hosts were last built 2025-02-27.
  - **New bands:** GALEX (now working after the prompt 6 `galex.py` fix) or 2MASS. No bands were lost.
- Tests: `test_frbhosts`, `test_galaxies`, `test_photom`, `test_frbsurveys::test_sdss` and `::test_galex` pass (12).
- `build_table_of_hosts()`: 103 hosts, 44 with `Mstar`.

### Prompt 8: Stellar-mass conflicts and non-Prospector masses (2026-10-05, Claude Opus 5.5)

**Prompt rewritten.** Prospector cannot be run, so prompt 8 no longer asks for new fits. It now asks to keep the literature masses, record Chang+2015 for 20191106C, run a colour-M/L sanity check on the non-Prospector masses, and look for published Prospector masses.

**Published Prospector masses.** I searched the PDFs in `papers/` for the KKO hosts with CIGALE/NED-LVS masses and for 20190303A. "Prospector" does not appear at all in the CHIME/KKO 2025, Leung2025b, CHIME 2023 or Ibik2024b papers, so there is nothing to add per Q26. FRB20181030A keeps Bhardwaj2021b (Prospector).

**Colour-M/L sanity check (a check only; no values were stored).** Method:
- Bell+2003 log(M/L_r) = −0.306 + 1.097(g−r), with the host's extinction-corrected g and r from the JSON (DECaL, else PS1, else SDSS). The distance modulus is Planck18, with no K-correction, and M_r,⊙ = 4.65.
- Calibrated on the 17 Gordon2023 (Prospector) hosts at z < 0.5: median offset (Prospector − Bell) = −0.51 dex, robust scatter 0.31 dex, range −1.15 to +0.01.
- Caveat: the calibration sample sits at higher z than the CHIME hosts and there is no K-correction, so absolute offsets are uncertain by a few tenths of a dex. The other literature Prospector masses (Ibik, Michilli, Ravi, Bhardwaj, Eftekhari) come out +0.1 to +0.7 dex above the calibrated estimate. Only outliers of about 1 dex are meaningful.

| Host | z | phot | g−r | log M* (repo) | log M* (M/L, cal.) | Δ | ref |
|---|---|---|---|---|---|---|---|
| 20191106C | 0.108 | DECaL | 0.82 | 9.47 | 10.46 (SDSS 10.32, PS1 10.30) | **−0.99** | Leung2025b_CIGALE |
| 20230926A | 0.055 | DECaL | 0.77 | 10.49 | 10.55 | −0.06 | Leung2025b_CIGALE |
| 20231011A | 0.078 | PS1 | 0.37 | 9.59 | 9.48 | +0.11 | Leung2025b_CIGALE |
| 20231123A | 0.073 | PS1 | 0.14 | 9.42 | 9.21 | +0.21 | Leung2025b_CIGALE |
| 20231201A | 0.112 | PS1 | 0.78 | 9.47 | 9.73 | −0.26 | Leung2025b_CIGALE |
| 20231229A | 0.019 | PS1 | 0.46 | 9.87 | 9.28 | +0.59 | Leung2025b_CIGALE |
| 20231230A | 0.030 | DECaL → **PS1** | 0.43 | 10.04 | 8.27 → **9.77** | +0.27 (with PS1) | Leung2025b_CIGALE |
| 20230222B | 0.110 | DECaL | 0.64 | 10.19 | 10.25 | −0.06 | Leung2025b_NEDLVS |
| 20231223C | 0.106 | DECaL | 0.84 | 10.40 | 10.81 | −0.41 | Leung2025b_NEDLVS |
| 20190303A | 0.064 | DECaL | 0.55 | 10.75 | 10.46 | +0.29 | Michilli2023_SDSS |

**Findings:**
- **20191106C:** Leung's 9.47 is about 1 dex below the photometric estimate in all three surveys. Chang+2015's 10.65 (Ibik2024a Tab 3, "previously published by Chang et al. 2015"; method not stated there) agrees. Leung lists the same 9.47 for 20231201A, the next row of Table 1. → **Q11**.
  - For now the repo keeps Leung (Q38 tiers). Chang+2015 (4.5 ± 1.2 e10) is recorded here and in Q11 for the CHIME-table `Refs` (prompt 9).
- **20231230A:** with DECaL the estimate is 1.8 dex below Leung, because DECaL_r = 19.89 is a shredded fragment (PS1_r = 15.00, DELVE_r = 14.56). With PS1 it agrees with Leung. Leung's mass is fine; the repo DECaL photometry is not. → **Q12**.
- All other CIGALE, NED-LVS and SDSS masses are within ±0.6 dex of the estimate, inside the calibration spread. No action is needed.

### Prompt 9: Regenerate the CHIME table (2026-10-05, Claude Opus 5.5)

**Answers applied first:**
- **Q11 (c):** TODO added (follow up with the Leung+2025b authors on 20191106C).
- **Q12 (a): new per-host photometry exclusion.**
  - `build_hosts.run` reads a new `public_hosts.csv` column, **`Bad_photom`** (comma-separated filter names). Those filters, with `_err` and `_ref`, are removed from the merged table before the extinction correction.
  - `public_hosts.csv` gained the column. Every row got a trailing empty field (CRLF kept); 20231230A has `"DECaL_g,DECaL_r,DECaL_z"`.
  - Rebuilt 20231230A: only the 12 DECaL entries were removed, and nothing else changed.
- **Q12 screening.** I compared DECaL/PS1/SDSS/DELVE g, r, z for every host (350 band pairs) and flagged |Δ| > 1 mag.
  - Only **20231230A** is a shred: DECaL is 4.7–5.6 mag *fainter* than PS1 and DELVE in all three bands.
  - The rest go the other way, with DECaL 1.0–1.9 mag *brighter*:
    - 20200120E (M81; PS1 saturated or shredded): g, r;
    - 20220920A: r, z vs PS1;
    - 20221219A: g vs SDSS (faint);
    - 20230307A: z vs PS1, r and z vs SDSS;
    - 20231230D: z vs PS1, g vs DELVE.
  - These look like extended or faint galaxies where the PS1/SDSS catalog magnitudes miss flux. I took no action. For the CHIME table only 20231230A matters (the others are not CHIME hosts).

**`zdm/papers/Mo_Repeaters/py/build_host_table.py`:**
- **Removed:**
  - `RENAME` (the repo now uses the source names);
  - `LIT_HOSTS` and the whole literature-only branch, with `catalog_mag`, the photometry cache and the now-unused imports (`nebular`, `photom`, `survey_utils`, `SkyCoord`, `units`, `warnings`). The script no longer imports `frb` at all, so the stale `astro` install no longer matters;
  - `LIT_SUPPLEMENT`. The only value it held that is not in the repo, the Chang+2015 note for 20191106C, moved to a small `NOTES` dict;
  - `FIX_181030A` and its fallback in `load_host_json`;
  - `LEUNG_L`.
- **Added:**
  - `PROSPECTOR_REFS` now includes Michilli2023, Ibik2024a, Moroianu2025, Ravi2023, Bhardwaj2025 and Eftekhari2024.
  - `REF_SUFFIX` maps the Q9 ref suffixes to tiers: `_CIGALE` → `Photometric:CIGALE:<ref>`, `_NEDLVS` → `Other:NED-LVS:<ref>`, `_SDSS` → `Other:SDSS:<ref>`.
  - An `SFR_photom` with an `_SDSS` ref gets `SFR_type = SDSS`.
  - `GTC_OSIRIS_r` added to `POINTED`.
  - **P(O|x)**: when `FRBs_base.csv` has no P(O|x), the script uses `public_hosts.csv` `P_Ox`. Without this the Q3 values (0.779, 0.899) would never trigger the 0.9 cut.
- Ran `python py/build_host_table.py`: 97 rows, the same FRB list as before. The previous output is saved as `CHIME_FRB_hosts.prev.csv`.

**Diff vs the previous table** (18 rows with value changes, 19 with `Refs` changes). Every change is explained:

| Change | Rows | Why |
|---|---|---|
| All host values blanked; `Refs = host_ignored:P(O\|x)=0.779<0.9` / `0.899` | 20190110C, 20200223B | Q3: P(U) = 0.1 PATH values, read from `public_hosts.csv` |
| Mag 19.89 (DECaL_r) → **14.56 (DELVE_r)** | 20231230A | Q12: shredded DECaL removed. DELVE ranks above PS1 in the Q18 order; PS1_r = 15.00 |
| `Mstar_err` blank → 0.16 / 0.3 dex | 7 Leung CIGALE / 2 NED-LVS hosts | Leung's stated uncertainties, now in the repo tables (prompt 7) |
| Mag −0.002 to −0.010 | 20231011A, 20231123A, 20231201A, 20231229A | Prompt 7 rebuild with the current G23 extinction law |
| Mag 27.164 → 27.172; Band GTC_r → GTC_OSIRIS_r | 20190208A | Extinction now uses the GTC/OSIRIS r curve (was the SDSS_r stand-in) |
| Mag 22.278 → 22.255 | 20190417A | Extinction now uses the repo GMOS_N_r curve and the repo EBV (was the GMOS_r stand-in) |
| SFR_err 0.6372 → 0.6381 | 20190303A | Repo lo/up errors symmetrized, instead of the log-derivative approximation |
| `Mstar_source` / `Refs` labels (`Michilli+2023` → `Michilli2023`, `Leung+2025(ApJL991,L25)` → `Leung2025b`, ...) | the literature hosts | Labels now come from the repo refs. `z:`/`mag:` refs now say `hostJSON`/`repo` instead of `catalog`. `name:repo_as_...` dropped (renamed in prompt 3) |

The values for 20180814A, 20220912A, 20240114A and 20240209A are unchanged; the repo photometry matches the old catalog/literature magnitudes. 20181119A is still a repeater row with every host column blank.

### Prompt 10: TODO updates (2026-10-06, Claude Opus 5.5)

All three TODO items are closed (checked above), per the answers under "TODO answers".

1. **20180814A ellipse.** Keep Michilli+2023 Table 1 (σRA = 18″, σDec = 20″). `FRBs_base.csv`, `FRB20180814A.json` and the host offsets already use it, so nothing changed.
2. **20180814A P(O|x) = 0.34** (primary candidate, from Dr. Michilli). Set `P_Ox = 0.34` in `public_hosts.csv` (CRLF kept). `FRBs_base.csv` P(O|x) stays blank, as for the other prompt-5 rows; `build_host_table.py` falls back to `public_hosts.csv`. No host rebuild was needed, because `P_Ox` is not stored in the host JSON.
3. **20191106C stellar mass.** C. Leung (priv. comm.): "FRB 20191106C = FRB 20231128A has M* = 10.43 ± 0.13 as measured by CIGALE. The 20231201A value is correct."
   - Edited the 20191106C row of `leung2025b_cigale_derived.csv` (text edit, so no other row changed): log M* 10.43 ± 0.13 → Mstar = 2.69e10 (+9.39e9/−6.96e9). The ref stays `Leung2025b_CIGALE`.
   - 20231201A (9.47) is unchanged.
   - Rebuilt 20191106C with `build_hosts.main`; only `Mstar`, `Mstar_loerr` and `Mstar_uperr` changed in the JSON.
   - The new value agrees with the colour-M/L estimate (10.3–10.5) and is 0.2 dex below Chang+2015 (10.65).

**CHIME table (zdm).** In `build_host_table.py`, `NOTES['FRB20191106C']` now records the correction (`Mstar:Leung2025b_Table1_typo_corrected(C.Leung,priv.comm.)`) instead of the old Chang `Mstar_alt` note. I reran it: 97 rows. The diff vs the previous output has exactly 2 rows:
- **FRB20180814A:** `P_Ox` blank → 0.34. Per Q41 (zdm doc), no P(O|x) cut is applied now, so the host values stay.
- **FRB20191106C:** `Stellar_Mass` 9.47 → 10.43, `Mstar_err` 0.16 → 0.13, and the Refs note as above.

### Prompt 11: Compare with the Google Drive copies (2026-10-07, Claude Opus 5.5)

- The Google Drive connector (signed in as jxp@ucsc.edu) could not find `Astronomy/Research/FRB` (My Drive › Astronomy › Research contains only DISCO and GitHub). `GDrive:` turned out to be the user's **rclone** remote. I downloaded `public_frbs` (Google Sheet → xlsx) and `Galaxies/Galaxy_DB/Public_Hosts.xlsx` with `rclone copy` into a temporary directory, then deleted them. The download was read-only; nothing on Drive was changed.
- `astro` has no `openpyxl`, so I used the existing `ocean14` env (pandas + openpyxl) for the comparison and installed nothing.
- Compared by key (`Name` / `FRB`): rows present only on one side, columns, row order, and every shared cell (strings exactly; numbers exactly, then with a tolerance to separate real changes from rounding). Each difference was checked against `main` with `git log -S` to label it as this PR's change or pre-existing.
- Results are in **PR Checklist → FRBs / Hosts** above. Every value difference is either a change made in prompts 3–10 or pre-existing. The pre-existing ones are: the empty, unheaded `P(O|x)` column in the Drive sheet; FRB20190711A `ee_b` 1.281 vs 1.28; and sub-mas or 5-digit rounding.

### Prompt 12: Sync the `public_frbs` Sheet (2026-10-07, Claude Opus 5.5). BLOCKED (Q14)

- **FRB20190711A `ee_b` = 1.281** (user's choice). Changed in `FRBs_base.csv` (binary edit, CRLF kept; one field) and rebuilt `FRB20190711A.json` with `build_frbs`. The JSON diff is `eellipse.b` 1.28 → 1.281, plus DMISM in the 15th decimal.
- `get_file_permissions` on the Sheet (`1nNwhYZWOnTcLq6Uv0KJebxMet4NzAnUKW7SFZ6n3GoY`) confirms that **jxp@ucsc.edu is a writer**.
- **Stopped before any Drive change.** This session has no Google Sheets editor tools (`update_values` etc.); a tool search found only Drive, Docs and Calendar tools. Drive alone cannot edit an existing Sheet in place. I did not make the backup copy, because no edit was attempted. See Q14.
- Once the Sheets connector is on, prompt 12 can run as written. B6 then disappears, since repo and Drive will both have 1.281.
