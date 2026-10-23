# Reconcile the FRB Repo with the Website

## Goals

Our goal is to sync the data between the FRB Repo and the data recorded in the
Italian spreadsheet.   When there are discrepancies, we need to identify them
and make a recommendation but a human will need to make the final decision.

## Context

- **FRB Repo Hosts:** all of `public_hosts.csv` (103) 
- **FRB Repo FRBs:** all of `frb/data/FRBs/FRBs_base.csv` 
- **Italian spreadsheets:** CSV files in `/mnt/tank/Astronomy/FRBs/website/` 
    - `FRB Hosts - FRBhosts.csv` :: Host observed properties
    - 'FRB Hosts - FRB_observed_parameters.csv' :: FRB properties
    - 'FRB Hosts - Hosts_derived_parameters.csv' :: Host derived properties
    - 'FRB Hosts - References.csv' :: References


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
- **Git:** Xavier runs all git write commands (branch, add, commit, push). Claude may run read-only git and `gh` commands, and opens the PR with `gh` only when asked.

- **Template:** The file `frb_host_template.json` is the "go-between" the Repo and the Website.  Code in the Repo reads the host and FRB JSON files and produces the web-site ready JSON file that we will provide the Italians.  The code that does this is `build_hosts_web.py`.

## Prompts

### Planning

1.  Please read the goals and context above.  Then discuss with me a plan of action.  
Put your questions in the Q&A/Planning section below.
Use Opus 5.5.  Log your work.

2.  I have updated the Context section to point at additional files and data.  Please
read that and let me know if you have any questions.
Our first task is to build a Markdown doc named `website/link_sheet_and_repo.md`
which describes the mapping between the sheets and the Repo.
I will ask Tarraneh and Matteo to vet this document.
Prepare an additional CSV file that they will modify to confirm the mapping.
If you have any additional questions, put them in the Q&A/Planning section below.
Use Opus 5.5.  Log your work.

## Q&A

### Planning

**Prompt 1 (2026-10-22, Claude Opus 5.5): proposed plan and questions.** Nothing was changed except this file.

#### What a first look shows

The sheet has 161 FRBs (one row each, plus a units row), with host position, E(B-V), z, P_host, morphology, *grizy JHKs* magnitudes with a reference per band, 9 line fluxes and notes. It has **no DM, RM, fluence or FRB localization**.

| Set | Count | Notes |
|---|---|---|
| In the sheet and in `public_hosts.csv` | 102 | The main comparison |
| In the sheet, not in `public_hosts.csv` | 59 | 16 are in `FRBs_base.csv` (no host entry); 43 are not in the repo at all |
| In `public_hosts.csv`, not in the sheet | 1 | 20231020A (CHIME, Leung2025, P_Ox 0.57, no z). Not 20231020B, which is a different FRB |

Spot checks of 4 hosts show that several columns differ **by convention**, not by error. A naive comparison would flag almost everything:
- **Magnitudes:** the repo stores them corrected for Galactic extinction; the sheet seems to store them as observed. Example: 20121102A, Bassa2017 r: repo `GMOS_N_r` = 23.68, sheet 25.46, with E(B-V) = 0.68.
- **Bands:** the sheet has generic bands (one value per band, chosen from one reference); the repo has instrument bands (`GMOS_N_r`, `DECaL_r`, `Pan-STARRS_r`, ...), often several per band.
- **Line fluxes:** the sheet uses units of 1e-17 erg/s/cm² from the cited paper; the repo values come from pPXF on our spectra, corrected for Galactic extinction (and flux scaling). Example: 20121102A Hα, repo 26.1, sheet 6.52.
- **Morphology:** the sheet gives semi-axes a, b, PA and Reff; the repo gives the Galfit/Galight `reff_ang`, `b/a`, `PA`, `n`. Example: 20121102A PA 24 vs 76.
- **Offsets:** the sheet `calc_ang_off` is computed by the sheet maintainers; the repo `ang_best` is computed from the FRB ellipse. Example: 20180301A 2.90″ vs 2.08″.
- **E(B-V):** the sheet uses Schlafly2011 and the repo uses the IRSA `meanValue`. They agree to about 0.01.
- **Host coordinates and z** agree closely in all 4 checks (z 0.3304 vs 0.33044 is rounding).

#### Proposed plan

1. **Column map (no code).** For each sheet column: the repo source (`public_hosts.csv`, `FRBs_base.csv`, host JSON key), units, convention, how to convert, and the tolerance. I write it as a table here for you to approve. Draft tolerances: position 0.5″ (or the larger quoted error); z 1e-4, or 3σ; P_host 0.01; E(B-V) 0.01; magnitudes 0.05 after converting to the same extinction convention; line fluxes 10%.
2. **Membership report.** List the 59 + 1 FRBs above, with the repo status of each (in `FRBs_base` or not, z, P(O|x), the reference in the sheet), and a recommendation: add to the repo, leave out (for example no z or low P_host), or a naming mismatch.
3. **Read-only comparison script** (`prompts/reconcile_website.py`, like `check_repo_csvs.py`). It reads the sheet, `public_hosts.csv`, `FRBs_base.csv` and the committed host JSONs. It writes one row per discrepancy: FRB, quantity, sheet value and reference, repo value and reference, difference, class. The classes:
   - **rounding**: within tolerance; reported in counts only;
   - **convention**: explained by the column map (extinction, units, definition);
   - **missing**: on one side only;
   - **conflict**: a real disagreement.
   It also checks the repo against itself (`public_hosts.csv` z and P_Ox against `FRBs_base.csv` z and `P(O|x)`, and against the JSON).
4. **Recommendations.** For each conflict: the most likely right value and why (newer paper, the reference we cite, a typo or a swapped column). I check the papers only where you ask. You decide.
5. **Apply the decisions** in a later prompt. Repo changes go through the usual path (the `public_hosts.csv` / `FRBs_base.csv` / literature tables, then a rebuild of the affected hosts). Sheet changes become an edit list for the sheet maintainers. You run all git write commands.

#### Questions

- **Q1. Direction.** Is the sheet the source for the website, with the repo checked against it? Or is the repo the reference, so that most fixes become edits to the sheet? Or each case on its merits (my default)?
>A. Neither is the gold-standard.  Each case on its merits
- **Q2. FRB-level data.** The sheet has no DM, RM, FRB coordinates or telescope, so `FRBs_base.csv` can be compared only on z and P(O|x). Does the Google Sheet have another tab (for example the FRBs) that you can export as CSV?
>A. Ignore the FRB level data for now
- **Q3. The 59 sheet-only FRBs.** Is adding them to the repo in scope (for example with the `add_frbs.md` workflow), or is a list with a recommendation enough for now? Do you want the 16 that are already in `FRBs_base.csv` treated first?
>A. We will generate a report before doing any work and we should include these 59 as potential hosts to add to the Repo.
- **Q4. 20231020A.** It is a repo host without a z and is not in the sheet. Expected (the website shows only hosts with a z)?
>A. This is fine
- **Q5. Magnitudes.** Do you agree that the sheet magnitudes are observed (not corrected for extinction)? My plan: match the sheet band to the repo filter with the same reference, apply the repo E(B-V) to the repo value, and compare. Where the sheet reference is not in the repo, report it as **missing**, not as a conflict. Should survey magnitudes (Pan-STARRS, DECaL, ...) that the sheet cites by survey name be matched the same way?
>A. This is a good plan
- **Q6. Line fluxes and morphology.** Since the definitions differ, should these be compared at all? My suggestion: compare line fluxes only where the sheet and the repo cite the same paper (literature values), and compare only morphology quantities defined the same way (PA, Reff in arcsec and kpc). Report the others side by side, unclassified.
>A. we wish to compare all of these properties.
- **Q7. Derived quantities.** The sheet has no stellar mass or SFR. Is the website meant to show repo-derived quantities (Mstar, SFR, Z, age, ...) that the sheet lacks? If so, the deliverable also needs an export of those.
>A. I have added those measurements to the Context now.
- **Q8. Output for the decisions.** A CSV plus a summary in this file? Or a private review page (claude.ai artifact) with one row per conflict, where you choose sheet / repo / other and add a note, and I read the choices back? The page could be shared with the sheet maintainers.
>A. We aren't ready to discuss output yet.
- **Q9. Which sheet copy.** The CSV is dated 2026-10-22 16:36. Is it the current one, and will you re-export it if the sheet changes during the work?
>A. Yes, that is the right one
- **Q10. Branch.** I am on `web_page_reconcile`. Is this the working branch, with the comparison script committed to `prompts/` like the earlier ones?
>A. Yes, that is the right branch

**Prompt 2 (2026-10-23, Claude Opus 5.5): the mapping document.** Written: `website/link_sheet_and_repo.md`, `website/link_sheet_and_repo.csv` (one row per sheet column, 110) and `website/link_sheet_refs.csv` (one row per sheet reference code, 157). Each mapping was tested on the 102 shared hosts; the agreement counts are in the CSV.

Findings that change what prompt 1 said:
- **Magnitudes:** the convention depends on the **reference**, not on the sheet. Sharma2024, Bassa2017 and Niu2022 sheet values are observed. Gordon2023, Heintz2020, Caleb2023 and Bhandari2023 values are extinction-corrected (they equal the repo JSON). The repo literature tables hold observed values, and the build corrects all of them.
- **Pan-STARRS:** in 32 of 39 matches the sheet is 0.2–1 mag fainter than the repo (Kron). The sheet probably uses a different magnitude type (doc Q4).
- **Line fluxes:** literature values mostly agree (11–12 of 12–14 per line). The exception is 20121102A: the repo `tendulkar2017_nebular.csv` holds values 4× the paper's, which looks like an extinction correction made when the table was entered. All other literature nebular tables are copied as they are. So prompt 1 was wrong to say the repo lines are all pPXF.
- **Offsets:** `calc_ang_off` is computed in the sheet from its own host and FRB positions. `Offset_ang` is the published value and matches `ang_best` better than `ang_avg`; the template uses `ang_avg`.
- **Morphology:** the sheet `Semiaxis_a/b` and `PosAngle` are not Galfit quantities (b/a 4/19, PA 2/18). Asked in the doc (Q6).
- **References:** 157 sheet codes vs 84 in the repo. Many differ by year (`Shannon2025`/`Shannon2024`, `Muller2026`/`Muller2025`, `Connor2024`/`Connor2023`) or by case.
- **`build_hosts_web.py`:** it multiplies `eellipse` `a`/`b` by 3600, but they are already in arcsec. It also writes the FRB `refs` list into `DM_ref`, `RM_ref` and `z_ref`, and has no source for `P_cc` or `offset_ref`.

Questions:
- **Q11. Two CSVs.** You asked for one CSV. I added a second (`link_sheet_refs.csv`), because the magnitude matching depends on the reference crosswalk and on each reference's extinction convention. Keep both, or merge them?
- **Q12. Location and sharing.** The three files are in the repo `website/` folder (not committed). How will Tarraneh and Matteo get them: the repo, Google Drive, or a private claude.ai page? If they edit the CSVs in Google Sheets, I will read back an export.
- **Q13. `tendulkar2017_nebular.csv`.** Was it corrected for Galactic extinction on purpose? If so, should other literature tables follow the same rule, or should this one be reverted to the published values?
- **Q14. `build_hosts_web.py`.** Fix the `eellipse` ×3600 bug now (small, separate), or wait until the template is revised after the review?
- **Q15. Reference matching by bibcode.** The repo `all_refs.csv` has DOIs and the `References` tab has ADS links. Should I resolve the year-drift cases through the ADS API (needs an ADS token), or leave them to the reviewers?
- **Q16. Scope of the doc.** It covers the 102 shared hosts and only mentions the 59 sheet-only FRBs. Should it also say which repo files a new host would need (a `public_hosts.csv` row, `FRBs_base.csv`, literature tables)?

## Logs

### Prompt 1: Planning (2026-10-22, Claude Opus 5.5)

Read-only. Read the three sources: the sheet CSV (161 rows plus a units row), `public_hosts.csv` (103) and `FRBs_base.csv` (183). Cross-matched the names (102 in common, 59 in the sheet only, of which 16 are in `FRBs_base.csv`, and 1 in the repo only: 20231020A). Compared 4 hosts (20121102A, 20180301A, 20190608B, 20240210A) with their committed JSONs. This showed the convention differences listed in the Q&A (extinction-corrected magnitudes, line-flux units and source, morphology definitions, offsets, E(B-V) map). Checked in `frb/builds/build_hosts.py` that the photometry and the spectra are corrected for Galactic extinction (lines 391-397, 447-449). Wrote the plan and Q1-Q10 under Q&A / Planning. No repo data, code or git state was changed.

### Prompt 2: Mapping document (2026-10-23, Claude Opus 5.5)

Read the four sheet CSVs (`FRBhosts`, `FRB_observed_parameters`, `Hosts_derived_parameters`, `References`), `frb_host_template.json`, `build_hosts_web.py`, and the 103 host JSONs and their FRB JSONs. Tested each candidate mapping on the 102 shared hosts with scratchpad scripts (`probe.py`, `probe_mag.py`, `probe_refs.py`), in the `astro` env (`/home/xavier/miniconda3/envs/astro/bin/python`; the default `python` has no `frb`). Matched the sheet magnitudes to repo filters by reference, or by survey name, and classified each match as observed or corrected using `frb.galaxies.photom.extinction_correction`. Checked the raw literature tables (`gordon2023_photom.csv`, `sharma2024_photom.csv`, `tendulkar2017_nebular.csv`) against the JSONs.

Wrote `website/link_sheet_and_repo.md`, and generated `website/link_sheet_and_repo.csv` and `website/link_sheet_refs.csv` (scratchpad `make_csv.py`, `make_refs.py`). Added Q11–Q16 under Q&A / Planning. No repo data or code was changed, and nothing was committed.
