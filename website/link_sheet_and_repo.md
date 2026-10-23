# Mapping between the FRB Hosts sheets and the FRB repo

Draft of 2026-10-23, for review by Tarraneh and Matteo.

This document maps each column of the FRB Hosts Google Sheet to where the same quantity lives in the FRB repo (<https://github.com/FRBs/FRB>), and says how to convert one to the other. Its purpose is to agree on the mapping **before** we compare the values. The comparison itself, and the decisions on each discrepancy, come later.

## How to review

There are two CSV files next to this document. Please edit them, not this document:

- **`link_sheet_and_repo.csv`**: one row per sheet column (110 rows). Fill in `confirm` (Y or N). If N, write the right mapping in `correct_mapping`. Add anything else in `comment`.
- **`link_sheet_refs.csv`**: one row per reference code used in the sheet (157 codes). Fill in `confirm`, and if needed `correct_repo_code` and `correct_mag_convention`. We need this file because the magnitudes are matched by reference (section 5).

Leave every other column as it is. The "agreement" numbers in the CSV come from a first test on the **102 FRBs that have a host in both the sheet and the repo**. They are there to show whether a mapping is plausible, not to report discrepancies yet.

The questions where we most need your help are collected in section 10.

## 1. The two sides

**The sheet** (exported as CSV on 2026-10-22/23). Each tab has a header row, a row of units, and one row per FRB (161 FRBs).

| Tab | Content |
|---|---|
| `FRBhosts` | Host position, E(B-V), z, P_host, morphology, *grizy JHKs* magnitudes, 9 line fluxes, notes |
| `FRB_observed_parameters` | FRB position and error ellipse, DM, DM_MW, RM, offsets, repeater flag, notes |
| `Hosts_derived_parameters` | D_L, scale, log M*, SFR (SED and nebular), gas metallicity, notes |
| `References` | Reference code, ADS link, title |

**The repo.** The FRB name is `FRB20121102A` in the repo and `20121102A` in the sheet (and in the host folder names).

| Source | Path in the repo | Content |
|---|---|---|
| `FRBs_base.csv` | `frb/data/FRBs/` | One row per FRB (183): position, ellipse, DM, RM, z, repeater, refs, P(O\|x) |
| FRB JSON | `frb/data/FRBs/FRB<name>.json` | Built from `FRBs_base.csv`. Adds `DMISM` (NE2001) |
| `public_hosts.csv` | `frb/data/Galaxies/` | One row per host (103): host position, P_Ox, z, references, build options |
| host JSON | `frb/data/Galaxies/<name>/FRB<name>_host.json` | Built from `public_hosts.csv`: `photom`, `redshift`, `offsets`, `morphology`, `neb_lines`, `derived` |
| Literature tables | `frb/data/Galaxies/Literature/` | Published photometry, line fluxes and derived values, listed with a DOI in `all_refs.csv` |
| Template | `frb/data/Galaxies/frb_host_template.json` | The fields of the website JSON, written by `frb/builds/build_hosts_web.py` |

When the host JSON is built, survey photometry (Pan-STARRS, DECaL, DES, DELVE, SDSS, 2MASS, WISE, VISTA, ...) is queried first. The literature tables are then merged in, and they win over the surveys. Every magnitude is then corrected for Galactic extinction.

**Which FRBs.** 102 FRBs are in both. 59 are only in the sheet: 16 of them are in `FRBs_base.csv` but have no host in the repo, and 43 are not in the repo at all. They will be listed as hosts we could add to the repo. One repo host, 20231020A, has no redshift and is not in the sheet, which is expected.

## 2. Conventions

| | Sheet | Repo |
|---|---|---|
| No value | empty; `-999`; `-99` in `Hosts_derived_parameters`; `999` in some magnitude columns | missing key; `-999` |
| Upper limit | `999` in some error columns? (to confirm) | magnitude error `999`; `SFR_nebular_err = -998` for a lower limit |
| Errors | one symmetric error, often `0` when unknown | `_err`, or `_loerr` and `_uperr` (asymmetric, in linear units for M*, SFR and Z) |
| References | one code per quantity (`ref_...` columns) | per value for literature photometry, lines and derived values (`_ref` keys); one list per FRB or host otherwise |
| Line fluxes | 1e-17 erg/s/cm² | erg/s/cm² |
| Stellar mass, metallicity | log10 M*/M_sun; log(Z/Z_sun) | linear M*/M_sun; linear Z/Z_sun |
| Cosmology | "Planck" | Planck18 (astropy) |
| E(B-V) | Schlafly & Finkbeiner 2011 | the same map, from IRSA (`meanValue`) |

Please confirm the meaning of the sheet sentinels (`999`, `99`, `-99`, `-999`, `0`, `?`). A misread sentinel would show up later as a false discrepancy.

## 3. Tab `FRBhosts`

The full detail is in the CSV. Here, "agree" is the number of the 102 shared hosts within the proposed tolerance.

| Sheet column | Repo | Conversion and notes | Agree |
|---|---|---|---|
| `GalID` | none | The repo has no galaxy ID. The template has `galID`, but nothing fills it | |
| `RA_deg`, `Dec_deg` | host JSON `ra`, `dec` | Within 0.5″ | 94, 92 / 102 |
| `err_RA`, `err_Dec` | host JSON `positional_error` | Only 8 repo hosts have this | |
| `ref_coords` | `public_hosts.csv` `References` | No reference per quantity in the repo | |
| `calc_ang_off` | host JSON `offsets/ang_best` | The sheet computes this from its own host and FRB positions (85/102 within 0.05″). So it follows those columns | 60 / 102 within 0.1″ |
| `calc_phys_off` | `offsets/physical` | The repo uses `ang_best` and Planck18 | 60 / 101 |
| `E_B-V` | `photom/EBV` | Same map; within 0.01 | 94 / 101 |
| `z` | host JSON `redshift/z` = `public_hosts.csv` `z` = `FRBs_base.csv` `z` | Within 1e-4 | 93 / 101 |
| `err_z`, `ref_z` | none | The repo stores no redshift error or redshift reference | |
| `is_zphot` | `redshift` has `z_phot` but no `z_spec` | All 103 repo hosts have `z_spec`. Of the 8 sheet `y` rows, 1 is a shared host (20240304A) | |
| `P_host`, `ref_P_host` | `public_hosts.csv` `P_Ox`; `FRBs_base.csv` `P(O\|x)` | Within 0.01. `FRBs_base.csv` has P(O\|x) for only 29 of these hosts | 81 / 101 |
| `Semiaxis_a`, `Semiaxis_b` | **unknown** | See section 7 | |
| `PosAngle` | `morphology/PA` | See section 7 | 2 / 18 |
| `Reff`, `Reff_kpc` | `morphology/reff_ang`, `reff_kpc` | Galfit or Galight; 19 repo hosts. See section 7 | 6 / 18, 11 / 18 |
| `g-mag` ... `K_s-mag` | `photom/<filter>` | See section 5 | |
| `OII_flux` ... `SII6731_flux` | `neb_lines/<line>` | See section 6 | |
| `Notes_host` | `public_hosts.csv` `Comments` | Free text; not compared | |

## 4. Tabs `FRB_observed_parameters` and `Hosts_derived_parameters`

| Sheet column | Repo | Conversion and notes | Agree |
|---|---|---|---|
| `FRB_RA_deg`, `FRB_Dec_deg` | FRB JSON `ra`, `dec` (= `FRBs_base.csv`) | Within 0.5″ | 93, 88 / 102 |
| `FRB_semi_a`, `FRB_semi_b`, `FRB_PA` | FRB JSON `eellipse/a`, `b`, `theta` | Arcsec on both sides. The repo keeps the statistical (`a`) and systematic (`a_sys`) parts separately; see Q5 | a: 64; b: 56; PA: 79 / 102 |
| `DM`, `err_DM` | FRB JSON `DM`, `DM_err` | Within 0.5 pc/cm³ | 77 / 102 |
| `DM_MW` | FRB JSON `DMISM` | Both NE2001; within 1 pc/cm³ | 58 / 102 |
| `RM`, `err_RM` | FRB JSON `RM`, `RM_err` | 33 sheet RMs are not in the repo | 19 / 28 |
| `ref_FRB_coord`, `ref_DM`, `ref_RM` | FRB JSON `refs` | One list per FRB in the repo | |
| `Offset_ang`, `err_Offset_ang` | host JSON `offsets/ang_best` | This is the literature value (`ref_Offset`); the repo computes its own offset. The template uses `ang_avg`, but `ang_best` agrees better (24 vs 12 of 49) | 24 / 49 within 0.1″ |
| `Offset_lin`, `err_Offset` | `offsets/physical` | | 35 / 56 |
| `Repeater` | FRB JSON `repeater` | y/n = true/false | 102 / 102 |
| `D_L_Planck`, `Scale_at_z` | none (computed) | Planck18 at the sheet z agrees to 0.3% (median). Computed, so not compared | 98, 95 / 101 within 1% |
| `logMstar`, `err_logMstar` | `derived/Mstar`, `Mstar_loerr`, `Mstar_uperr` | log10 of the linear repo value. Proposed error: (log10(M+up) − log10(M−lo))/2 | 30 / 36 within 0.05 dex |
| `Mstar_UL` | none | `?` in all 137 filled rows | |
| `SFR-SED` | `derived/SFR_SED` | | 51 / 57 within 5% |
| `SFR_Nebular` | `derived/SFR_nebular` | From Hα (or Hβ). The repo corrects for dust when it has A_V from Hα/Hβ | 14 / 23 within 10% |
| `SFR-Neb_UL-Dust` | `derived/AV_nebular`; `SFR_nebular_err = -998` | Meaning to confirm (Q8) | |
| `Z-gas_emiss` | `derived/Z_gas` | log10 of the repo value (repo is linear Z/Z_sun) | 21 / 22 within 0.05 dex |
| `Notes_FRB`, `Notes_param` | none | Free text | |

## 5. Photometry

This is the hardest part of the mapping, for three reasons.

**Bands.** The sheet has one value per generic band (*g r i z y J H K_s*). The repo has one value per instrument filter (`GMOS_N_r`, `DECaL_r`, `Pan-STARRS_r`, `VLT_FORS2_R`, `2MASS_j`, ...), often several per band. We propose to pick the repo filter from the sheet reference:
- **A paper** (for example `Gordon2023`): the repo filter in that band whose `_ref` is the same paper, after the reference crosswalk (`link_sheet_refs.csv`).
- **A survey** (`PanSTARRS`, `DES-DR2`, `Legacy Survey`, `SDSS`, `2MASS`, `2MASS-XSC`): that survey's filter in the repo. Survey values in the repo carry no `_ref`.

With this rule, 250 sheet magnitudes on the shared hosts have a repo match: 208 by paper and 42 by survey. Another 155 have a repo filter in the same band but no matching reference; for them, the sheet reference is not in the repo, or its code differs. The last 64 have no repo filter in that band.

**Galactic extinction.** The repo literature tables store **observed** magnitudes. The host JSONs store magnitudes **corrected** for Galactic extinction (Gordon et al. 2023 curve, R_V = 3.1). The sheet follows each paper, so its convention depends on the reference:

| Sheet reference | Sheet value equals | Evidence (number of magnitudes) |
|---|---|---|
| Sharma2024, Bassa2017, Niu2022 | observed | 75 observed / 5 corrected (Sharma2024); 4/0; 1/0 |
| Gordon2023, Heintz2020, Caleb2023, Bhandari2023 | corrected | 31/1; 7/0; 6/0; 1/0 |
| PanSTARRS | unclear (see below) | 4 observed, 2 corrected, 32 neither |

Example: 20121102A, Bassa2017 *r*: sheet 25.46, repo literature table 25.46, repo JSON 23.68 (E(B-V) = 0.68). The proposed convention for each reference is the `proposed_mag_convention` column of `link_sheet_refs.csv`. Please confirm or correct it there.

**PanSTARRS.** In 32 of the 39 matches, the sheet value is 0.2–1 mag **fainter** than the repo value, even after the extinction correction. The repo uses the Pan-STARRS **Kron** magnitude. The sheet may use the PSF magnitude, or a different catalogue. See Q4.

**Website.** `build_hosts_web.py` currently writes only an *r* magnitude, from DECaL or else Pan-STARRS, and it is extinction-corrected. Which bands, and which convention, the website should show is a decision for later.

## 6. Line fluxes

| Sheet | Repo `neb_lines` | Agree (10%) |
|---|---|---|
| `Halpha_flux` | `Halpha` | 11 / 14 |
| `Hbeta_flux` | `Hbeta` | 11 / 12 |
| `OIII5007_flux` | `[OIII] 5007` | 12 / 14 |
| `OIII4959_flux` | `[OIII] 4959` | 1 / 1 |
| `NII6583_flux` | `[NII] 6584` | 11 / 14 |
| `NII6549_flux` | `[NII] 6548` | no overlap |
| `SII6717_flux`, `SII6731_flux` | `[SII] 6716`, `[SII] 6731` | 0 / 1 each |
| `OII_flux` | `[OII] 3726` + `[OII] 3729` (the repo keeps the doublet separately) | no overlap |

The repo values are in erg/s/cm². Multiply by 1e17 to get the sheet units. The repo has line fluxes for only about 15 hosts; the sheet has them for about 45. A repo value comes from one of two places:
- a literature nebular table, when the line has a `_ref` key. The build copies these values as they are in the table, with no extinction correction;
- the repo's own pPXF fit of the spectrum, corrected for Galactic extinction (no `_ref` key).

The values from literature tables usually agree. One exception is 20121102A Hα, where both sides cite Tendulkar2017: the repo has 26.1 and the sheet 6.52. The sheet notes that these fluxes are not corrected for Galactic absorption. The repo table seems to have been corrected when it was entered: the factor of 4 matches A(Hα) ≈ 1.5 mag at E(B-V) = 0.68. Please confirm (Q7) whether the sheet fluxes are as published or corrected for Galactic extinction.

## 7. Morphology

The repo morphology comes from Galfit or Galight, for 19 hosts: `reff_ang` (arcsec), `reff_kpc`, Sérsic `n`, `PA`, `b/a`. The sheet has `Semiaxis_a`, `Semiaxis_b`, `PosAngle`, `Reff` and `Reff_kpc`. On the 18–19 hosts in both:
- `Semiaxis_b / Semiaxis_a` does not match the Galfit `b/a` (4/19), and `PosAngle` does not match `PA` (2/18; median difference 58°). So `Semiaxis_a/b` and `PosAngle` are probably not a Sérsic fit. What are they (an aperture, an isophote, the PATH galaxy size)? Is the PA measured east of north?
- `Reff_kpc` matches better (11/18) than `Reff` (6/18). The two sheet columns may come from different sources.
- The repo also has `derived/reff_iso` (isophotal radius, Mannings2021, 8 hosts). It does not match `Reff` either (1/8).

## 8. Reference codes

The sheet uses 157 reference codes (132 of them in at least one `ref_` column) and the repo 84. The same paper often has a different code on the two sides:
- **Year:** sheet `Shannon2025`, repo `Shannon2024`. Likewise `Muller2026`/`Muller2025`, `Connor2024`/`Connor2023`, `Law2024`/`Law2023`, `Bhandari2022`/`Bhandari2021`, `Gordon2025`/`gordon2024`. These are probably the arXiv year vs the publication year.
- **Case and suffix:** `gordon2023`/`Gordon2023`, `Ibik2024`/`Ibik2024a`, `Hewitt2024b`/`Hewitt2024`. The repo is not consistent about case either.
- **Survey names:** `PanSTARRS` (also spelled `PanSTARRs`, `PanSTARSS`, `PannSTARRS`), `Legacy Survey`, `DES-DR2`, `2MASS-XSC`, `SkyMapper`. These map to repo survey names. The repo has no SkyMapper photometry.
- **Not references:** `?`, `TODO`, `_computed_`, `this_work`, `Photo-z`, `CIGALE`, `NED`, `TNS`.
- **Not in the `References` tab:** 34 codes used in the sheet, for example `Lee-Waddel2023` (vs `Lee-Waddell2023`), `Bruni25`, `Chen25`, `Xu22`, `CHME2025`, `Bassa2023`.

`link_sheet_refs.csv` gives a proposed repo code for each sheet code: exact, case only, the same first author with a year within 1 (**to confirm**), or a survey name. 92 codes have no proposal; most are papers that the repo does not use yet. The repo `all_refs.csv` has a DOI for each literature table, and the `References` tab has ADS links. A later step can match the two by bibcode, which would settle the year-drift cases.

## 9. Repo quantities that the sheet does not have

These could go on the website later; they are not compared now.
- **Host:** `Mtotal`, `age_mass`, `Z_stellar`, `lg_sSFR`, `AV_old`, `AV_young`, `AV_nebular`, `f_AGN`, `agn_tau`, `M_r`, `u-r`, `z_SED`, `SFR_photom`, `Mstar_spec`, Sérsic `n`. Also the Mannings2021 local quantities (`UV_SB`, `IR_SB`, `SSFR`, `SMStar`, `UVff`, `IRff`, `IRfe`).
- **FRB:** fluence, pulse width and scattering (9 FRBs), the systematic part of the error ellipse and its confidence level.

**The template.** `frb_host_template.json` has 39 fields; each one is in the `template_key` column of the CSV. While mapping, I noticed three things in `build_hosts_web.py`, for Xavier rather than for this review:
- It multiplies `eellipse` `a` and `b` by 3600, as if they were in degrees. They are already in arcsec.
- It writes the full FRB `refs` list into `DM_ref`, `RM_ref` and `z_ref`.
- `P_cc` and `offset_ref` have no source in the repo.

## 10. Questions for the reviewers

- **Q1. Sentinels.** What do `999`, `99`, `-99`, `-999`, `0` (in an error column) and `?` mean in each tab? In particular: is a magnitude error of `999` an upper limit, and is a magnitude of `999` "no data"?
- **Q2. Errors.** When the published error is asymmetric, what goes into the sheet's single error column: the mean, the larger of the two, or something else?
- **Q3. Magnitude convention.** Please confirm the `proposed_mag_convention` for each reference in `link_sheet_refs.csv`. Is the rule "as published in the paper" (sometimes observed, sometimes corrected)?
- **Q4. PanSTARRS.** Which Pan-STARRS magnitude does the sheet use (PSF, Kron, aperture), and from which data release? The repo uses the DR2 Kron magnitude.
- **Q5. FRB error ellipse.** Is `FRB_semi_a/b` statistical only, or statistical and systematic combined? At what confidence level? For the 38 hosts that disagree, we will check the convention before calling a discrepancy.
- **Q6. Morphology.** What are `Semiaxis_a`, `Semiaxis_b` and `PosAngle`, and how were they measured? Where do `Reff` and `Reff_kpc` come from?
- **Q7. Line fluxes.** Are they as published, or corrected for Galactic extinction? Is `OII_flux` the sum of the doublet?
- **Q8. Flags.** What do `Mstar_UL` and `SFR-Neb_UL-Dust` mean? Both are mostly `?`.
- **Q9. Photometric redshifts.** 20240304A has `is_zphot = y` in the sheet, but the repo stores its redshift as spectroscopic (`z_spec`), as for all 103 repo hosts. Which is right?
- **Q10. Offsets.** Is `calc_ang_off` meant to be the offset between the two sheet positions, and `Offset_ang` the published one? Which of the two should the website show?
