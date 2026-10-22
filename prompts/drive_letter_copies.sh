#!/bin/bash
# Copy Galaxy_DB products stored under the pre-TNS-letter host names
# (HG20190523_...) to the lettered names that build_hosts looks for
# (HG20190523A_...), on Google Drive (rclone remote GDrive:).
# Server-side copies; the original files are kept.  See rebuild_hosts.md,
# prompt 3 and Q2.
#
#   bash prompts/drive_letter_copies.sh --dry-run   # show what would be copied
#   bash prompts/drive_letter_copies.sh             # copy
#
# Afterwards, sync the local Galaxy_DB from the Drive, e.g.
#   rclone copy GDrive:Astronomy/Research/FRB/Galaxies/Galaxy_DB /mnt/tank/Astronomy/Galaxy_DB --dry-run

DRIVE="GDrive:Astronomy/Research/FRB/Galaxies/Galaxy_DB"
DRY="$1"; [ "$1" == "--remove-originals" ] && DRY="--dry-run"

copy() {   # copy <dir> <old name> <new name>
    rclone copyto $DRY "$DRIVE/$1/$2" "$DRIVE/$1/$3" -v
}

# 20180301A (Realfast/Bhandari2021)
for f in CIGALE.fits CIGALE_SFH.fits CIGALE_model.fits DEIMOS_ppxf.ecsv DEIMOS_ppxf.fits GMOS-S_r_galfit.fits; do
    copy Realfast/Bhandari2021 HG20180301_$f HG20180301A_$f
done
# 20190523A (Realfast/Bhandari2021 and CRAFT/Heintz2020)
for f in CIGALE.fits CIGALE_SFH.fits CIGALE_model.fits CIGALE_photo.dat; do
    copy Realfast/Bhandari2021 HG20190523_$f HG20190523A_$f
done
copy CRAFT/Heintz2020 HG20190523_LRIS_r_galfit.fits HG20190523A_LRIS_r_galfit.fits
# 20191228A (Realfast/Bhandari2021)
for f in CIGALE.fits CIGALE_SFH.fits CIGALE_model.fits LRISr_ppxf.ecsv LRISr_ppxf.fits; do
    copy Realfast/Bhandari2021 HG20191228_$f HG20191228A_$f
done
# 20201123A (MeerTRAP/rajwade2022)
for f in Goodman_ppxf.ecsv Goodman_ppxf.fits; do
    copy MeerTRAP/rajwade2022 HG20201123_$f HG20201123A_$f
done
# Not copied: HG20190611_* and HG20190714_* (Realfast/Bhandari2021). Those hosts
# already get CIGALE values from other files, and a lettered copy here could
# change which CIGALE fit the build picks (search_for_file takes the last match).

# ---------------------------------------------------------------------------
# Step 2 (2026-10-21, at Xavier's request): remove the letterless originals,
# after checking that the lettered copy exists with the same size.  Google
# Drive moves deleted files to the trash (rclone's default --drive-use-trash).
#   bash prompts/drive_letter_copies.sh --remove-originals
if [ "$1" == "--remove-originals" ]; then
    for spec in \
        "Realfast/Bhandari2021 HG20180301 CIGALE.fits CIGALE_SFH.fits CIGALE_model.fits DEIMOS_ppxf.ecsv DEIMOS_ppxf.fits GMOS-S_r_galfit.fits" \
        "Realfast/Bhandari2021 HG20190523 CIGALE.fits CIGALE_SFH.fits CIGALE_model.fits CIGALE_photo.dat" \
        "CRAFT/Heintz2020 HG20190523 LRIS_r_galfit.fits" \
        "Realfast/Bhandari2021 HG20191228 CIGALE.fits CIGALE_SFH.fits CIGALE_model.fits LRISr_ppxf.ecsv LRISr_ppxf.fits" \
        "MeerTRAP/rajwade2022 HG20201123 Goodman_ppxf.ecsv Goodman_ppxf.fits"; do
        set -- $spec; dir=$1; base=$2; shift 2
        for f in "$@"; do
            old="$DRIVE/$dir/${base}_$f"; new="$DRIVE/$dir/${base}A_$f"
            so=$(rclone size --json "$old" 2>/dev/null | python3 -c "import json,sys;print(json.load(sys.stdin)['bytes'])")
            sn=$(rclone size --json "$new" 2>/dev/null | python3 -c "import json,sys;print(json.load(sys.stdin)['bytes'])")
            if [ -n "$sn" ] && [ "$so" == "$sn" ] && [ "$sn" != "0" ]; then
                rclone deletefile "$old" -v && echo "REMOVED $dir/${base}_$f"
            else
                echo "KEPT $dir/${base}_$f (copy missing or size mismatch: $so vs $sn)"
            fi
        done
    done
fi
