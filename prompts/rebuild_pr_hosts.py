""" Rebuild the host JSONs touched by PR #240 with the current code and
report what changed (add_frbs.md, prompts PR-2 to PR-11).

Run from the repo root, in the `astro` env, on a machine with $FRB_GDB
(and ideally $NEDLVS) set:

    python prompts/rebuild_pr_hosts.py --dry-run          # check the setup
    python prompts/rebuild_pr_hosts.py                    # all 55 PR hosts
    python prompts/rebuild_pr_hosts.py --hosts 20190711A,20200223B
    python prompts/rebuild_pr_hosts.py --retry 2          # re-run hosts that lost a band
    python prompts/rebuild_pr_hosts.py --hosts 20180916B --keep-morphology 20180916B
                                    # build without Galfit, keep the committed morphology

For every host the committed JSON (git HEAD) is compared with the rebuilt
one.  A host that *lost* photometry bands not listed in its Bad_photom
column is flagged as a possible flaky-survey failure and, with --retry,
rebuilt again.  Lost neb_lines, derived or morphology keys are always
flagged.  Nothing is committed; inspect `git diff` afterwards.
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np

# The 55 host JSONs added, modified or renamed by PR #240
# (git diff --name-only --diff-filter=AMR origin/main...HEAD -- 'frb/data/Galaxies/*host.json')
PR_HOSTS = """20171020A 20180301A 20180814A 20181030A 20190110C 20190208A 20190303A
20190417A 20190711A 20191106C 20200223B 20210807D 20211203C 20220204A 20220207C
20220208A 20220307B 20220310F 20220319D 20220330D 20220418A 20220506D 20220509G
20220726A 20220825A 20220912A 20220914A 20220920A 20221012A 20221027A 20221029A
20221101B 20221113A 20221116A 20221219A 20230124A 20230216A 20230222B 20230307A
20230501A 20230626A 20230628A 20230712A 20230926A 20231011A 20231020A 20231120A
20231123A 20231123B 20231201A 20231223C 20231229A 20231230A 20240114A 20240209A""".split()

GAL_PATH = 'frb/data/Galaxies/'
SECTIONS = ['photom', 'neb_lines', 'derived', 'morphology', 'offsets', 'redshift']


def json_path(host):
    return os.path.join(GAL_PATH, host, f'FRB{host}_host.json')


def committed_json(host):
    """The JSON as committed at HEAD (None if it is new)."""
    try:
        txt = subprocess.run(['git', 'show', f'HEAD:{json_path(host)}'],
                             capture_output=True, text=True, check=True).stdout
    except subprocess.CalledProcessError:
        return None
    return json.loads(txt)


def bands_of(photom):
    return sorted(k for k in photom if not k.endswith(('_err', '_ref', '_flux', '_flux_err')))


def compare(old, new, bad_photom):
    """Compare two host dicts.  Returns (flags, details)."""
    flags, details = [], []
    if old is None:
        return ['new'], ['no committed JSON to compare with']
    # Photometry bands
    ob, nb = set(bands_of(old.get('photom', {}))), set(bands_of(new.get('photom', {})))
    lost = sorted(ob - nb)
    unexpected = [b for b in lost if b not in bad_photom]
    if lost:
        details.append(f"lost bands: {', '.join(lost)}")
    if unexpected:
        flags.append('lost_bands')
    if nb - ob:
        details.append(f"new bands: {', '.join(sorted(nb - ob))}")
    # Other sections: lost keys
    for sec in ['neb_lines', 'derived', 'morphology']:
        ok, nk = set(old.get(sec, {})), set(new.get(sec, {}))
        if ok - nk:
            flags.append(f'lost_{sec}')
            details.append(f"lost {sec}: {', '.join(sorted(ok - nk))}")
        if nk - ok:
            details.append(f"new {sec}: {', '.join(sorted(nk - ok))}")
    # Changed values
    changed = []
    for sec in SECTIONS:
        for k, v in new.get(sec, {}).items():
            if k not in old.get(sec, {}):
                continue
            o = old[sec][k]
            if isinstance(v, (int, float)) and isinstance(o, (int, float)):
                if not np.isclose(v, o, rtol=1e-5, atol=1e-6, equal_nan=True):
                    changed.append(f'{sec}.{k}: {o:g} -> {v:g}')
            elif v != o:
                changed.append(f'{sec}.{k}: {o!r} -> {v!r}')
    if changed:
        details.append(f"changed ({len(changed)}): " + '; '.join(changed[:12])
                       + (' ...' if len(changed) > 12 else ''))
    return flags, details


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--hosts', help='Comma-separated host names (default: the 55 PR hosts)')
    ap.add_argument('--retry', type=int, default=0,
                    help='Rebuild a host up to N more times if it lost photometry bands')
    ap.add_argument('--keep-nedlvs', action='store_true',
                    help='Query NED-LVS even if $NEDLVS is unset (the build will crash)')
    ap.add_argument('--keep-morphology', default='',
                    help='Comma-separated hosts to build without Galfit (Galfit file missing); '
                         'their committed morphology block is copied into the new JSON')
    ap.add_argument('--dry-run', action='store_true', help='Check the setup and list the hosts')
    ap.add_argument('--report', default='rebuild_report.md', help='Markdown report file')
    args = ap.parse_args()

    hosts = args.hosts.split(',') if args.hosts else PR_HOSTS
    keep_morph = [h for h in args.keep_morphology.split(',') if h]

    # ---- Setup checks
    import frb
    repo = os.path.realpath('.')
    frb_file = os.path.realpath(frb.__file__)
    if not frb_file.startswith(repo):
        print(f"WARNING: frb is imported from {frb_file}, not from this repo. "
              "Run `pip install -e .` in this env.")
    gdb = os.getenv('FRB_GDB')
    if not gdb or not os.path.isdir(gdb):
        print(f"WARNING: $FRB_GDB = {gdb!r} is not a directory. pPXF/CIGALE/Galfit "
              "products will not be found (pPXF/CIGALE silently; Galfit raises).")
    from frb.surveys import survey_utils
    if not os.getenv('NEDLVS') and not args.keep_nedlvs:
        if 'NEDLVS' in survey_utils.optical_surveys:
            survey_utils.optical_surveys.remove('NEDLVS')
        print("$NEDLVS is not set: NED-LVS removed from the survey list (adds no photometry).")
    from frb.builds import build_hosts
    from frb.galaxies import hosts as frbhosts
    host_tbl = frbhosts.load_host_tbl()
    missing = [h for h in hosts if h not in host_tbl.FRB.values]
    if missing:
        print(f"ERROR: not in public_hosts.csv: {missing}")
        sys.exit(1)
    print(f"{len(hosts)} hosts to rebuild; frb from {frb_file}")
    if args.dry_run:
        print(' '.join(hosts))
        return

    # ---- Build and compare
    rows = []
    for host in hosts:
        row = host_tbl[host_tbl.FRB == host].iloc[0]
        if host in keep_morph:
            row = row.copy()
            row['Galfit_filter'] = np.nan
        bad_photom = [s.strip() for s in str(row.Bad_photom).split(',')] \
            if isinstance(row.Bad_photom, str) else []
        old = committed_json(host)
        attempt, flags, details, error = 0, [], [], None
        while True:
            attempt += 1
            print(f"\n===== {host} (attempt {attempt}) =====")
            try:
                build_hosts.run(row)
                error = None
            except Exception as e:  # keep going with the other hosts
                error = f'{type(e).__name__}: {e}'
                print(f"FAILED: {error}")
            if error is None:
                with open(json_path(host), encoding='utf-8') as f:
                    new = json.load(f)
                if host in keep_morph and old is not None and 'morphology' in old:
                    from frb import utils as frb_utils
                    new['morphology'] = old['morphology']
                    frb_utils.savejson(json_path(host), new, easy_to_read=True, overwrite=True)
                    print(f"Kept the committed morphology for {host}")
                flags, details = compare(old, new, bad_photom)
            if error is None and 'lost_bands' not in flags:
                break
            if attempt > args.retry:
                break
            print("Retrying after a short pause ...")
            time.sleep(20)
        rows.append((host, attempt, error, flags, details))

    # ---- Report
    lines = ['# Rebuild report', '',
             '| Host | Tries | Status | Flags | Details |', '|---|---|---|---|---|']
    for host, attempt, error, flags, details in rows:
        status = 'FAILED' if error else 'ok'
        lines.append(f"| {host} | {attempt} | {status} | {', '.join(flags)} | "
                     f"{error or ' / '.join(details)} |")
    failed = [r[0] for r in rows if r[2]]
    flagged = [r[0] for r in rows if r[3] and not r[2]]
    lines += ['', f"Failed: {', '.join(failed) or 'none'}",
              f"Flagged (check before committing): {', '.join(flagged) or 'none'}"]
    with open(args.report, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')
    print('\n'.join(lines[-2:]))
    print(f"Report written to {args.report}")


if __name__ == '__main__':
    main()
