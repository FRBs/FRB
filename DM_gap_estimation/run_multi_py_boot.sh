for j in {1..500}; do
  python DM_kde_runs_boot.py &
  sleep .25
done 2>/dev/null
