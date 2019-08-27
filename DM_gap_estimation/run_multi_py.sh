for j in {1..500}; do
  python save_kde_redraws.py &
  sleep .25
done 2>/dev/null
