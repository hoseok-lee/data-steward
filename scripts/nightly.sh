#!/usr/bin/env bash
# Nightly upkeep of a data-steward catalogue: back it up, then check every dataset.
# Run it from cron or a systemd timer (see README, "Nightly checks"), e.g.
#   0 3 * * *  /path/to/nightly.sh
# Uses the catalogue stew finds by itself; set STEW_ROOT (or STEW) to override.
set -uo pipefail
STEW=${STEW:-stew}
ROOT=$("$STEW" config --get root) || { echo "stew: no catalogue found" >&2; exit 1; }
LOGS="$ROOT/.data-steward/logs"
mkdir -p "$LOGS"
day=$(date +%F)
{
  echo "== $(date '+%F %T') backup"
  "$STEW" backup --keep 30
  echo "== $(date '+%F %T') check"
  "$STEW" check --report "$LOGS/check-latest.txt"; rc=$?
  echo "== $(date '+%F %T') done (check exit $rc)"
} >> "$LOGS/nightly-$day.log" 2>&1
ls -1t "$LOGS"/nightly-*.log | tail -n +61 | xargs -r rm --
