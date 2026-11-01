#!/usr/bin/env bash
# For machines without cron or systemd user timers: runs nightly.sh every day at
# 03:00 in a loop. Start it in tmux or screen:
#   tmux new-session -d -s stew-cron /path/to/scheduler.sh
# It writes a heartbeat file so you can tell it is alive.
here=$(dirname "$(realpath "$0")")
STEW=${STEW:-stew}
ROOT=$("$STEW" config --get root) || { echo "stew: no catalogue found" >&2; exit 1; }
beat="$ROOT/.data-steward/logs/scheduler-heartbeat"
mkdir -p "$(dirname "$beat")"
while true; do
  next=$(date -d 'tomorrow 03:00' +%s)
  [[ $(date +%H%M) < 0300 ]] && next=$(date -d 'today 03:00' +%s)
  while (( $(date +%s) < next )); do date '+%F %T' > "$beat"; sleep 300; done
  STEW="$STEW" bash "$here/nightly.sh"
done
