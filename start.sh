#!/bin/bash

set -u

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
API_PORT="${DDM_API_PORT:-8000}"
UI_PORT="${DDM_UI_PORT:-5173}"
API_URL="http://127.0.0.1:${API_PORT}/api/health"
API_PID=""
UI_PID=""

cleanup() {
  if [[ -n "$API_PID" ]]; then kill "$API_PID" 2>/dev/null || true; fi
  if [[ -n "$UI_PID" ]]; then kill "$UI_PID" 2>/dev/null || true; fi
}
trap cleanup EXIT INT TERM

if curl --silent --max-time 1 "$API_URL" >/dev/null 2>&1; then
  echo "[ERROR] Port ${API_PORT} already serves an API. Stop it or set DDM_API_PORT."
  exit 1
fi

echo "[API] Starting DDM on 127.0.0.1:${API_PORT}..."
(
  cd "$ROOT_DIR/api" || exit 1
  Rscript -e "plumber::plumb('plumber.R')\$run(host='127.0.0.1', port=${API_PORT})"
) &
API_PID=$!

API_READY=false
for _ in {1..40}; do
  if ! kill -0 "$API_PID" 2>/dev/null; then
    echo "[ERROR] The R API exited during startup."
    exit 1
  fi
  HEALTH="$(curl --silent --max-time 1 "$API_URL" 2>/dev/null || true)"
  if [[ "$HEALTH" == *'"DDM"'* && "$HEALTH" == *'"dtd"'* ]]; then
    API_READY=true
    break
  fi
  sleep 0.25
done

if [[ "$API_READY" != true ]]; then
  echo "[ERROR] Health check did not confirm the DDM engine."
  exit 1
fi

echo "[UI] Starting Vite on 127.0.0.1:${UI_PORT}..."
(
  cd "$ROOT_DIR/frontend" || exit 1
  VITE_API_TARGET="http://127.0.0.1:${API_PORT}" npm run dev -- --host 127.0.0.1 --port "$UI_PORT"
) &
UI_PID=$!

echo "DDM-UI is ready: http://127.0.0.1:${UI_PORT}"
echo "Press Ctrl+C to stop both processes."

wait "$UI_PID"
