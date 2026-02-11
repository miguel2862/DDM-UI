#!/bin/bash
echo "Starting DDM-UI..."
echo ""

# Start R Plumber API in background
echo "[API] Starting R Plumber API on port 8000..."
cd "$(dirname "$0")/api"
Rscript -e "plumber::plumb('plumber.R')\$run(host='0.0.0.0', port=8000)" &
API_PID=$!
echo "[API] PID: $API_PID"

# Wait for API to be ready
echo "[API] Waiting for API to start..."
sleep 3

# Start React frontend
echo "[UI] Starting React frontend on port 5173..."
cd "$(dirname "$0")/frontend"
npm run dev &
UI_PID=$!
echo "[UI] PID: $UI_PID"

echo ""
echo "========================================="
echo "  DDM-UI is running!"
echo "  Frontend: http://localhost:5173"
echo "  API:      http://localhost:8000"
echo "========================================="
echo ""
echo "Press Ctrl+C to stop both servers."

# Wait for either process and clean up on exit
trap "kill $API_PID $UI_PID 2>/dev/null; exit" INT TERM
wait
