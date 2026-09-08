# Compatibility entry point.
#
# Existing scripts that source api/simulation.R continue to receive the stable
# DDM implementation, without loading any alternative engine.

.simulation_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.simulation_dir <- if (is.null(.simulation_file)) getwd() else dirname(normalizePath(.simulation_file, mustWork = FALSE))
source(file.path(.simulation_dir, "simulation.R.dtd-backup"))
rm(.simulation_file, .simulation_dir)
