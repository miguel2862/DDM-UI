# Publication architecture and step-by-step observation

## Network architecture

Open **Network → Publication**. Drag the units to arrange the figure, then choose **Export PNG** or **SVG**. Exports include the complete network without editor controls, selection outlines or locks.

Diffuse fields indicate learning modulation; they do not add synaptic projections. Only explicitly configured connections are drawn. The notation distinguishes excitatory and inhibitory units and connections.

## Step-by-step observation

1. Configure the network, trials and phases.
2. In **Simulate**, enable **Record equations step by step**.
3. Choose a recording limit and run the simulation.
4. Use timestep/trial navigation and select a unit or connection to inspect its calculations.
5. Choose **Save trace** to retain the numerical record as JSON.

The detailed record belongs to the first simulated network. It contains the actual random threshold, asynchronous update order, activation operands, diffuse learning signals and weights before and after each update.

The network drawing shows the state at the end of the timestep. The equations show operands from the exact point when each unit was updated. The logistic plot is an analytical function with the observed input points; the weight plot is a recorded time series.

Navigation never reruns the simulation or draws new random values. Recording limits restrict the trace, not the simulation. Values unused by the active equation branch are shown as unavailable, not invented zeros.

Existing result files cannot supply operands that were never recorded. Enable recording before execution; editing the network, protocol or relevant parameters invalidates the previous trace. The large trace is not included in autosave, so export its JSON separately.
