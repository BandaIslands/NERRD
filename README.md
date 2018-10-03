# NERRD
Near-Rotary-Resonance Relaxation Dispersion (Solid-state NMR)
Relaxation-dispersion (RD) NMR approaches have proven very successful for studying ms–ms motion. In solution-state NMR,
only fluctuations of the isotropic chemical shift are relevant when analyzing millisecond to microsecond conformational
dynamics. In contrast in solid-state NMR, the anisotropic interactions (dipolar couplings and chemical shift
anisotropy-CSA) do contribute. In MAS-ssNMR, the various time-dependent processes (sample rotation, rf irradiation and
conformational dynamics) may interfere and the fluctuation of the anisotropic relaxation induces relaxation.
Their influence on R1rho rates becomes more pronounced when the spin-lock field strength is close to the sample rotation
frequency, under the conditions of the so-called rotary resonance.
In order to analyze the relaxation profiles obtained in MAS-ssNMR, the evolution of a two-spin system that undergoes
exchange between two-states has been simulated using the GAMMA framework (1). The spin-evolution is computed first over
a rotor period and then over the spin lock period. As the experiments are recorded on powders, the computation is repeated
100 of crystal orientations. Simulation of the coherence decay is performed on a grid of values for the population of the
two states, the chemical shift difference, the angular fluctuation of the anisotropic interactions, the exchange rate
between the states, and the rf spin-lock.
Using a python program, the experimental profiles are then compared to the simulated profiles in the grid and evaluated
using a Chi2 value. The spectral and conformational parameters associated with the lowest values are considered.

(1) S. Smith, T. Levante, B. Meier, and R. Ernst, “Computer simulations in magnetic resonance. An object-oriented programming approach,” J. Magn. Reson., vol. 106, pp. 75–105, 1994
