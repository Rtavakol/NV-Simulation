# Electromagnet Studio

Open electromagnet-studio.html from the simulator navigation. Current diamond, Hamiltonian constants, linewidth and bias direction transfer via sessionStorage. Applying a candidate transfers its field and model settings back and suppresses the auto-demo for that navigation.

## Calculation chain

1. Screen the full specified magnitude range with an electronic S=1 Hamiltonian over the specified direction or a Fibonacci hemisphere plus the current direction. Refine six shortlisted directions with the full S=1,I=1 model; probe larger fields if the separation minimum is not usable.
2. Calculate central-transition derivatives A (MHz/µT), temperature column kT (MHz/°C), and A+ (µT/MHz) via eigendecomposition of A-transpose A. Require rank three for the selected reconstruction subset. All eight profiles enter separation and capture checks.
3. Bound frequency excursions using A over the component box/ellipsoid, plus full-model field/temperature samples. Multiply the larger drift estimate by the reserve factor. Report the sampled reconstruction error, not a noise sensitivity.
4. Extract an isolated triplet's sinusoidal first-harmonic FM response. Linear interval ends at the first local-slope deviation beyond ±10% (editable); use its symmetric half-range. Distinguish the full width from the half-range. Optional measured half-range overrides capture checks.
5. Require both reserved drift inside each half-range and disjoint adjacent linear intervals. Independently require disjoint swept allocations with hyperfine span, tone offsets, FM deviation and linewidth guard. A separation-only candidate is explicitly not feasible.
6. Determine winding capacity from the rectangular envelope and insulated wire pitch, then either select the minimum turns satisfying wire/driver current and voltage limits or fill the envelope to minimize current. Sum the central on-axis Biot–Savart contribution of discrete circular turns. Compute current, copper resistance at the entered temperature, supply compliance, copper power, current density and axial endpoint variation.

## Limits

Search is restricted to 0.1–30 mT. Numerical direction and magnitude sampling does not certify a global optimum or coverage of every point of a continuous nonlinear domain. Firmware, instrument dynamics, strain, power broadening, coil thermal equilibrium and transverse homogeneity are not modeled. Temperature feed-forward assumes exact temperature and per-branch carrier control. All calculations are browser-side. The physics module is a model snapshot derived from index.html blob e354868701d753ce320d2720d3cd58bafafded24; regenerate/review it if the main Hamiltonian changes.

## Validation

Run: node tests/electromagnet-studio.test.cjs

Covers AWG diameter, single-loop analytical limit, field/current scaling, Ohm law and copper-temperature correction, pseudoinverse identity (analytic and actual model), slope-defined interval, default full-model separation-only design, zero-range feasible design, overlapping linear-interval rejection, and invalid insulated diameter. Browser checks cover settings transfer, search display, design download and bias application.

## Automatic coil sizing and wire limits

Successful field searches explicitly transfer the fitted magnitude into the disabled target field box, select the search-result source, and recalculate the winding. A manual-source option remains available. Minimum-turns mode selects turns under the effective current and voltage constraints; fill mode uses the full envelope. A partial outer layer is wound nearest the axial centre first.

The AWG reference current is diameter_in_mils² / 700, following PowerStream's conservative power-transmission wiring guideline (https://www.powerstream.com/Wire_Size.htm). This is not a universal wound-coil ampacity. Apply the entered derating factor (default 0.8), current-density cap, and supply-current limit. A verified manufacturer coil-current rating can replace the AWG reference. The code does not predict equilibrium temperature.

Wire resistance updates from the selected AWG area, actual selected turn-by-turn wire length plus leads, and entered copper temperature. The live wire panel shows reference/derated current and Ω/m at 20 °C and at the entered temperature. Coil input changes automatically recalculate the design. JSON export refreshes the coil calculation first.

Additional tests cover minimum-turn selection against the preceding rejected turn, AWG-dependent current caps and resistance per metre, and temperature-corrected total wire resistance.
