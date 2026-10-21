# Feeding and temperature (P1.4)

`dilute!(state, spec; fraction, substrate_addition)` mixes retained culture with
fresh substrate. It multiplies every physical state by the retained fraction,
then adds each substrate's contribution to the **final mixture concentration**.
For fresh-feed concentration `Sfresh`, that contribution is
`(1-fraction) .* Sfresh`. Products and organisms in the fresh feed are assumed
zero. Concentrations and populations must share a consistent dilution basis.

For a starter:flour:water mass ratio of 1:N:N, `feeding_fraction(N)` returns
`1/(1+2N)`. Using this fraction with volumetric concentrations requires a
consistent mass-to-volume convention. The helper alone does not convert flour
sugar content, hydration, CFU/g, or CFU/mL. Those conversions belong to the
food-specific model and need measured or published inputs.

`dilution_callback(spec, times; fraction, substrate_addition)` schedules identical
feeds using DiffEqCallbacks' `PresetTimeCallback`. Times must be unique and inside
the simulation interval. A feed at the initial time acts during initialization;
a feed at the final time is included in the final state. Both sides of interior
jumps are saved, so duplicate output times are intentional. The callback can be
reused for another solve. For different feed compositions or fractions, combine
separate callbacks with OrdinaryDiffEq's `CallbackSet` and use distinct times.

Logarithmic populations receive `log(fraction)` directly, so tiny populations do
not underflow during dilution. A zero retained fraction is valid in physical
coordinates and rejected for logarithmic populations.

```julia
temperature = TemperatureSchedule([0., 12., 24.], [28., 8., 28.])
stops = temperature_stops(temperature, problem.tspan)
feed = dilution_callback(spec, [12., 24.];
    fraction=feeding_fraction(1), substrate_addition=[20/3])
solution = solve(problem, Tsit5(); callback=feed,
    tstops=stops, d_discontinuities=stops)
```

An organism modifier must explicitly call `temperature(t)` to connect this input
to growth. `TemperatureSchedule` supports right-continuous `:step` (default) and
`:linear` interpolation; both hold their endpoint values outside the knot range.
Pass the interior knots to **both** `tstops` and `d_discontinuities`, so the solver
lands on each change and refreshes cached derivatives. For linear schedules this
also resolves changes in slope.

Tests check exact mixing, forced event times even with large solver steps,
endpoint feeds, callback reuse, convergence of repeated synthetic 1:1:1 feeds,
and the analytic exponential-growth integral across a temperature step. These
are numerical checks with synthetic constants, not calibrated starter results.

API references: [preset-time callbacks](https://docs.sciml.ai/DiffEqCallbacks/stable/timed_callbacks/)
and [solver time-step options](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/).
