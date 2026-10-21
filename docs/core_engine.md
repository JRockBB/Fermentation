# Consumer-resource engine (P1.3)

`ModelSpec` declares organism, substrate, and product names; a stoichiometric
matrix; uptake kinetics; and optional organism-specific response functions.
`build_rhs(spec)` produces the in-place OrdinaryDiffEq function. All states use
the order `[populations; substrates; products]`, and time is in hours.

Each organism has a growth flux `v_i = X_i * mu_max_i * modifier_i * uptake_i`.
The population rows of the stoichiometric matrix are the identity; the remaining
rows map growth to substrate consumption and product formation. Negative
substrate coefficients and positive product coefficients are required. For every
consumed substrate, uptake includes `S/(Ks+S)` (Monod) or `S/(Ks*X+S)` (Contois).
Multiplying these terms is an explicit co-limitation assumption. Cross-feeding
(a metabolite produced by one organism and consumed by another), alternative
substrates with additive uptake, maintenance, death, and variable yields need
separate model extensions.

References: [Monod (1949)](https://www.annualreviews.org/content/journals/10.1146/annurev.mi.03.100149.002103)
and [Contois (1959)](https://www.microbiologyresearch.org/content/journal/micro/10.1099/00221287-21-1-40).

Parameters are a NamedTuple with `mu_max` (one entry per organism, 1/h), `Ks`
(substrates × organisms), and any parameters needed by modifiers. A modifier
receives `(physical_state, parameters, time)` and must return a factor in [0,1].
For example, multiply a temperature response, a pH response, and acid inhibition
for the organism to which they apply. The engine does not infer species biology.

Monod Ks uses substrate concentration units. Contois Ks uses
substrate-concentration/population units. Yields must connect the stated population
basis to concentrations; CFU counts cannot be treated as mass without a sourced
conversion. Optional `conserved_weights` check that each stoichiometric column
conserves the declared mass or carbon quantity, rather than assuming equal weights.

`fermentation_problem(spec, physical_u0, tspan, parameters)` builds an ODEProblem
and rejects negative accepted states. RHS evaluation treats negative intermediate
trial values as zero without modifying the state. Set `log_populations=true` to
integrate natural logarithms of positive populations; substrates and products stay
in physical units. Use `physical_state` to decode results. Exact zero population
is supported in ordinary coordinates; logarithmic coordinates require positive
initial populations and do not invent a numerical population floor.

For ForwardDiff through parameters, promote physical initial values to the
derivative scalar type before constructing the problem, for example
`physical_u0 .* one(mu)` inside a function differentiated with respect to `mu`.
This lets the solver allocate derivative-compatible state storage. With log
populations, observations must decode the state before comparing physical units.

Tests compare numerical trajectories to the integrated Monod relation, verify
positivity and declared conservation, compare ordinary and logarithmic states,
check co-limitation and Contois values, and differentiate through an ODE solution.
