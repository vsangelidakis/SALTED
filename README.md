# SALTED
 Sequential Particle Deposition

---

## What is SALTED?
SALTED is an efficient, open-source, event-driven algorithm for sequential ballistic deposition of complex-shaped rigid particles. Each of the particles is constructed as an agglomeration of hard spheres of variable radii, where the sizes and relative positions of the spheres may mutually overlap. In the sequential deposition process, the particles move along the steepest descent in a landscape formed by the boundaries and previously deposited particles, by performing steps of both rolling and translational motion. SALTED generalizes the Visscher–Bolsterli algorithm, which is frequently used for packing of spheres, to non-spherical particles. This code allows for the simulation of multi-million particle systems using desktop computers with reasonable time-runs, unlike typical soft-sphere algorithms, where the integration of Newton's equation of motion requires small timesteps and increased computational cost.


## Examples
SALTED can generate packings with complex particles of any shape, as demonstrated in Topic and Poeschel (2016).

<figure>
	<img src="figures/Topic_and_Poeschel_2016 - Figure_13.png" alt="drawing" width="1000"/>
	<figcaption> Fig. Packing of complex particles simulated with SALTED. Case studies for a) coins, b) ellipsoids, c) tetrapods, d) spirals and e) rods.
</figcaption>
</figure>


[//]: <> (## Architectural features)


## Acknowledging SALTED
Topic, N. and Pöschel, T., 2016. Steepest descent ballistic deposition of complex shaped particles. Journal of Computational Physics, 308, pp.421-437.

[//]: <> (<h4 align="center">2022 © Vasileios Angelidakis, Nikola Topic, Thorsten Poeschel. Institute for Multiscale Simulation, FAU, Germany </a></h4>)


