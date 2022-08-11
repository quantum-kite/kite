KITE evaluates generic electronic response functions and spectral properties of large-scale molecular and solid-state systems by means of extremely accurate spectral expansions of products of single-particle Green's functions[^1]. KITE's code uses as input lattice models (tight-binding matrices) of arbitrary complexity that can be imported from standard formats or defined directly via its versatile and user-friendly Pybinding interface.

At the heart of the KITE software is an exact spectral expansion of broadened lattice Green's functions discovered independently by A. Ferreira (KITE's team) in collaboration with E. Mucciolo (U Central Florida)[^2] and by A. Braun and P. Schmitteckert (Karlsruhe Institute of Technology)[^3]. A large-RAM "single-shot" recursive algorithm developed by by A. Ferreira (with technical assistance from M. D. Costa, National University of Singapore) enables the evaluation of zero-temperature response functions in systems with multi billions of orbitals $N\sim 10^{10}$.

In lattice models with small coordination number $[Z=O(1)]$, evaluations of response functions at fixed Fermi energy in large-memory nodes take only a few hours even when billions of spectral coefficients (Chebyshev moments) are retained. This gives access to accuracy and energy resolutions several orders of magnitude beyond previous approaches[^4]. To assess generic response functions at finite temperature/frequency, KITE implements the Green's function spectral approach to carry out a direct evaluation of the Kubo-Bastin formula as proposed by L. Covaci and T. G. Rappoport (KITE's team) in collaboration with J. H. García (ICN2)[^5].

The pre-release of KITE contains the following functionalities:

* Average density of states (DOS) and local DOS;
* Generic multi-orbital local (on-site) and bond disorder;
* Generic linear response functions for generic orbital observables;
* Linear and non-linear optical (AC) conductivity;

To optimize multi-threading and speed up spectral expansions,  KITE provides the option to thread pre-defined partitions in real space (i.e., lattice domains) by means of a domain decomposition algorithm developed by J. Lopes (KITE's team).

For more details about the current pre-release (including a to-do list) refer to the documentation section.

[^1]: Kernel polynomial method, A. Weiße, G. Wellein, A. Alvermann and H. Fehske, [Rev. Mod. Phys. 78, 275 (2016)](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.78.275).

[^2]: Critical delocalization of chiral zero energy modes in graphene, A. Ferreira and E. Mucciolo, [Phys. Rev. Lett. 115, 106601 (2015)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.115.106601).

[^3]: Numerical evaluation of Green's functions based on the Chebyshev expansion, A. Braun and P. Schmitteckert, Phys. [Rev. B 90, 165112 (2014)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.165112).

[^4]: Efficient multiscale lattice simulations of strained and disordered graphene, N. Leconte, A. Ferreira, and J. Jung. [Semiconductors and Semimetals 95, 35 (2016)](https://www.sciencedirect.com/science/article/abs/pii/S0080878416300047).

[^5]: Real-Space Calculation of the Conductivity Tensor for Disordered Topological Matter, J. H. García, L. Covaci, and T. G. Rappoport, [Phys. Rev. Lett. 114, 116602 (2015)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.114.116602).
