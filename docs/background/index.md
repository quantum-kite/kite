KITE is a general purpose open-source tight-binding software for accurate real-space simulations of electronic structure and quantum transport properties of large-scale molecular and condensed systems with up to tens of billions of atomic orbitals ($N\sim 10^{10}$). In a nutshell, KITE takes real-space lattice models (tight-binding hamiltonians) of arbitrary complexity as an input that can be promptly defined by the user through its versatile Pybinding interface. Then, its memory-efficient and heavily-parallelised C++ code employs extremely accurate Chebyshev spectral expansions[^1] in order to extract *static electronic properties* (DoS, LDoS and spectral functions), *perturbative response functions* (linear and nonlinear conductivities) or even *dynamical effects* arising from the time-evolution of electronic wave-packets entirely in real-space. Since it is based upon real-space hamiltonians, KITE’s scope is not limited to periodic single-particle hamiltonians but, instead, its true power in unveiled through the study of more realistic models which may include disordered potentials, dilute impurities, structural defects, ad-atoms, mechanical strain and even external magnetic fields. Some illustrative examples may be found in the KITE’s presentation paper[^2]. See the [Tutorial][tutorial] section for a quick-start guide to the main features of KITE v1.1.

KITE's latest release (version 1.1) contains the following functionalities:

* Average density of states (DOS) and local DOS;
* $\mathbf{k}$-space Spectral Functions and ARPES Response;
* Linear DC conductivity (using the Kubo-Greenwood formula);
* First and Second-Order Optical (AC) Conductivities;
* Spin-Relaxation by Time-Evolution of Gaussian Wave-packets.

These calculations can now be applied to arbitrary two- and three-dimensional tight-binding hamiltonians that have:

* Generic multi-orbital local (on-site) and bond disorder;
* User-defined local potential profile and structural disorder;
* Different boundary conditions (periodic, open and twisted);
* Applied perpendicular magnetic field (limited use in 3D);

For more details about the current release refer to the documentation section.

# A Short Background Story

At the heart of the KITE software is an exact spectral expansion of broadened lattice Green's functions discovered independently by A. Ferreira (KITE's team) in collaboration with E. Mucciolo (U Central Florida)[^3] and by A. Braun and P. Schmitteckert (Karlsruhe Institute of Technology)[^2]. A large-RAM "single-shot" recursive algorithm developed by A. Ferreira (with technical assistance from M. D. Costa, National University of Singapore) enables the evaluation of zero-temperature response functions in systems with multi billions of orbitals $N\sim 10^{10}$.

In lattice models with small coordination number $[Z=O(1)]$, evaluations of response functions at fixed Fermi energy in large-memory nodes take only a few hours even when billions of spectral coefficients (Chebyshev moments) are retained. This gives access to accuracy and energy resolutions several orders of magnitude beyond previous approaches[^5]. To assess generic response functions at finite temperature/frequency, KITE implements the Green's function spectral approach to carry out a direct evaluation of the Kubo-Bastin formula as proposed by L. Covaci and T. G. Rappoport (KITE's team) in collaboration with J. H. García (ICN2)[^6].

To optimize multi-threading and speed up spectral expansions,  KITE provides the option to thread pre-defined partitions in real space (i.e., lattice domains) by means of a domain decomposition algorithm developed by J. Lopes (KITE's team).

[^1]: Kernel polynomial method, A. Weiße, G. Wellein, A. Alvermann and H. Fehske, [Rev. Mod. Phys. 78, 275 (2016)](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.78.275).

[^2]: KITE: high-performance accurate modelling of electronic structure and response functions of large molecules, disordered crystals and heterostructures, S. M. João, M. Anđelković, L. Covaci, T. G. Rappoport, João M. Viana Parente Lopes, and A. Ferreira, [R. Soc. open sci. 7, 191809 (2020)](https://royalsocietypublishing.org/doi/10.1098/rsos.191809).

[^3]: Critical delocalization of chiral zero energy modes in graphene, A. Ferreira and E. Mucciolo, [Phys. Rev. Lett. 115, 106601 (2015)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.115.106601).

[^4]: Numerical evaluation of Green's functions based on the Chebyshev expansion, A. Braun and P. Schmitteckert, Phys. [Rev. B 90, 165112 (2014)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.165112).

[^5]: Efficient multiscale lattice simulations of strained and disordered graphene, N. Leconte, A. Ferreira, and J. Jung. [Semiconductors and Semimetals 95, 35 (2016)](https://www.sciencedirect.com/science/article/abs/pii/S0080878416300047).

[^6]: Real-Space Calculation of the Conductivity Tensor for Disordered Topological Matter, J. H. García, L. Covaci, and T. G. Rappoport, [Phys. Rev. Lett. 114, 116602 (2015)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.114.116602).

[tutorial]: ../documentation/index.md
