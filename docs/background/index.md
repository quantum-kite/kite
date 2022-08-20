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

The seeds for KITE’s project were laid in 2014, when an exact spectral expansion of the broadened lattice Green's function was discovered by **Aires Ferreira** in collaboration with Eduardo R. Mucciolo (University of Central Florida)[^3] and, independently, by A. Braun and P. Schmitteckert (Karlsruhe Institute of Technology)[^4]. **Aires Ferreira** then developed a large-RAM "single-shot" recursive algorithm (with the technical assistance of M. D. Costa, National University of Singapore) that enabled the efficient evaluation of zero-temperature DC conductivity of huge tight-binding systems (more than 10 billion orbitals) entirely in real-space. At that time, this method proved essential to numerically show that the DC conductivity of graphene indeed attains finite value in the presence of dilute vacancies[^3].

In the following years, the usefulness of the above method in generic disordered lattices[5] with small coordination number $[Z=O(1)]$ became ever more clear by the development of similar techniques to  evaluate generic response functions, most notably the linear conductivity tensor  at finite temperature/frequency (proposed by **L. Covaci** and **T. G. Rappoport** in collaboration with José H. García (ICN2)[^6] [^7]) and the non-linear optical response (proposed by **S. M. João** and **J. M. Viana Parente Lopes**[^8]). It was the conjunction of all these proposals that put forward the joint venture that led to the first release of KITE in 2018[^2].

From its beginnings, KITE was sought to handle real-space lattice models of arbitrary complexity and realistically large sizes. Thereby, its architecture allies a versatile and user-friendly $\texttt{python}$ interface, with a highly efficient and CPU-scalable $\texttt{C++}$ code (developed by **J. M. Viana Parente Lopes**) that handles the heavy spectral computations. The interface is based on Pybinding’s syntax[^9] that allows the user to input an arbitrary (2D or 3D) lattice model decorated with a myriad of non-periodic perturbations, such as on-site disorder, personalized structural defects, strain, etc... This model hamiltonian is then passed to the $\texttt{C++}$ high-performance code ($\texttt{KITEx}$) that implements a *matrix-free* Chebyshev iteration which combines a domain-decomposition of the lattice with a *"tile-by-tile"* matrix-vector multiplication strategy in order to minimize memory-transfer overheads and thus boost parallelization and calculational efficiency[^2]. Such an approach allows the spectral calculations to be carried out for unheard-of huge system sizes, particularly by using modern large-RAM multi-core computers. Finally, a convenient *post-processing tool* ($\texttt{KITE-tools}$), developed by **S. M. João**, was also included in the package thereby turning KITE into a ready-to-use tool for practical applications. Teste

[^1]: Kernel polynomial method, A. Weiße, G. Wellein, A. Alvermann and H. Fehske, [Rev. Mod. Phys. 78, 275 (2016)](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.78.275).

[^2]: **KITE:** high-performance accurate modelling of electronic structure and response functions of large molecules, disordered crystals and heterostructures, S. M. João, M. Anđelković, L. Covaci, T. G. Rappoport, João M. Viana Parente Lopes, and A. Ferreira, [R. Soc. open sci. 7, 191809 (2020)](https://royalsocietypublishing.org/doi/10.1098/rsos.191809).

[^3]: Critical delocalization of chiral zero energy modes in graphene, A. Ferreira and E. Mucciolo, [Phys. Rev. Lett. 115, 106601 (2015)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.115.106601).

[^4]: Numerical evaluation of Green's functions based on the Chebyshev expansion, A. Braun and P. Schmitteckert, Phys. [Rev. B 90, 165112 (2014)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.165112).

[^5]: Efficient multiscale lattice simulations of strained and disordered graphene, N. Leconte, A. Ferreira, and J. Jung. [Semiconductors and Semimetals 95, 35 (2016)](https://www.sciencedirect.com/science/article/abs/pii/S0080878416300047).

[^6]: Real-Space Calculation of the Conductivity Tensor for Disordered Topological Matter, J. H. García, L. Covaci, and T. G. Rappoport, [Phys. Rev. Lett. 114, 116602 (2015)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.114.116602).

[^7]: Numerical calculation of the Casimir-Polder interaction between a graphene sheet with vacancies and an atom, T. P. Cysne, T. G. Rappoport, Aires Ferreira, J. M. Viana Parente Lopes, and N. M. R. Peres, [Phys. Rev. B 94, 235405 (2016)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.94.235405)

[^8]: Basis-independent spectral methods for non-linear optical response in arbitrary tight-binding models, S. M. João and J. M. Viana Parente Lopes, [J. Phys.: Condens. Mat. 32 (12), 125901 (2019)](https://iopscience.iop.org/article/10.1088/1361-648X/ab59ec/meta)

[^9]: D. Moldovan, M. Anđelković, F. M. Peeters, **Pybinding V0.9.4:** a python package for tight-binding calculations [Zenodo (2017)](doi:10.5281/zenodo.826942)

[tutorial]: ../documentation/index.md
