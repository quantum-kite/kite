The tight-binding (TB) method (and the closely related linear combination of atomic orbitals (LCAOs) method in quantum chemistry)
is a computationally fast and robust approach to handle large-scale molecular and condensed matter systems.
In the TB approximation, electrons are assumed to be strongly bound to the nuclei.
The one-particle wavefunctions $\{\psi_{\alpha}(\mathbf{x})\}$ are approximated by linear combinations
of Slater-Koster-type states (i.e., LCAOs) for isolated atoms, i.e.,

$$
    \psi_{\alpha}(\mathbf{x}) = \frac{1}{\sqrt{N}}\sum\limits_{i=1}^N a_\alpha^{(i)} \phi_\text{SK}(\mathbf{x}-\mathbf{x}_i),
$$

where $i=\{1, \ldots, N\}$ runs over all sites and orbitals.
The one-particle states $| \psi_\alpha \rangle$ are eigenvectors of the Hamiltonian matrix (the TB Hamiltonian),
$\hat{H} = \sum_{i,j} t_{i,j} |i\rangle \langle j |$.
The TB matrix elements — encoding on-site energies $(i = j)$ and hopping integrals between different atomic orbitals
$(i \neq j)$ — can be estimated for example by means of the Slater-Koster approach or by matching the spectrum
to that obtained by first-principles calculations in a suitable reference system[^1][^2][^3][^4].

Suitably parameterized TB models provide an accurate description of molecular orbitals in molecules and Bloch wavefunctions in many solids.
The complexity of TB models only grows linearly with the number of atomic orbitals,
providing a basis for large-scale calculations of a plethora of equilibrium and non-equilibrium physical properties,
including optical absorption spectra, simulations of amorphous solids, and wave-packet propagation.
Disorder, interfaces, and defects can be conveniently added to a TB model by modifying on-site energies and hopping integrals,
and adding auxiliary sites. Such a multi-scale approach has proven very successful in describing impurity scattering[^5][^6],
moiré patterns[^7], complex scattering potentials induced by ad-atoms[^8], optical conductivity of disordered 2D materials[^9], and geometrical properties, vibrational frequencies and non-covalent interactions of large molecular systems[^10].


[^1]: Simplified LCAO method for the periodic potential problem, J. C. Slater and G. F. Koster, [Phys. Rev. 94, 1498 (1954)](https://journals.aps.org/pr/abstract/10.1103/PhysRev.94.1498)

[^2]: Elementary prediction of linear combination of atomic orbitals matrix elements, S. Froyen and W.A. Harrison, [Phys. Rev. B 20, 2420 (1979)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.20.2420)

[^3]: Tight-binding modelling of materials, C. M. Goringe, D. R. Bowler, and E. Hernández, [Rep. Prog. Phys. 60, 1447 (1997)](https://iopscience.iop.org/article/10.1088/0034-4885/60/12/001/pdf)

[^4]: The Slater–Koster tight-binding method: a computationally efficient and accurate approach, D. A. Papaconstantopoulos and M. J. Mehl, [Journal of Physics: Condensed Matter 15, R413 (2003)](https://iopscience.iop.org/article/10.1088/0953-8984/15/10/201)

[^5]: Resonant scattering by realistic impurities in graphene, T. O. Wehling et al. [Phys. Rev. Lett. 105, 056802 (2010)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.105.056802)

[^6]: Unified description of the dc conductivity of monolayer and bilayer graphene at finite densities based on resonant scatterers, A. Ferreira et al., [Phys. Rev. B 83, 165402 (2011)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.83.165402)

[^7]: Ab initio theory of moiré superlattice bands in layered two-dimensional materials, J. Jung, A. Raoux, Z. Qiao, and A. H. MacDonald, [Phys. Rev. B 89, 205414 (2014)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.89.205414)

[^8]: Impact of complex adatom-induced interactions on quantum spin Hall phases. F. J. dos Santos et al., [Phys. Rev. B 98, 081407(R) (2018)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.98.081407)

[^9]: Numerical calculation of the Casimir-Polder interaction between a graphene sheet with vacancies and an atom. T. Cysne et al., [Phys. Rev. B 94, 235405 (2016)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.94.235405)

[^10]: A robust and accurate tight-binding quantum chemical method for structures, vibrational frequencies, and noncovalent interactions of large molecular systems parametrized for all spd-block elements (Z = 1−86), S. Grimme , C. Bannwarth, and P. Shushkov, [J. Chem. Theory Comput., 13 , 1989 (2017)](https://pubs.acs.org/doi/full/10.1021/acs.jctc.7b00118)
