!!! Info
    
    Here, we briefly introduce the **main concepts** underlying all spectral methods currently implemented in KITE.
    The focus is placed on the role of the simulation parameters that can be adjusted by the end user in order to suit their specific purposes.

# The Ground Rules for Spectral Methods

The central object to any calculation done in KITE is the lattice single-particle Hamiltonian (SPH) — $H$ — which is
always a sparse $D\!\times\!D$ hermitian matrix ($D$ being the total number of orbitals in the lattice).
This fully embodies the simulated system, from its underlying bravais lattice structure and local orbital basis,
to all the non-periodic terms that realize disordered potentials, specific boundary conditions,
complex structural defects and external magnetic fields.

Any finite-dimensional SPH has a bounded real-valued spectrum which must be shifted and rescaled to suitably fit within $[-1,1]$, the convergence interval of the method. This conversion is performed internally by KITE, which transforms $H\to\mathcal{H}=(H\!-\!\varepsilon_{0})/\delta\varepsilon$ and rescales all energy variables by $\delta\varepsilon$.
This step requires an early (over)estimation of $H$'s spectral bandwidth, which can be set manually (using [`#!python spectrum_range`][configuration-spectrum_range] = $[\,\varepsilon_{0}\!-\!\delta\varepsilon/2$,$\,\varepsilon_{0}\!+\!\delta\varepsilon/2]$ in [`#!python kite.Configuration()`][configuration]), or be automatically done by KITE upon generation of the hdf5 configuration file. At this stage, the user must also specify:

1. The dimensions of the simulated lattice, `#!python length=[lx,ly,(lz)]`;
2. The number of subdivisions for parallelization of the matrix-vector operation, `#!python divisions=[nx,ny,(nz)]`,
   where $n_{x}n_{y}n_{z}$ (or $n_{x}n_{y}$ for 2D models) is the available CPU-cores;
3. The type of boundary conditions;
4. The type of data to be handled in internal arithmetic operations (for further information see [Settings](settings.md)).

Generally, a target function $\mathcal{Q}$ evaluated by [`#!bash KITEx`](../api/kitex.md) fits one of the following forms:

$$
    \mathcal{Q}\left(\left\{ \lambda_{i}^{j}\right\} \right)=\begin{cases}
\text{Tr}\left[F_{1}\left(\left\{ \lambda_{i}^{1}\right\};\mathcal{H}\right)\mathcal{O}_{1}F_{2}\left(\left\{ \lambda_{i}^{2}\right\};\mathcal{H}\right)\mathcal{O}_{2}\cdots\mathcal{O}_{N}F_{N}\left(\left\{ \lambda_{i}^{N}\right\};\mathcal{H}\right)\right]\\
\left\langle \Psi\right|F_{1}\left(\left\{ \lambda_{i}^{1}\right\};\mathcal{H}\right)\mathcal{O}_{1}F_{2}\left(\left\{ \lambda_{i}^{2}\right\};\mathcal{H}\right)\mathcal{O}_{2}\cdots\mathcal{O}_{N}F_{N}\left(\left\{ \lambda_{i}^{N}\right\};\mathcal{H}\right)\left|\Psi\right\rangle 
\end{cases}\,\,\, ,
$$

where $\mathcal{O}_{j}$ are sparse lattice operators (e.g., identities, velocity operators or spin operators) and $F_{j}$ are functions of $\mathcal{H}$, as well as other scalar parameters 
$\left\{ \lambda_{i}^{j}\right\} _{i=i_{min},...,i_{max}}^{j=1,...,N}$, such as the temperature or frequency. Furthermore, the $\left|\Psi\right\rangle$ in the second line of the above equation is a specific state/basis vector that depends on the observable that is being computed. Hence, there are two distinct categories of observables that can be computed with  [`#!bash KITEx`](../api/kitex.md), as shown below.

## Traces Over The Full Hilbert Space

:  Complete traces of operators are needed to investigate global observables, such as the Density of States (DoS) and DC-conductivity, as well as the $1^{\text{st}}$- and $2^{\text{nd}}$-order optical conductivity (see João *et al.*[^1] for further details). For all these cases, the trace is evaluated stochastically as an average of expectation values for $R$ normalized random vectors[^2], i.e.,

$$
  \text{Tr}\left[\cdots\right]\approx\frac{1}{R}\sum_{r=1}^{R}\left\langle\xi_{r}\right|\cdots\left|\xi_{r}\right\rangle .
$$

:  Within the user interface, the number of independent random vectors is specified by the parameter `#!python num_random`, which must be large enough to ensure a well-estimated trace. The associated error scales as $1/\sqrt{R\,D}$, and thus requires very few random vectors if the simulated system is very large[^2]. On top of this averaging, if $\mathcal{H}$ has some random component (e.g., by hosting disorder or featuring randomly twisted boundaries), it is often the case that the results need to be averaged over an ensemble of random Hamiltonians. Such averaging is also done inside  [`#!bash KITEx`](../api/kitex.md) and the number of random configurations is specified by user with the parameter `#!python num_disorder`.

## Diagonal Matrix Elements

: This class of target functions includes local observables such as the local density of states (LDoS) and the $\mathbf{k}$-space spectral function (for ARPES's response), as well as the time-evolution of Gaussian wave-packets. Note that `#!python num_random` is no longer a relevant parameter for these target functions.

: In both classes, the core of the method is to expand the functions $F_{j}$ as a truncated Chebyshev series of $\mathcal{H}$, which allows one to write 

$$
\mathcal{Q}\left(\left\{ \lambda_{i}^{j}\right\} \right)\!=\!\begin{cases}
\sum_{n_{1}=0}^{M-1}\!\cdots\!\sum_{n_{N}=0}^{M-1}G_{1}\!\left(\left\{ \lambda_{i}^{1}\right\};n_1\right)\cdots G_{N}\!\left(\left\{ \lambda_{i}^{N}\right\};n_N\right)\,\,\text{Tr}\left[T_{n_{1}}\!\left(\mathcal{H}\right)\mathcal{O}_{1}\cdots\mathcal{O}_{N}T_{n_{N}}\!\left(\mathcal{H}\right)\right]\\
\sum_{n_{1}=0}^{M-1}\!\cdots\!\sum_{n_{N}=0}^{M-1}G_{1}\!\left(\left\{ \lambda_{i}^{1}\right\};n_1\right)\cdots G_{N}\!\left(\left\{ \lambda_{i}^{N}\right\};n_N\right)\left\langle\Psi\right|T_{n_{1}}\!\left(\mathcal{H}\right)\mathcal{O}_{1}\cdots\mathcal{O}_{N}T_{n_{N}}\!\left(\mathcal{H}\right)\left|\Psi\right\rangle
\end{cases}\,\,,
$$

: where $T_{n}(x)$ are Chebyshev polynomials of the $1^{\text{st}}$-kind and $G_{j}$ are the expansion coefficients of $F_{j}$. The advantage gained by using the Chebyshev expansion is that both $\text{Tr}\left[\cdots\right]$ and $\left\langle\Psi\right|\cdots \left|\Psi\right\rangle$ can be evaluated recursively using only matrix-vector operations[^1]. For all observables implemented in KITE, the functions $F_{j}$ are of three types: 

: 1. Single-Parameter Dirac-$\delta$ Functions.—$F_{j}(\lambda_{1}^{j};\mathcal{H})\to\delta\left(\lambda-\mathcal{H}\right)$
  2. Broadened Single-Particle Green's Functions.—$F_{j}(\lambda_{1}^{j},\lambda_{2}^{j};\mathcal{H})\to\left[\lambda+i\eta-\mathcal{H}\right]^{-1}$
  3. Quantum Time-Evolution Operators.—$F_{j}\left(\lambda_{1}^{j},\mathcal{H}\right)\to\exp\left(\frac{i\,t\,\mathcal{H}}{\hbar}\right)$

: For these functions, analytical forms of the Chebyshev expansion coefficients are known[^2][^3][^4][^5][^6][^7] and used in KITE. In the user interface, the truncation order $M$ is specified by the parameter `#!python num_moments`, and always impacts the validity of the expanded results. Nevertheless, its precise effect depends crucially on the specific case, as shown in Fig. 1 below. We will now discuss three common cases. 

## Dirac-delta Function

: An order-$M$ expansion (regularized by the Jackson kernel) produces a Gaussian approximation of $\delta(\lambda\!-\!\mathcal{H})$ endowed by a width $\sigma_{\lambda}\!\approx\!\delta \varepsilon\,\pi/M$ in $\lambda$[^2][^7]. 
The choice of $M$ then fixes the effective spectral width, $\sigma_{\lambda}$, which must be sufficiently narrow to accurately describe all relevant features of the calculated property. 
However, if $\sigma_{\lambda}$ becomes too narrow ($M$ too high), the discrete eigenvalues of the SPH are well-resolved and the obtained data suffers from large (finite-size) fluctuations. For information on other available kernels see Weisse *et al.*[^2].

    !!! Info "Rule of Thumb"

        If $\Delta\varepsilon$ is the mean-level spacing of the simulated system, then $M$ must be kept **smaller than** $\frac{\pi\,\delta\varepsilon}{\Delta\varepsilon}$ in order to avoid resolving individual energy levels. Simultaneously, for obtaining high-resolution results, the artificial broadening must remain much smaller than the lattice bandwidth.

## Single-Particle Green's Functions

: No kernel is required for lattice Green's functions [^3][^4]. Note that these functions _must_ be broadened by a finite $\eta$ due to the discrete nature of the energy spectrum of finite systems. An exact spectral decomposition of broadened lattice Green's functions exists and is the basis of the Chebyshev polynomial Green's function method [^3] implemented in KITE. Provided $\eta$ exceeds the spacing between eigenvalues of the SPH, the truncation order may be arbitrarily increased, and convergence is achieved **when the data ceases to depend on M**. 

    !!! Info "Rule of Thumb"

        Since the energy resolution is fixed by $\eta$, the number of polynomials must be larger enough to resolve such a broadening. As shown below in Fig. 1(b), an apt rule of thumb is to have $M \gtrsim 10*\delta \varepsilon / \eta$.

## Quantum Time-Evolution Operators

: In time evolution problems, Chebyshev truncation errors propagate and limit the computation accuracy after some finite time has elapsed. 
In Fig. 1(c), it is demonstrated that the time evolution operator is converged as long as $t\!\lesssim\!\hbar\,M/\delta\varepsilon$.

    !!! Info "Rule of Thumb"

        A safe empirical rule of thumb (used in Santos Pires *et al.*[^6]) is to have $M\!\gtrsim\!8\,\delta\varepsilon\,\hbar^{-1}t_{\text{max}}\!\!$, where $t\!<\!t_{\text{max}}$ is the length of the time-interval intended for the evolution.

: For some target functions, the output of [`#!bash KITEx`](../api/kitex.md) will span energy and spatial coordinates (as is the case with the LDoS or the spectral function).
Note that some response functions (such as the transverse conductivity) will reflect the properties of both the Fermi surface and the Fermi sea of occupied states, and therefore require the raw output of KITE to be numerically integrated over energy. This integration procedure is always done at the post-processing level by [`#!bash KITE-tools`](../api/kite-tools.md).

	!!! Warning "Post-Processing Integration"
	
		For target functions that require energy integrations, a key parameter is the **number of energy points** for which the integrand is evaluated. In order to adjust it at the post-processing level, the user can use the `#!python -E` flag of the [`#!python KITE-tools`](../api/kite-tools.md) executable.

    <div>
      <figure>
        <img src="../../assets/images/tutorial/conv.png" style="width: 40em;" />
        <figcaption>Figure 1: Convergence of the Chebyshev series for (a) a Dirac-$\delta$ function, (b) a single-particle Green's function (real part as inset), and (c)  the time-evolution operator at two different energies (top: 0.5 and bottom: 0.8).</figcaption>
      </figure>
    </div>

[^1]: **KITE:** high-performance accurate modelling of electronic structure and response functions of large molecules, disordered crystals and heterostructures, S. M. João, M. Anđelković, L. Covaci, T. G. Rappoport, João M. Viana Parente Lopes, and A. Ferreira, [R. Soc. open sci. 7, 191809 (2020)](https://royalsocietypublishing.org/doi/10.1098/rsos.191809).

[^2]: Kernel polynomial method, Alexander Weiße, Gerhard Wellein, Andreas Alvermann, and Holger Fehske. [Rev. Mod. Phys. 78, 275 (2006)](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.78.275).

[^3]: Critical delocalization of chiral zero energy modes in graphene, A. Ferreira and E. Mucciolo, [Phys. Rev. Lett. 115, 106601 (2015)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.115.106601).

[^4]: Numerical evaluation of Green's functions based on the Chebyshev expansion, A. Braun and P. Schmitteckert, [Phys. Rev. B 90, 165112 (2014)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.165112).

[^5]: An accurate and efficient scheme for propagating the time dependent Schrödinger equation, H. Tal-Ezer and R. Kosloff, [J. Chem. Phys. 81, 3967-3971 (1984)](https://aip.scitation.org/doi/10.1063/1.448136).

[^6]: Landauer transport as a quasisteady state on finite chains under unitary quantum dynamics, J. P. Santos Pires, B. Amorim, and J. M. Viana Parente Lopes, [Phys. Rev. B 101, 104203 (2020)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.101.104203).

[^7]: Spectral functions of one-dimensional systems with correlated disorder, N. A. Khan, J. M. Viana Parente Lopes, J. P. Santos Pires, J. M. B. Lopes dos Santos, [J. Phys.: Condens. Matt. 31, 175501 (2019)](https://iopscience.iop.org/article/10.1088/1361-648X/ab03ad/meta).

[configuration]: ../api/kite.md#configuration
[configuration-spectrum_range]: ../api/kite.md#configuration-spectrum_range
