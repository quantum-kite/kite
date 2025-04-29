Static magnetic fields are an important case of lattice modifications that can be performed automatically by KITE. 
This feature is of particular interest to the study of Landau levels, magneto-transport and magneto-optical effects.

## Uniform B-fields in KITE: overview

The automated $\mathbf{B}$-field functionality works by the addition of Peierls phases in the Hamiltonian and can be
used in conjunction with other lattice modifications, including disorder.
This is a new feature under development that currently allows for

* uniform $\mathbf{B}$ fields in 2D lattices (with the B-field perpendicular to the plane)
* uniform $\mathbf{B}$ fields in 3D lattices (with the B-field collinear to the third primitive lattice vector)

The $\mathbf{B}$-field is added by means of the following [KITE modification][modification-par-magnetic_field]:

``` py
modification = kite.Modification(magnetic_field = mag)
``` 

where (`#!python mag`) is the magnetic field strength (in Tesla).
When used with [periodic boundary conditions][configuration-boundaries],
$|\mathbf{B}|$ is restricted to be a multiple of a minimum magnetic field,
which is determined internally when generating the configuration file (see details below). 

For example, to compute the DOS of a disordered system subject to the
$\mathbf{B}$-field modification discussed above one may use

``` py linenums="1"
calculation.dos(
    num_points=100,
    num_moments=5000,
    num_random=10,
    num_disorder=10
)
mod_mag_field = kite.Modification(magnetic_field = mag)
kite.config_system(
    lattice,
    configuration,
    calculation,
    modification=mod_mag_field,
    disorder=disorder,
    filename="B_field.h5"
)
``` 
Here, the configuration file requests a DOS calculation with `#!python 5000` Chebyshev moments, `#!python 10`
random vectors and `#!python 10` disorder realizations.
Note the addition of the lattice modification (`#!python modification=mod_mag_field`) with respect to the settings
as discussed in [the section about disorder](disorder.md).

## Implementation details

!!! Info "Units"  

    Lattice parameters must be given in nanometers (see example `#!python dos_dccond_square_lattice.py`
    in the [`#! bash kite/examples/`-folder][examples-folder]).

The magnetic fields considered in KITE are uniform, so the corresponding vector potential is linear. It is naturally expressed in terms of the primitive reciprocal lattice vectors ($\mathbf{b}_{i=1,2,3}$) in the Landau gauge

$$
\mathbf{A}\left(\mathbf{r}\right)=\frac{h}{\left(2\pi\right)^{2}e}\frac{n}{N_{2}}\left(\mathbf{r}\cdot\mathbf{b}_{2}\right)\mathbf{b}_{1}
$$

where h is Planck's constant, $e>0$ is the elementary charge, $N_{2}$ is the number of unit cells along the $\mathbf{a}_{2}$ direction (primitive vector of the direct lattice) and $n$ is an integer. The corresponding magnetic field points along the $\mathbf{a}_{3}$ direction for 3D systems and perpendicularly to the basal plane ($\mathbf{e}_{\perp}\equiv\hat{\mathbf{z}}$) for 2D systems:


$$
\mathbf{B}=\begin{cases}
\frac{h}{e\Omega_{c}}\frac{n}{N_{2}}\mathbf{a}_{3} & (\textrm{3D})\\
\\
\frac{h}{e\Omega_{c}}\frac{n}{N_{2}}\hat{\mathbf{z}} & (\textrm{2D})
\end{cases}
$$

and is restricted to be a multiple of a minimum field 

$$
B_{\textrm{min}}=\frac{h}{e\Omega_{c}}\frac{1}{N_{2}}\times\begin{cases}
|\mathbf{a}_{3}| & (\textrm{3D})\\
\\
1 & (\textrm{2D})
\end{cases}
$$

where $\Omega_{c}$ is the 3D/2D volume of the unit cell.  When the user requests a magnetic field strength $|\mathbf{B}|$ (in Tesla), KITE calculates $B_{\textrm{min}}$ first and then uses that to determine the required $n$ to achieve the closest possible value of $|\mathbf{B}|$ by rounding $|\mathbf{B}|/\mathbf{B}_{\textrm{min}}=n$ to the nearest integer. If $n$ rounds down to zero, it means that the system is too small to support the requested magnetic field. When determining $B_{\textrm{min}}$, _KITE assumes that the primitive vectors in the Python configuration script are given in nanometers_.



[modification-par-magnetic_field]: ../api/kite.md#modification-par-magnetic_field
[configuration-boundaries]: ../api/kite.md#configuration-boundaries
[examples-folder]: more_examples/additional_examples.md
