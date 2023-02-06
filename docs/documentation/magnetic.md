Static magnetic fields are an important case of lattice modifications that can be performed automatically by KITE. This feature is of particular interest to the study of Landau levels, magneto-transport and magneto-optical effects.

## Uniform B-fields in KITE: overview

The automated B-field functionality works by the addition of Peierls phases in the Hamiltonian and can be used in conjunction with other lattice modifications, including disorder. This is a new feature under development that currently allows

* uniform B-field in generic 2D lattices (with the B-field perpendicular to the plane)
* uniform B-field in 3D lattices with cubic symmetry (with the B-field parallel to the second primitive lattice vector)

The B-field is added by means of the following KITE modification:

``` python
modification = kite.Modification(magnetic_field = mag)
``` 

where (`#!python mag` ) is the magnetic field strength. When used with periodic boundary conditions, B is restricted to be a multiple of a minimum magnetic field, which is determined internally when generating the configuration file (see details below). 

For example, to compute the DOS of a disordered system subject to the B-field modification discussed above we use

``` python
calculation.dos(num_points=100, num_moments= 5000, num_random=10, num_disorder=10)
kite.config_system(lattice, configuration, calculation, modification=modification, disorder=disorder,filename='B_field.h5")
``` 

## Implementation details

!!! Info "Units"  
    When using the B-field feature lattice parameters must be given in nanometers (see example "..." in the KITE folder (link here)). 
!!!

The magnetic fields considered in KITE are uniform, so the corresponding vector potential is linear. It is naturally expressed in terms of the primitive reciprocal lattice vectors ($\mathbf{b}_{i=1,2,3}$) in the Landau gauge

$$
\mathbf{A}\left(\mathbf{r}\right)=\frac{h}{\left(2\pi\right)^{2}e}\frac{n}{N_{2}}\left(\mathbf{r}\cdot\mathbf{b}_{2}\right)\mathbf{b}_{1}
$$

where h is Planck's constant, $e>0$ is the elementary charge and $N_{2}$ is the number of unit cells along the $\mathbf{a}_{2}$ direction (primitive vector of the direct lattice) and $n$ is an integer. The corresponding magnetic field points in the $\mathbf{a}_{3}$ direction.

$$
\mathbf{B}=\frac{h}{e\Omega_{c}}\frac{n}{N_{2}}\mathbf{a}_{3}
$$

and is restricted to be a multiple of a minimum field 

$$
B_{\textrm{min}}=\frac{h}{e\Omega_{c}}\frac{1}{N_{2}}\left|\mathbf{a}_{3}\right|,
$$

where $\Omega_{c}$ is the volume of the unit cell. When the user requests a magnetic field strength $|\mathbf{B}|$, KITE calculates $B_{\textrm{min}}$ first and then uses that to determine the required n to achieve the closest possible value of $\mathbf{B}$ by rounding $\mathbf{B}|/B_{\textrm{min}}=n$ to the nearest integer. If n rounds down to zero, it means that the system is too small to support the requested magnetic field. When determining $B_{\textrm{min}}$, KITE assumes that the primitive vectors in the Python configuration script are given in nanometers.


