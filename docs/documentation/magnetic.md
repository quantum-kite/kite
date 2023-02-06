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

!!! Info "More examples"
    
    In the section [Examples], some applications to different structures are given, including the caluclation.
    More examples can be found in the KITE-repository under
    [kite/examples](https://github.com/quantum-kite/kite/tree/master/examples/readme.md).
