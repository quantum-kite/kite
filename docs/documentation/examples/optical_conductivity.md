## Gaussian disorder
KITE also calculates the optical conductivity of a given lattice for a given Fermi energy.
To illustrate this capability, we calculate the optical conductivity of disordered graphene, that can be compared qualitatively with previous results[^1].

### lattice
Instead of defining the lattice in our python script, we can use of one the pre-defined lattices from Pybinding:

```python
    from pybinding.repository import graphene
    lattice = graphene.monolayer()
```

## KITE-part
### Disorder
To illustrate a different type of disorder, we add random on-site energies that follow a Gaussian distribution

```python 
        disorder = kite.Disorder(lattice)
    disorder.add_disorder('B', 'Gaussian', 0.0,1.5)
    disorder.add_disorder('A', 'Gaussian',  0.0, 1.5)
```

where we define the type of statistical distribution, the sublattices they are located, the mean value of the distribution and its width.

### Settings
After configuring the system, as presented in section [Getting Started](../index.md), it is time to set the calculation:

```python 
    calculation = kite.Calculation(configuration)
    calculation.conductivity_optical(num_points=1000, num_disorder=1,
                     num_random=20, num_moments=512, direction='xx')
```

The optical conductivity can also be calculated in different directions, which can be quite interesting in the case of Hamiltonians with non-trivial topology that also present transverse optical conductivity.
However, in this example we focus on longitudinal optical conductivity.
The other quantities that can be set in the python script are the same as for the calculation of the density of states: number of energy points used in KITE-tools, moments in the expansion, random vectors and disorder realisations.

### Calculation
When calculating the optical conductivity, it is also necessary to set the Fermi energy.
However, in this BETA version, the executable KITE-tools only calculates the optical conductivity in one pre-defined Fermi energy.
To change it, it is necessary to edit the source code and recompile KITE-tools.

### Visualization
The results of the real and imaginary parts of the optical conductivity presented in the first figure were obtained on a normal desktop with calculations that took 8 minutes for a system with `N=512 x 512` units cells and 512 expansion moments.

![optcond](../../assets/images/optical_conductivity/optcond.png)

For systems sizes of `N=1536 x 1536` unit cells and two different Gaussian widths, we show $\Re [\sigma_{xx}(\omega)]$ for low frequencies, where we can see Drude's peak for $\omega\rightarrow 0$ and the onset of interband transitions at $\hbar \omega>2 E_F$ [^2].

![optcond_Re](../../assets/images/optical_conductivity/optcond_Re.png)

The complete python script for this calculation can be found [here](https://gist.github.com/quantum-kite/bab11b55ad6672f59f15a85aac800a3a).


[^1]: Shengjun Yuan, Rafael Roldán, Hans De Raedt, Mikhail I. Katsnelson, [Phys. Rev. B **84**, 195418 (2011)](https://link.aps.org/doi/10.1103/PhysRevB.84.195418).

[^2]: T. Stauber, N. M. R. Peres, A. K. Geim, [Phys. Rev. B **78**, 085432 (2008)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.78.085432).

