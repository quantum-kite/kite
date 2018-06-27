# Getting started

## Code workflow

KITE has 3 different layers: a user interface (Python); a main program (C++); and a post-processing tool (C++). The TB model is built and exported together with the calculation settings to a [HDF5][1] file (`*.h5`), which is then used as I/O to the main program (called **KITEx**). The workflow is as follows:

* Export model and calculation settings to a `*.h5` file for **KITEx**;
* Run **KITEx**;
* Run the post-processing tools and visualise the data;

Below, we illustrate the use of the [*Pybinding*][2] package to build a TB model for a crystal and illustrate the basic funcionalities of KITE. Those already familiar with Pybinding may consider skipping to [Examples][3].

## Building and exporting a TB model

KITE's intuitive user interface is based on [Pybinding][2] with be-spoke features for complex disorder/fields modifications and target functions e.g., DoS, conductivity, etc.

### Importing the Pybinding package

If all the installation requirements are fulfilled, Pybinding package can be imported in the python script. In this tutorial, the required packages will be included with the following aliases:
``` python
import kite
import pybinding as pb
import numpy as np
import matplotlib.pyplot as plt
```
If you want to use pybinding predefined styles for visualizing the lattice, simply add to the script:

``` python
pb.pltutils.use_style()
```

### Building the model

The most important object for building a TB model is **pb.Lattice**, which carries the information about the crystal structure (lattice and basis) and hopping parameters. Pybinding also provides additional features based on the real-space information, as for example, the reciprocal vectors and the Brillouin zone. As a simple example, let us consider a square lattice with a single lattice site.

First, import all the packages (including KITE special features):
``` python
import kite
import pybinding as pb
import numpy as np
import matplotlib.pyplot as plt
```
The following syntax can be used to define primitive lattice vectors:
``` python
a1 = np.array([1,0]) # [nm] define the first lattice vector
a2 = np.array([0, 1]) # [nm] define the second lattice vector

lat = pb.Lattice(a1=a1, a2=a2) # define a lattice object
```
Add the desired lattice sites inside the unit cell:
``` python
lat.add_sublattices(
    # make a lattice site (sublattice) with a tuple
    # (name, position, and onsite potential)
    ('A', [0, 0], onsite[0])
)
```
By default the main cell has the index `[n1,n2] = [0, 0]`. The hoppings between neighboring sites can be added with the simple syntax:
``` python
lat.add_hoppings(
    # make an hopping between lattice site with a tuple
    # (relative unit cell index, site from, site to, hopping energy)
    ([1, 0], 'A', 'A', - 1 ),
    ([0, 1], 'A', 'A', - 1 )
)
```
Here, the relative indices `n1,n2` represent the number of integer steps needed to reach a neighboring cell starting from the main one.

*Important note*: When adding the hopping `(n, m)` between sites `n` and `m`, the conjugate hopping term `(m, n)` is added automatically and it is not allowed to add them twice.

Now we can plot the lattice:
``` python
lat.plot()
plt.show()
```
and visualize the Brillouin zone:
``` python
lat.plot_brillouin_zone()
plt.show()
```
For a crystal with two atoms per unit cell see our [graphene example][4]. For other examples and pre-defined lattices consult the [Pybinding][2] documentation.

## Calculation settings and KITE target functions

The following Python classes from KITE are used to define the target functions and calculation settings:

1. `Configuration`
2. `Calculation`
3. `Modification`

### Configuration

The class `Configuration` carries the following information:

* `divisions` - integer number that defines the number of decomposition parts in each spatial direction. **KITEx** implements a domain decomposition technique to divide the lattice into various partitions that are computed in parallel. To activate this feature set a number of decomposition parts larger than unit `nx * ny > 1`. For example:
``` python
nx = ny = 2
```
decomposes a 2D lattice into four regions of equal size. The product `nx * ny` equals the number of threads used by **KITEx** and thus must not exceed the number of avaliable cores in the computer.

The domain decomposition is optimized at the design level and allows a substantial speed up of multithreaded calculations. We recommend its usage.

* `length` - integer number of unit cells along the direction of lattice vectors:
``` python
lx = 256
ly = 256
```
The lateral size of the decomposed parts are given by **lx/nx** and **ly/ny** that need to be integer numbers.

* `boundaries` - boolean value. **True** for periodic boundary conditions and **False** for open boundary conditions. Currently, **KITEx** only accepts *periodic boundary conditions*.

* `is_complex` - boolean value. For optimisation purposes, **KITEx** only considers and stores complex data with the setting **is_complex=True**. **False** should be used for real symmetric Hamiltonians.

* `precision` - integer identifier of data type. **KITEx** allows users to define the precision of the calculation. Use **0** for float, **1** for double, and **2** for long double.

* `spectrum_range` - array of reals (OPTIONAL). By default KITE executes an automated rescaling of the Hamiltonian; see [Resources][5]. Advanced users can override this feature using `spectrum_range=[Emin,Emax]`, where `Emin(Emax)` are the minimum (maximum) eigenvalues of the TB matrix.

As a result, a `Configuration` object is structured in the following way:
``` python
configuration = ex.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True], is_complex=False, precision=1)
```
### Calculation

Finally, the `Calculation` object carries out the information about the quantities that are going to be calculated i.e., the CPGF target functions. For this part, we still need to include more parameters, related to the Chebyshev expansion (our examples already have optimized parameters for a standard desktop computer). All target functions require the following parameters:

1. **num_moments** defines the number of moments of the Chebyshev expansion and hence the energy resolution of the calculation; see [Resources][5].

2. **num_random** defines the number of random vectors for the stochastic evaluation of target functions; see [Resources][5].

3. **num_disorder** defines the number of disorder realisations.

The target function functions currently available are:

* `dos`
* `conductivity_optical`
* `conductivity_dc`
* `conductivity_optical_nonlinear`
* `singleshot_conductivity_dc`

with the following parameters:

* `direction` - direction along which the conductivity is calculated (longitudinal: 'xx', 'yy', transversal: 'xy', 'yx')
* `temperature` - temperature used in Fermi Dirac distribution that is used for the calculation of optical and DC conductivities.
* `num_points` - number of points the in energy axis that is going to be used by the post-processing tool to output the density of states.
* `special` - simplified form of nonlinear optical conductivity hBN example
* `energy` - selected value of energy at which we want to calculate the singleshot_conductivity_dc
* `eta` - imaginary term in the denominator of the Green function's that provides a controlled broadening / inelastic energy scale (for technical details, see [Resources][5]).

The **calculation** is structured in the following way:
``` python
calculation = Calculation(configuration)

calculation.dos(num_points=1000, num_random=10, num_disorder=1, num_moments=512)

calculation.conductivity_optical(num_points=1000, num_random=1, num_disorder=1, num_moments=512, direction='xx')

calculation.conductivity_dc(num_points=1000, num_moments=256, num_random=1, num_disorder=1,direction='xy', temperature=1)

calculation.singleshot_conductivity_dc(energy=[(n/100.0 - 0.5)*2 for n in range(101)], num_moments=256, num_random=1, num_disorder=1,direction='xx', gamma=0.02)

calculation.conductivity_optical_nonlinear(num_points=1000, num_moments=256, num_random=1, num_disorder=1,direction='xxx', temperature=1.0, special=1)
```
*Important note*: The user can decide what functions are used in a calculation. However, it is not possible to configure the same function twice in the same Python script (HDF5 file).

When these objects are defined, we can export the file that will contain set of input instructions for **KITEx**:
``` python
kite.export_lattice(lattice, configuration, calculation, 'test.h5')
```
### Running the code

To run the code and the postprocess it, use
``` python
./KITEx test.h5
./tools/KITEtools test.h5
```
# Visualizing the data

After calculating the quantity of interest and post-processing the data, we can plot the resulting data with the following script:

https://gist.github.com/quantum-kite/9a935269845eae3f8590f364be12cb49
 
!![image][6]

If you want to make these steps more automatic, you can use the following Bash script

https://gist.github.com/quantum-kite/c002610a4d43a478cf0f967129f97da7
 
[1]: https://www.hdfgroup.org/
[2]: http://docs.pybinding.site/en/stable/
[3]: https://quantum-kite.com/category/examples/
[4]: https://gist.github.com/quantum-kite/4b7593e9aa082b1d242c5e6b2361c3f3
[5]: http://quantum-kite.com/resources/
[6]: https://user-images.githubusercontent.com/39924384/41257230-fb584956-6dc3-11e8-9b4d-530b0cf3be6c.png
