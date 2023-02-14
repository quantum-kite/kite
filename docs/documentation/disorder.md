The disorder implementation — general-purpose and user-friendly — is one of the main features of the KITE code.
The inclusion of disorder in a given system follows a very simple procedure:
the user specifies one or more _disorder patterns_ in the [Python interface][kite-script]
(patterns are local modifications of the Hamiltonian that can be constricted to one unit cell or can
connect neighboring unit cells) together with the desired disorder statistics.
The disorder patterns are then replicated across the lattice automatically when running [KITEx].
KITE handles both standard uncorrelated disorder (e.g., random on-site energies) and realistic
short-range disorder (e.g., vacancies or impurity scattering centers distributed randomly over the
lattice sites with a specified concentration).
Several types of disorder can be added simultaneously.

After defining a [regular lattice][tutorial-lattice], disorder can be added to the system.
KITE allows the user to select between on-site and structural disorder by choosing between predefined classes
in the [python interface][kite-script].
The interface provides two different classes of disorder:

* [`#!python kite.Disorder`][disorder] - onsite disorder with three possible statistical distributions
* [`#!python kite.StructuralDisorder`][structural_disorder] - generic short-range disorder, i.e. multi-orbital local disorder (including bond disorder) with a given concentration

## Onsite disorder

[`#!python kite.Disorder`][disorder] adds randomly generated onsite terms at the sites of a desired sublattice based on a certain statistical distribution:

* Gaussian
* Uniform
* Deterministic

Beside the type of statistical distribution, we can select a sublattice type in which the disorder will appear, the mean
value and the standard deviation of the selected distribution.
To include onsite disorder following a given statistical distribution, we build the [`#!python pb.lattice`][lattice] and
use the following procedure:

``` python
# define an object based on the lattice
disorder = kite.Disorder(lattice)
# add Gaussian distributed disorder at all sites of a selected sublattice
disorder.add_disorder('A', 'Gaussian', 0.1, 0.1)
```

In a single object it is possible to select multiple sublattices, each of one with different disorder distributions
following the rule [`#!python kite.Disorder.add_disorder('sublattice', 'type', mean, std)`][disorder-add_disorder] :

``` python
disorder.add_disorder('A', 'Gaussian', 0.1, 0.1)
disorder.add_disorder('B', 'Uniform', 0.2, 0.1)
disorder.add_disorder('C', 'Deterministic', 0.1)
```

In the case of deterministic disorder, the standard deviation is not set.

After defining the desired disorder, it can be added to the configuration file as an additional parameter in the
[`#!python kite.config_system`][config_system] function:

``` python
kite.config_system(..., disorder=disorder)
```

A complete example that calculates the density of states of graphene with different on-site disorder
distributions for each sublattice can be seen here:

``` python linenums="1"
import kite
from pybinding.repository import graphene

""" Onsite disorder
    Lattice : Monolayer graphene;
    Disorder : Disorder class Deterministic and Uniform at different sublattices,
    Configuration : size of the system 512×512, without domain decomposition (nx=ny=1),
                    periodic boundary conditions,
                    double precision, automatic scaling;
    Calculation : dos;
    Modification : magnetic field is off;
"""

# load graphene lattice and structural_disorder
lattice = graphene.monolayer()
# add Disorder
disorder = kite.Disorder(lattice)
disorder.add_disorder('B', 'Deterministic', -1.0)
disorder.add_disorder('A', 'Uniform', +1.5, 1.0)
    # number of decomposition parts [nx,ny] in each direction of matrix.
# This divides the lattice into various sections, each of which is calculated in parallel
nx = ny = 2
# number of unit cells in each direction.
lx = ly = 512

# make config object which caries info about
# - the number of decomposition parts [nx, ny],
# - lengths of structure [lx, ly]
# - boundary conditions [mode,mode, ... ] with modes:
#   . "periodic"
#   . "open"
#   . "twisted" -- this option needs the extra argument angles
#   . "random"

# Boundary Mode
mode = "periodic"
configuration = kite.Configuration(
    divisions=[nx, ny],
    length=[lx, ly],
    boundaries=[mode, mode],
    is_complex=False,
    precision=1
)
# require the calculation of DOS
calculation = kite.Calculation(configuration)
calculation.dos(num_points=5000, num_moments=512, num_random=1, num_disorder=1)
# configure the *.h5 file
kite.config_system(lattice, configuration, calculation, filename='on_site_disorder.h5', disorder=disorder)
``` 
<div>
  <figure>
    <img src="../../assets/images/disorder/on_site_disorder_dos.png" width="300" style="display: inline-block;"/>
    <figcaption>DOS for the on-site disorder.</figcaption>
  </figure>
</div>


## Structural disorder

[`#!python kite.StructuralDisorder`][structural_disorder] class adds the possibility of selecting between two different
structural disorder types; vacancy, randomly distributed with a certain concentration in sites of a selected sublattice,
and a more generic structural disorder which is a combination of onsite terms and bond disorder
(also distributed with a certain concentration).

### Vacancy disorder

The vacant site distribution can be selected from a single sublattice with a
[`#!python concentration`][structuraldisorder-concentration] defined in a parent object:

``` python
struc_disorder = kite.StructuralDisorder(lattice, concentration=0.2)
struc_disorder.add_vacancy('B') # add a vacancy to a selected sublattice
```

!!! note
    To distribute the vacancies in both sublattices, one needs to add the vacancies on each sublattice as a separate
    object of the class [`#!python kite.StructuralDisorder`][structural_disorder]
    (unless you want precisely the same pattern of disorder in both sublattices).

    ``` python
    struc_disorder_A = kite.StructuralDisorder(lattice, concentration=0.1)
    struc_disorder_A.add_vacancy('A')
    struc_disorder_B = kite.StructuralDisorder(lattice, concentration=0.1)
    struc_disorder_B.add_vacancy('B')
    disorder_structural = [struc_disorder_A, struc_disorder_B]
    ```


### Structural disorder

To manyally set the [`#!python spectrum_range`][configuration-spectrum_range], it is necessary to add an extra parameter
to the [`#!python kite.Configuration`][configuration] class:

``` python
configuration = kite.Configuration(
    divisions=[nx, ny],
    length=[lx, ly],
    boundaries=["periodic", "periodic"],
    is_complex=False,
    precision=1,
    spectrum_range=[-10, 10]
)
```

The following example shows a definition of our most general type of disorder, which includes both onsite disorder
terms and bond modifications.
This type of disorder can be added as an object of the class [`#!python kite.StructuralDisorder`][structural_disorder].
The procedure for adding the structural disorder is the same of adding a hopping term to the
[Pybinding lattice object][lattice],
with a single difference that the bond disorder is not bounded to the hopping term starting from the `#!python [0, 0]` unit cell,
which is the case of the hopping term in [Pybinding].

For the sake of clarity, let us first define sublattices that will compose the disorder.
In this case we are not restricted to a single unit cell:

``` python
#  define a node in a unit cell [i, j] selecting a single sublattice
node0 = [[+0, +0], 'A']
node1 = [[+0, +0], 'B']
node2 = [[+1, +0], 'A']
node3 = [[+0, +1], 'B']
node4 = [[+0, +1], 'A']
node5 = [[-1, +1], 'B']
```

After the definition of a parent [`#!python kite.StructuralDisorder`][structural_disorder] object, we can select the desired pattern:

``` python
 # define an object based on the lattice with a certain concentration
struc_disorder = kite.StructuralDisorder(lattice, concentration=0.2)

struc_disorder.add_structural_disorder(
    # add bond disorder in the form
    #  [from unit cell], 'sublattice_from', [to_unit_cell], 'sublattice_to', value:
    (*node0, *node1, 0.5),
    (*node1, *node2, 0.1),
    (*node2, *node3, 0.5),
    (*node3, *node4, 0.3),
    (*node4, *node5, 0.4),
    (*node5, *node0, 0.8),
    # in this way we can add onsite disorder in the form [unit cell], 'sublattice', value
    ([+0, +0], 'B', 0.1)
)

# It is possible to add multiple different disorder types which
#  should be forwarded to the config_system function as a list.
another_struc_disorder = kite.StructuralDisorder(lat, concentration=0.6)
another_struc_disorder.add_structural_disorder(
    (*node0, *node1, 0.05),
    (*node4, *node5, 0.4),
    (*node5, *node0, 0.02),
    ([+0, +0], 'A', 0.3)
)
```

Before exporting the settings to the [HDF5]-file, it is possible to define multiple disorder realizations which will be
superimposed to the clean system.

The following script has a minimal example of how to configure the structural disorder

``` python linenums="1"
""" Bond disorder
    Lattice : Honeycomb 1[nm] interatomic distance and t=1[eV] hopping;
    Disorder : StructuralDisorder class bond and vacancy disorder;
    Configuration : size of the system 512x512, without domain decomposition (nx=ny=1),
                    periodic boundary conditions,
                    double precision, manual scaling;
    Calculation : dos;
    Modification : magnetic field is off;
"""

import numpy as np
import pybinding as pb
import kite


def honeycomb_lattice(onsite=(0, 0)):
    """Make a honeycomb lattice with nearest neighbor hopping"""

    theta = np.pi / 3
    a1 = np.array([1 + np.cos(theta), np.sin(theta)])
    a2 = np.array([0, 2 * np.sin(theta)])

    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(
        a1=a1,
        a2=a2
    )

    # Add sublattices
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [0, 0], onsite[0]),
        ('B', [1, 0], onsite[1])
    )

    # Add hoppings
    lat.add_hoppings(
        # inside the main cell, between which atoms, and the value
        ([0, 0], 'A', 'B', - 1),
        # between neighboring cells, between which atoms, and the value
        ([-1, 0], 'A', 'B', - 1),
        ([-1, 1], 'A', 'B', - 1),
    )

    # Add bond disorder as an object of a class StructuralDisorder.
    # In this manner we can add onsite and bond defects
    # with a specific concentration, which will be added to the simulated system.
    # The procedure for adding is same as adding the hopping, with the difference
    # that the bond disorded is not bounded to one site in the [0, 0]
    # unit cell.
    node0 = [[+0, +0], 'A']
    node1 = [[+0, +0], 'B']
    node2 = [[+1, +0], 'A']
    node3 = [[+0, +1], 'B']
    node4 = [[+0, +1], 'A']
    node5 = [[-1, +1], 'B']

    struc_disorder_one = kite.StructuralDisorder(lat, concentration=0.05)
    struc_disorder_one.add_structural_disorder(
        # add bond disorder in the form
        #  [from unit cell], 'sublattice_from', [to_unit_cell], 'sublattice_to', value:
        (*node0, *node1, 1),
        (*node1, *node2, 1),
        (*node2, *node3, 1),
        (*node3, *node4, 1),
        (*node4, *node5, 1),
        (*node5, *node0, 1),
        # in this way we can add onsite disorder in the form [unit cell], 'sublattice', value
        ([+0, +0], 'B', 0.3)
    )
    # It is possible to add multiple different disorder type which
    # should be forwarded to the export_lattice function as a list.
    struc_disorder_two = kite.StructuralDisorder(lat, concentration=0.2)
    struc_disorder_two.add_structural_disorder(
        (*node0, *node1, 0.4),
        (*node4, *node5, 0.4),
        (*node5, *node0, 0.4),
        ([+0, +0], 'B', 0.4)
    )
    struc_disorder_two.add_vacancy('B')

    struc_disorder_three = kite.StructuralDisorder(lat, concentration=0.01)
    struc_disorder_three.add_vacancy('A')

    # if there is disorder it should be returned separately from the lattice
    return lat, [struc_disorder_one, struc_disorder_two, struc_disorder_three]


# load a honeycomb lattice and structural_disorder
lattice, disorder_structural = honeycomb_lattice()
# number of decomposition parts in each direction of matrix.
# This divides the lattice into various sections, each of which is calculated in parallel
nx = ny = 1
# number of unit cells in each direction.
lx = ly = 512
# make config object which caries info about
# - the number of decomposition parts [nx, ny],
# - lengths of structure [lx, ly]
# - boundary conditions, setting True as periodic boundary conditions, and False elsewise,
# - info if the exported hopping and onsite data should be complex,
# - info of the precision of the exported hopping and onsite data,
#      0 - float, 1 - double, and 2 - long double.
# - scaling, if None it's automatic, if present select spectrum_bound=[e_min, e_max]
configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                   is_complex=False, precision=1, spectrum_range=[-15, 15])
# require the calculation of DOS
calculation = kite.Calculation(configuration)
calculation.dos(num_moments=1024, num_random=1, num_disorder=1, num_points=1000)
# configure the *.h5 file
kite.config_system(lattice, configuration, calculation, filename='structural_disorder.h5',
                   disorder_structural=disorder_structural)
```

with the resulting density of states:

<div>
  <figure>
    <img src="../../assets/images/disorder/dos.png" width="300" style="display: inline-block;"/>
    <figcaption>DOS for the structural disorder.</figcaption>
  </figure>
</div>

[HDF5]: https://www.hdfgroup.org
[Pybinding]: https://docs.pybinding.site/en/stable
[lattice]: https://docs.pybinding.site/en/stable/_api/pybinding.Lattice.html
[tutorial-lattice]: tb_model.md
[comment]: <> (Class StructuralDisorder)
[structural_disorder]: ../api/kite.md#structuraldisorder
[comment]: <> (Class Parameters)
[structuraldisorder-lattice]: ../api/kite.md#structuraldisorder-lattice
[structuraldisorder-concentration]: ../api/kite.md#structuraldisorder-concentration
[structuraldisorder-position]: ../api/kite.md#structuraldisorder-position
[comment]: <> (Class Methods)
[structuraldisorder-add_vacancy]: ../api/kite.md#structuraldisorder-add_vacancy
[structuraldisorder-add_structural_disorder]: ../api/kite.md#structuraldisorder-add_structural_disorder
[structuraldisorder-add_local_vacancy_disorder]: ../api/kite.md#structuraldisorder-add_local_vacancy_disorder
[structuraldisorder-add_local_bond_disorder]: ../api/kite.md#structuraldisorder-add_local_bond_disorder
[structuraldisorder-add_local_onsite_disorder]: ../api/kite.md#structuraldisorder-add_local_onsite_disorder
[structuraldisorder-map_the_orbital]: ../api/kite.md#structuraldisorder-map_the_orbital

[kite-script]: ../api/kite.md
[KITEx]: ../api/kitex.md

[comment]: <> (Class Disorder)
[disorder]: ../api/kite.md#disorder
[comment]: <> (Class Parameters)
[disorder-lattice]: ../api/kite.md#disorder-lattice
[comment]: <> (Class Methods)
[disorder-add_disorder]: ../api/kite.md#disorder-add_disorder
[disorder-add_local_disorder]: ../api/kite.md#disorder-add_local_disorder
[config_system]: ../api/kite.md#config_system

[comment]: <> (Class Configuration)
[configuration]:  ../api/kite.md#configuration
[comment]: <> (Class Parameters)
[configuration-divisions]:  ../api/kite.md#configuration-divisions
[configuration-length]:  ../api/kite.md#configuration-length
[configuration-boundaries]:  ../api/kite.md#configuration-boundaries
[configuration-is_complex]:  ../api/kite.md#configuration-is_complex
[configuration-precision]:  ../api/kite.md#configuration-precision
[configuration-spectrum_range]:  ../api/kite.md#configuration-spectrum_range
[configuration-angles]:  ../api/kite.md#configuration-angles
[configuration-custom_local]:  ../api/kite.md#configuration-custom_local
[configuration-custom_local_print]:  ../api/kite.md#configuration-custom_local_print
[comment]: <> (Class Attributes)
[configuration-energy_scale]:  ../api/kite.md#configuration-energy_scale
[configuration-energy_shift]:  ../api/kite.md#configuration-energy_shift
[configuration-comp]:  ../api/kite.md#configuration-comp
[configuration-prec]:  ../api/kite.md#configuration-prec
[configuration-div]:  ../api/kite.md#configuration-div
[configuration-bound]:  ../api/kite.md#configuration-bound
[configuration-leng]:  ../api/kite.md#configuration-leng
[configuration-type]:  ../api/kite.md#configuration-type
[configuration-custom_pot]:  ../api/kite.md#configuration-custom_pot
[configuration-print_custom_pot]:  ../api/kite.md#configuration-print_custom_pot
[comment]: <> (Class Methods)
[configuration-set_type]:  ../api/kite.md#configuration-set_type
