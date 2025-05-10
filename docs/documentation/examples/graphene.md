## The TB-model for Graphene

The electronic structure of graphene is well-described by a simple tight-binding model that only uses one $p_z$ orbital
in a hexagonal unit cell with two equivalent carbon atoms.
These atoms are located on the different sublattices, *A* and *B*, and don't have an on-site energy term.

Although this model is relatively simple, it is used often within the literature. The example below uses graphene to show how to construct a lattice with sublattices and perform a basic calculation using KITE.
_The script for this example can be found [here](https://github.com/quantum-kite/kite/blob/313a00e54a9f9aa33b22886eaf97ce62aaec3996/examples/dos_graphene.py)._

### Lattice
We start by building the [`#!python pb.lattice`][lattice] for graphene:

* Define the parameter ($t$ in eV)
* Define the vectors of the unit-cell ($\vec a_1$ and $\vec a_2$ in units of $a$, length of the unit-cell)
* Create a [`#!python pb.lattice`][lattice]-object
* Define the on-site energies
* Define the hopping parameters
    * 1 *normal* hopping within the unit cell
    * 2 *rotated* hopping $\pm 2 \pi/3$ to neighbouring cells
* Return the [`#!python pb.lattice`][lattice]-object to be used by KITE

``` python linenums="1"
def graphene_lattice(onsite=(0, 0)):
    """Return lattice specification for a honeycomb lattice with nearest neighbor hoppings"""

    # parameters
    a = 0.24595  # [nm] unit cell length
    a_cc = 0.142  # [nm] carbon-carbon distance
    t = 2.8  # eV

    # define lattice vectors
    a1 = a * np.array([1, 0])
    a2 = a * np.array([1 / 2, 1 / 2 * np.sqrt(3)])

    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(a1=a1, a2=a2)

    # add sublattices
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [0, -a_cc/2], onsite[0]),
        ('B', [0,  a_cc/2], onsite[1])
    )

    # Add hoppings
    lat.add_hoppings(
        # inside the main cell, between which atoms, and the value
        ([0, 0], 'A', 'B', -t),
        # between neighboring cells, between which atoms, and the value
        ([1, -1], 'A', 'B', -t),
        ([0, -1], 'A', 'B', -t)
    )
    return lat
```

We can visualize this lattice using the following code:

``` python linenums="1"
lat = graphene_lattice()
lat.plot()
plt.show()
```

<div>
  <figure>
    <img src="../../../assets/images/getting_started/graphene_example.png" style="width: 40em;"/>
    <figcaption>A visualization for the defined lattice of graphene.</figcaption>
  </figure>
</div>

## KITEx calculation
### Settings and calculation
We can make the [`#!python kite.Calculation`][calculation] object
``` python linenums="1"
configuration = kite.Configuration(
  divisions=[64, 64],
  length=[512, 512],
  boundaries=["periodic", "periodic"],
  is_complex=False,
  precision=1
)
calculation = kite.Calculation(configuration)
calculation.dos(
  num_points=4000,
  num_moments=256,
  num_random=256,
  num_disorder=1
)
```
Set up the `#!python kite.config_system` as detailed in Sec. [Calculation][calculation] and export the KITE model running
``` bash
python3 script_name_here.py
```
which then creates the necessary HDF file. Next, run the [KITEx][kitex] program and the [KITE-tools][kitetools].

### Visualization

``` python linenums="1"
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('dos.dat')
plt.plot(data[:,0], data[:,1])
plt.show()
```

<div>
  <figure>
    <img src="../../../assets/images/getting_started/first_calculation.png" style="width: 40em;" />
    <figcaption>The DOS for graphene.</figcaption>
  </figure>
</div>

[tutorial]: ../index.md
[calculation]: ../calculation.md
[getting_started]: ../index.md
[lattice]: https://docs.pybinding.site/en/stable/_api/pybinding.Lattice.html
[kitex]: ../../api/kitex.md
[kitetools]: ../../api/kite-tools.md