!!! Note

    In this and next tutorials, the required packages will be included with the following aliases
    
    ``` python
    import kite
    import pybinding as pb
    import numpy as np
    import matplotlib.pyplot as plt
    ```

We will construct a [`#!python pb.Lattice`][lattice] using [Pybinding]
and calculate its DOS using [`#!python pb.Solver`][lattice] from [Pybinding].

!!! Info

    If you are familiar with [Pybinding], you can go directly to the next tutorial page.

## Making a `#!python pb.Lattice`
The [`#!python pb.Lattice`][lattice]-class from [Pybinding] carries the information about the [TB model][tightbinding].
This includes

* [crystal structure (*lattice* and *basis*)][unitcell]
* [onsite energies][onsite]
* [hopping parameters][hopping]

[Pybinding] also provides additional features based on the *real-space* information,
like the *reciprocal vectors* and the *Brillouin zone*.

### Unit-cell
As an example, we construct a square lattice with a singe lattice site.
The following syntax can be used to define primitive lattice vectors:

``` python linenums="1"
a1 = np.array([1, 0]) # [nm] define the first lattice vector
a2 = np.array([0, 1]) # [nm] define the second lattice vector

lat = pb.Lattice(a1=a1, a2=a2) # define a lattice object
```

### Lattice sites
Add the desired lattice sites inside the unit cell:

``` python linenums="1"
onsite = 0 # onsite potential
lat.add_sublattices(
    # make a lattice site (sublattice) with a tuple
    # (name, position, and onsite potential)
    ('A', [0, 0], onsite)
)
```

### Adding Hoppings
By default, the main cell has the index `#!python [n1,n2] = [0, 0]`.
The hoppings between neighboring sites can be added with the simple syntax:

``` python linenums="1"
lat.add_hoppings(
    # make an hopping between lattice site with a tuple
    # (relative unit cell index, site from, site to, hopping energy)
    ([1, 0], 'A', 'A', - 1 ),
    ([0, 1], 'A', 'A', - 1 )
)
```

Here, the relative indices `#!python n1,n2` represent the number of integer steps needed to reach a neighboring cell starting from the main one.

!!! note

    When adding the hopping `#!python (n, m)` between sites `#!python n` and `#!python m`,
    the conjugate hopping term `#!python (m, n)` is added automatically and it is not allowed to add them twice.

### Visualization
Now we can plot the [`#!python pb.lattice`][lattice] and visualize the Brillouin zone:

``` python linenums="1"
lat.plot()
plt.show()

lat.plot_brillouin_zone()
plt.show()
```

<div>
  <figure>
    <img src="../../assets/images/getting_started/lattice.png" width="300" style="display: inline-block;"/>
    <img src="../../assets/images/getting_started/brillouin.png" width="300" style="display: inline-block;"/>
    <figcaption>The visualization of the lattice and its Brillouin zone.</figcaption>
  </figure>
</div>

!!! Example "Examples"

    For a crystal with two atoms per unit cell, look in the [Examples] section.
    For other examples and pre-defined lattices consult the [Pybinding] documentation.


## Using Pybinding's solver
[Pybinding] as build-in solvers for

* [LAPACK] *(exact diagonalization)* and
* [ARPACK] (targeted diagonalization of sparse matrices).

### Making a `#!python pb.Model`
The [`#!python pb.Model`][model]-class contains all the information of the structure we want to calculate.
This structure can be larger than the unit-cell (*stored in the [`#!python pb.Lattice`][lattice]-class).
We will double the unit cell in both direction in the [`#!python pb.Model`][model] and add periodic bounadry conditions:
``` python linenums="1"
model = pb.Model(
    lat,  # pb.Lattice, use the previously defined unit-cell
    pb.primitive(2, 2),  # double the unit-cell in both directions
    pb.translational_symmetry(a1=2, a2=2)  # periodic boundary conditions with period '2'
)
```
We can visualise this [`#!python pb.Model`][model] with
``` python linenums="1"
model.plot()
plt.show()
```
<div>
  <figure>
    <img src="../../assets/images/getting_started/model.png" width="300" style="display: inline-block;"/>
    <figcaption>The visualization of the model.</figcaption>
  </figure>
</div>

### Make a `#!python pb.Solver`
The [`#!python pb.Solver`][solver]-class takes a [`#!python pb.Model`][model]-class as input and prepares the system
to do numerical calculation. We will use the [LAPACK]-solver:
``` python linenums="1"
solver = pb.solver.lapack(
    model  # pb.Model, use the previously defined system
)
```

### Band structure
As an example, the band structure is calculated using the [`#!python pb.Solver`][solver] defined above.

First, we must define the path in the reciprocal space to follow. Using the [`#!python pb.Lattice`][solver] build-in
method, the high symmetry point for the corners of a path can be found easily:
``` python linenums="1"
bz = lat.brillouin_zone()
gamma = np.array([0, 0]) 
x = (bz[1] + bz[2]) / 2
s = bz[2]
```

 Then, we just pass these corners to the [`#!python pb.Solver`][solver] and visualize the result
``` python linenums="1"
bands = solver.calc_bands(gamma, x, s, gamma, step=0.01)
bands.plot(point_labels=[r"$\Gamma$", "X", "S", r"$\Gamma$"])
plt.show()

lat.plot_brillouin_zone(decorate=False)
bands.k_path.plot(point_labels=[r"$\Gamma$", "X", "S", r"$\Gamma$"])
plt.show()
```
<div>
  <figure>
    <img src="../../assets/images/getting_started/bands.png" width="300" style="display: inline-block;"/>
    <img src="../../assets/images/getting_started/k_path.png" width="300" style="display: inline-block;"/>
    <figcaption>The visualization of the band structure and its path in the reciprocal space.</figcaption>
  </figure>
</div>

For more info about [Pybinding]'s capabilities, look at its [tutorial][tutorial-pb] or [API guide][api-pb].

!!! Example "Summary of the code from this section"

    ``` python linenums="1"
    import kite
    import pybinding as pb
    import numpy as np
    import matplotlib.pyplot as plt
    
    a1 = np.array([1, 0]) # [nm] define the first lattice vector
    a2 = np.array([0, 1]) # [nm] define the second lattice vector
    
    lat = pb.Lattice(a1=a1, a2=a2) # define a lattice object
    
    
    onsite = 0 # onsite potential
    lat.add_sublattices(
        # make a lattice site (sublattice) with a tuple
        # (name, position, and onsite potential)
        ('A', [0, 0], onsite)
    )

    lat.add_hoppings(
        # make an hopping between lattice site with a tuple
        # (relative unit cell index, site from, site to, hopping energy)
    ([1, 0], 'A', 'A', - 1 ),
    ([0, 1], 'A', 'A', - 1 )
    )
    
    model = pb.Model(
        lat,  # pb.Lattice, use the previously defined unit-cell
        pb.primitive(2, 2),  # double the unit-cell in both directions
        pb.translational_symmetry(a1=2, a2=2)  # periodic boundary conditions with period '2'
    )
    
    solver = pb.solver.lapack(
        model  # pb.Model, use the previously defined system
    )
    
    bz = lat.brillouin_zone()
    gamma = np.array([0, 0]) 
    x = (bz[1] + bz[2]) / 2
    s = bz[2]
    
    bands = solver.calc_bands(gamma, x, s, gamma, step=0.01)
    bands.plot(point_labels=[r"$\Gamma$", "X", "S", r"$\Gamma$"])
    plt.show()

    lat.plot_brillouin_zone(decorate=False)
    bands.k_path.plot(point_labels=[r"$\Gamma$", "X", "S", r"$\Gamma$"])
    plt.show()
    ```

[unitcell]: #unit-cell
[onsite]: #lattice-sites
[hopping]: #adding-hoppings
[hdf5]: https://www.hdfgroup.org
[Pybinding]: https://docs.pybinding.site/en/stable
[lattice]: https://docs.pybinding.site/en/stable/_api/pybinding.Lattice.html
[model]: https://docs.pybinding.site/en/stable/_api/pybinding.Model.html
[solver]: https://docs.pybinding.site/en/stable/_api/pybinding.solver.html#module-pybinding.solver
[ARPACK]: https://docs.pybinding.site/en/stable/_api/pybinding.solver.html#pybinding.solver.arpack
[LAPACK]: https://docs.pybinding.site/en/stable/_api/pybinding.solver.html#pybinding.solver.lapack
[tutorial-pb]: https://docs.pybinding.site/en/stable/tutorial/index.html
[api-pb]: https://docs.pybinding.site/en/stable/api.html

[tightbinding]: ../documentation/tight_binding.md
[Examples]: examples/graphene.md

[kitepython]: ../api/kite.md