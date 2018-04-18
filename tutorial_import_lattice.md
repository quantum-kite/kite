Quantum Kite is a multithreaded C++ package for efficient evaluation of
 spectral properties of Large-Scale Tight-Binding (TB) Hamiltonians.
 Under the assumptions of the basic knowledge of TB approximation, inside
 this tutorial we will present the code functionalities through different
 examples in form of inline code and gists which you can simply copy or
 download from our Github page _add_link_ and run. For more information on
 the TB method we suggest you to take a look at the section Resources _add_link_.
 If you encountered a peculiar problem, found a bug, or have any further question,
 we encourage you to contact us _add_link_, and we will try to respond
 as soon as possible.


# Code workflow

At present, code is divided into three independent levels. Starting
point is making a TB model using the Pybinding package. In the next section
we will show a couple of simple examples, and introduce basic functionalities
of Pybinding package, which will be used for making the model. More
advanced examples you can check in the section _add_link_.

Quantum Kite is shiped as a source code which can be compiled following
the instruction provided in the section Installation _add_link_. The code
consists of the transport code, and the post-processing script.

The interconnections between different part of the code is done using
the Hierarchical Data Format (HDF5 package). The model is made and exported
as a *.h5 file, which is used laser as a TB input to Quantum Kite, and
posprocessing tools (changing the name??).

As a last step, after obtaining the desired quantity, we can visualize
it the resulting data using a simple script _add_gist_.

In short, the code workflow is the following:

* Making and exporting a TB model from Pybinding.
* Running Quantum Kite with an input TB model.
* Running PostProcessingTools.
* Visualizing the data.

If you want to make this more automatic, you can use the following Bash
script (_add_gist_).

# Making and exporting a TB model from Pybinding

Before going to examples, let's see how to load Pybinding first.

## Importing the package

If all the requirements are fulfilled, Pybinding package can be imported
in the python script. In all the scripts in this tutorial, required packages
will be included with following aliases.
```
import pybinding as pb
import numpy as np
import matplotlib.pyplot as plt
```

If you want to use predefined styles for visualization, you can simply
write:

 ```
 pb.pltutils.use_style()
 ```

inside the script.

## Making the model

The most important object for making a TB model is pb.Lattice, which
carries the full information about the unit cell (position of lattice sites,
info on the type of sites (different atoms etc.), lattice vectors,
hopping parameters). These are all input parameters for the Lattice.
Additional features available based on the real-space information are
reciprocal vectors and the Brillouin zone.

Let's make a simple square lattice with a single lattice site.
_add_gist_
First, import all the packages:

```
import pybinding as pb
import numpy as np
import matplotlib.pyplot as plt
```

Defining the lattice vectors and making the lattice out of them
has following syntax:

```
a1 = np.array([1,0]) # [nm] define the first lattice vector
a2 = np.array([0, 1]) # [nm] define the second lattice vector

lat = pb.Lattice(a1=a1, a2=a2) # define a lattice object
```

Now we can add the desired lattice site:

```
lat.add_sublattices(
    # make a lattice site (sublattice) with a tuple
    # (name, position, and onsite potential)
    ('A', [0, 0], onsite[0])
)
```

and adding hoppings between the neighboring sites:

```
lat.add_hoppings(
    # make an hopping between lattice site with a tuple
    # (relative unit cell index, site from, site to, hopping energy)
    ([1, 0], 'A', 'A', - 1 / energy_scale),
    ([0, 1], 'A', 'A', - 1 / energy_scale)
)
```

Now we can plot the lattice:

```
lat.plot()
plt.show()
```

or visualize the Brillouin zone:

```
lat.plot_brillouin_zone()
plt.show()
```

Relative unit cell index ```[n, m]``` is a parameter of the unit cell
(in the notation ```n * a1 , m * a2```) to which the hopping occurs. The index
```[0, 0]``` is a reference hopping inside the unit cell, while other indexes mark
the periodic hopping.

It's important to emphasize that by adding the hopping between sites
```(i, j)``` the hopping therm ```(j, i)``` is added automatically, and it's not
allowed to add them twice. Also, it's not allowed to add hopping
```(i, i)``` inside the cell ```[0,0]``` because these terms are actually
onsite energies that can be added when adding a lattice site (sublattice).

Now let's try to make just a slightly advanced example, a graphene lattice:

_add_gist_
https://gist.github.com/MAndelkovic/a5f085ce48b5d28de68b03b08008b57f

## Exporting the model

After making the Lattice object, we can output both the model
and information about the quantities that we want to calculate.
For this, we need an additional functionalities that can be imported:

```
import export_lattice as ex
```

The following script displays how we can export the lattice:

https://gist.github.com/MAndelkovic/633949a568f5c7842381be265b49c02a

## TODO: Tutorial about the types of disorder and more on exporting the lattice...
