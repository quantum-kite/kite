# Disorder

Purpose of this tutorial is the introduction of different disorder realizations that can be added to the KITE Tight-Binding model. Due to the specific nature of the low memory KPM expansion implemented in KITE, full Hamiltonian matrix is never made. Therefore, although possible, it's not efficient to modify each term in a full Hamiltonian matrix. It is more suitable to for example define a rather small disorder pattern which could appear randomly.
<!-- Due to the specific nature of the low memory KPM expansion implemented in KITE, all of the realizations are randomly distributed through the domains of the simulated sample. This happens inside the C++ code and the selected pattern in the interface is randomly generated. -->

<!-- After defining a lattice (procedure is explained in the import lattice tutorial), we can add disorder to our system.  -->
Usually, disorder can be modeled either as a modification of onsite potential appearing on the lattice sites, or as a combination of onsite potential and bond disorder. Hence, KITE allows it's user to select between the two types by choosing between predefined classes in the python interface. At the present point, two specific available disorder classes are:
 * Disorder - distributed onsite disorder,
 * StructuralDisorder - generic structural disorder, the combination of onsite potential and bond disorder.

# Distributed onsite disorder

```Disorder``` realization adds randomly generated onsite terms at the sites of a desired sublattice based on a certain distribution:

* Gaussian,
* Uniform,
* Deterministic.

Beside the type of distribution, we can select a sublattice type in which the disorder will appear, and mean value and the standard deviation (CHANGE TO POTENTIAL DEPTH!!!) of the selected distribution. The procedure after making a ```lattice``` is the following:
```python
from export_lattice import Disorder # importing the Disorder class

disorder = Disorder(lattice) # define an object based on the lattice
disorder.add_disorder('A', 'Gaussian', 0.1, 0.1) # add Gaussian distributed disorder at all sites of a selected sublattice
```
In a single object it is possible to select multiple sublattices, having different disorder distributions:
```python
disorder.add_disorder('A', 'Gaussian', 0.1, 0.1)
disorder.add_disorder('B', 'Uniform', 0.2, 0.1)
disorder.add_disorder('C', 'Deterministic', 0.1)
```
After defining the desired disorder, it can be added to the configuration file as an additional parameter in the ```export_lattice``` function:
```python
export_lattice(..., disorder=disorder)
```
A full example of adding a sublattice distributed disorder and selecting a density of states before exporting a configuration file can be seen bellow:
https://gist.github.com/MAndelkovic/2a5004f5736b7d8d3d97696551d66a9a

# System disorder
```StructuralDisorder``` class adds a possibility to select between two different structural disorder types, vacancy, distributed with a certain concentration from sites in a selected sublattice, and a more generic structural disorder which is a combination of onsite terms and bond disorder (also distributed throughout the structure with a certain concentration).

## Vacancy disorder
The vacant site distribution can be selected from a single sublattice with a concentration defined in a parent object:
```python
from export_lattice import StructuralDisorder # importing the StructuralDisorder class

struc_disorder = StructuralDisorder(lattice, concentration=0.2) # define an object based on the lattice with a certain concentration
struc_disorder.add_vacancy('B') # add a vacancy to a selected sublattice with previously chosen concentration

```

## Structural disorder

The following example shows a definition of currently most generic type of disorder, which includes both onsite disorder terms and bond modifications. This type of disorder can be added as an object of a class ```StructuralDisorder```. The procedure for adding is the same as adding a hopping term to the Pybinding lattice object, with a single difference that the bond disorded is not bounded to the hopping term starting from the [0, 0] unit cell.

For the sake of clarity, let's first define sublattices that will compose the disorder. In this case we're not bounded to a single unit cell:
```python
node0 = [[+0, +0], 'A'] # define a node in a unit cell [i, j] selecting a single sublattice
node1 = [[+0, +0], 'B']
node2 = [[+1, +0], 'A']
node3 = [[+0, +1], 'B']
node4 = [[+0, +1], 'A']
node5 = [[-1, +1], 'B']
```

After the definition of a parent ```StructuralDisorder``` object, we can select the desired pattern:
```python
from export_lattice import StructuralDisorder # importing the StructuralDisorder class

struc_disorder = StructuralDisorder(lattice, concentration=0.2) # define an object based on the lattice with a certain concentration

struc_disorder.add_structural_disorder(
    # add bond disorder in the form [from unit cell], 'sublattice_from', [to_unit_cell], 'sublattice_to', value:
    (*node0, *node1, 0.5),
    (*node1, *node2, 0.1),
    (*node2, *node3, 0.5),
    (*node3, *node4, 0.3),
    (*node4, *node5, 0.4),
    (*node5, *node0, 0.8),
    # in this way we can add onsite disorder in the form [unit cell], 'sublattice', value
    ([+0, +0], 'B', 0.1)
)
# It is possible to add multiple different disorder type which should be forwarded to the export_lattice function
# as a list.
another_struc_disorder = StructuralDisorder(lat, concentration=0.6)
another_struc_disorder.add_structural_disorder(
    (*node0, *node1, 0.05),
    (*node4, *node5, 0.4),
    (*node5, *node0, 0.02),
    ([+0, +0], 'A', 0.3)
)
```
Before exporting, it's possible to define multiple disorder realizations which will be superimposed to the clean system.
The following scripts shows a minimal example of configuring above discussed disorder:
https://gist.github.com/MAndelkovic/366bcb493915c4fc0ae601c55478544f

__To show the flexibility of these predefined disorder types, we'll model an effect of functionalized graphene.__
