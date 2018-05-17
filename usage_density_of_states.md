# Basic usage: density of states

 Our first example, a density of state calculation of a clean system,  introduces a serie of concepts that are in the core of **kite ** and can be used in the calculation of other quantities. After exporting the lattice and the disorder of choice, covered in the [tutorial for lattice importing](later), one is ready to perform some calculations with **kite**. Our first python script can be separated in the following steps:

1. Define the number of decomposition parts in each direction of matrix. This divides the lattice into various sections that are computed, each of which is calculated in parallel:

   ```python
   nx = ny = 2
   ```

   Thes decomposition allows a great speed up in of the calculation that scales with the the number of decomposed parts. We recommend its usage. However, the product of **nx** and **ny** cannot exceed the number of cores available in your computer. One must also notice that it is not efficient to decompose small systems with lateral sizes smaller than 128 unit cells of a common lattice.

2. Define the total number of units cells in each direction:

   ```python
   lx = 256
   ly = 256
   ```

   The lateral size of the decomposed parts will be given by **lx/nx** and **ly/ny** that need to be integer numbers.

3. Define an object called **configuration** that carriers several informations

   1.  The number of decomposition parts **[nx, ny]**
   2. The lateral sizes  **[lx, ly]** given in terms of unit cells
   3.  The boundary conditions: **True** for periodic boundary conditions, and **False** for open boundary conditions in each direction,
   4.  For optimisation purposes, **kite** only considers and stores complex data with the setting  **is_complex=True**. **False** indicates real values.
   5. For optimisation purposes, **kite** also allows the used to define the precision of the calculation. Use  **0** for float, **1** for double, and **2** for long double.
   6. You should define the energy scale of your hamiltonian, once Chebyshev expantions needs renomalized hamiltonians where the energy spectrum is bounded [-1,1]. If you need more details about this point, refere to *Resources* where we discuss the method in details.

   As a result, **configuration** is structured in the following way:

   ```python
   configuration = Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],is_complex=False, precision=1, energy_scale=energy_scale)
   
   ```

4. Finally it is time to write the calculation object that carries out the information about the quantity that it is going to be calculated. For that part, we still need to include more paramenters,  related to the Chebyshev expansion  (our examples already have optimized parameters for a normal desktop computer):

   1. **num_moments** defines the number of moments of the Chebyshev expansion. This number can be varied, dependening on the energy resolution you expect. Tipically we use **num_moments>max(lx,ly) **,  so it should scales with the size of your system. However, the optimal number depends on the lattice. The user should also avoid an excessive number of moments that exceed the desired energy resolution. Otherise the density of states will begin to converge to the discrete energy levels of the finite system. We provide an example of  convergence script to optimise this value.
   2.  **num_random** defines the number of random vectors involved in the stochastic calculation of quantities. This number also depends of the size of the system. For large systems that are self-averaged, it can be very small. We provide an example of  convergence script to optimise this value.
   3.  **num_disorder** defines the number of disorder realisations that are useful for disordered systems.
   4. **num_points** is the number of points the in energy axis that is going to be used in the output of the density of states.

   

   As a result, **calculation** is structured in the following way:

   ```python
   calculation = Calculation(configuration)
   calculation.dos(num_points=1000, num_random=10, num_disorder=1, num_moments=512)
   ```

   

   You can still add another object names modification to include magnetic field, for example. However, it is covered in another example.

   ```python
   modification = Modification(magnetic_field=False)
   ```

5. Finally, it is time to export all the settings to a hdf5 that is the input of our code:

   ```python
   export_lattice(lattice, configuration, calculation, modification, 'example_dos.h5')
   ```

   

You can find the source of this code [here](link)