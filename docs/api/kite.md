The KITE package for pre-processing is split up in various subclasses and contains several functions:

* [*class* `#!python kite.StructuralDisorder(lattice, concentration=0, position=None)` - Add disorder to the lattice.][structural_disorder]
* [*class* `#!python kite.Disorder(lattice)` - Add Guassian disorder to the lattice.][disorder]
* [*class* `#!python kite.Modification(magnetic_field=None, flux=None)` - Add a magnetic field to the lattice.][modification]
* [*class* `#!python kite.Configuration([...])` - Define the basic parameters used in the calculation.][configuration]
* [*class* `#!python kite.Calculation(configuration=None)` - Describe the required target functions.][calculation]
* [*function make_pybinding_model*][make_pybinding_model]
* [*function estimate_bounds*][estimate_bounds]
* [*function config_system*][config_system]
* [*warning LoudDeprecationWarning*][loud_deprecation_warning]

## StructuralDisorder

!!! declaration-class "*class* `#!python kite.StructuralDisorder(lattice, concentration=0, position=None)`"
    
     
:   Class that introduces Structural Disorder into the initially built [`#!python pb.Lattice`][lattice].

:   **Parameters**
    : <span id='structuraldisorder-lattice'>`#!python lattice`: *[`#!python pb.Lattice`][lattice]*</span>
        : The lattice used to build the structural disorder.
    : <span id='structuraldisorder-concentration'>`#!python concentration`: *`#!python float`*</span>
        : Concentration of disorder *(can only be different from `#!python 0` if `#!python position=None`)*.
    : <span id='structuraldisorder-position'>`#!python position`: *`#!python array_like`*</span>
        : Exact position of disorder *(can only be different from `#!python None` if `#!python concentration=0`)*.

:   **Methods**
    :   | Method                                                                                                         | Description                                                          |
        |----------------------------------------------------------------------------------------------------------------| -------------------------------------------------------------------- |
        | [`#!python add_vacancy(*disorder)`][structuraldisorder-add_vacancy]                                            | Add vacancy disorder to the lattice.                                 |
        | [`#!python add_structural_disorder(*disorder)`][structuraldisorder-add_structural_disorder]                    | Add structural disorder to the lattice.                              |
        | [`#!python add_local_vacancy_disorder(relative_index, sub)`][structuraldisorder-add_local_vacancy_disorder]    | Internal function to add one vacancy disorder to chosen position.    |
        | [`#!python add_local_bond_disorder(relative_index_from, [, ...])`][structuraldisorder-add_local_bond_disorder] | Internal function to add one bond disorder between chosen positions. |
        | [`#!python add_local_onsite_disorder(relative_index, [, ...])`][structuraldisorder-add_local_onsite_disorder]  | Internal function to add one onsite disorder to chosen position.     |
        | [`#!python map_the_orbital(orb, nodes_map)`][structuraldisorder-map_the_orbital]                               | Internal function to map the orbitals to the NodesMap.               |

    :   !!! declaration-function "<span id='structuraldisorder-add_vacancy'>*function* `#!python add_vacancy(*disorder)`"
            
            
        :   Add vacancy disorder to the lattice.
            
             **Parameters**

             :  | Parameter                                                                  | Description                                                                                                       |
                |----------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------|
                | `#!python *disorder`:*`#!python str` or `#!python tuple(array_like, str)`* | Vacancy disorder, in the form of *`#!python sublattice_name`* or *`#!python ([relatice_index], sublattice_name)`* |
    
    :   !!! declaration-function "<span id='structuraldisorder-add_structural_disorder'>*function*`#!python add_structural_disorder(*disorder)`</span>"
            
            
        :   Add structural disorder to the lattice.
            
             **Parameters**

             :  | Parameter                                                                                                                                                                                                                                       | Description                                                                                                                                                                                                            |
                |-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
                | `#!python *disorder`:*`#!python tuple(array_like, str, float)` or `#!python tuple(array_like, str, array_like, str, array_like)` or `#!python tuple(array_like, str, array_like, str, float)` or `#!python tuple(array_like, str, array_like)`* | Vacancy disorder, in the form of *`#!python ([relatice_index], sublattice_name, onsite_energy)`* or *`#!python ([relatice_index_from], sublattice_name_from, [relatice_index_to], sublattice_name_to, onsite_energy)`* |
    
    

    :   !!! declaration-function "<span id='structuraldisorder-add_local_vacancy_disorder'>*function*`#!python add_local_vacancy_disorder(relative_index, sub)`</span>"
            
            
        :   Internal function to add one vacancy disorder to chosen position.
            
             **Parameters**

             :  | Parameter                                         | Description                                                 |
                |---------------------------------------------------|-------------------------------------------------------------|
                | `#!python relative_index`:*`#!python array_like`* | Relative index of the position to change the onsite energy. |
                | `#!python sub`:*`#!python str`*                   | Name of the sublattice to change the onsite energy.         |

    :   !!! declaration-function "<span id='structuraldisorder-add_local_bond_disorder'>*function* `#!python add_local_bond_disorder(relative_index_from, from_sub, relative_index_to, to_sub, hoppings)`</span>"
            
            
        :   Internal function to add one bond disorder between chosen positions.
            
             **Parameters**

             :  | Parameter                                                       | Description                                                                                                                        |
                |-----------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------|
                | `#!python relative_index_from`:*`#!python array_like`*          | Relative index of the position from wich the bond connects to change the onsite energy.                                            |
                | `#!python from_sub`:*`#!python str`*                            | Name of the sublattice from wich the bond connects to change the onsite energy.                                                    |
                | `#!python relative_index_to`:*`#!python array_like`*            | Relative index of the position to wich the bond connects to change the onsite energy.                                              |
                | `#!python to_sub`:*`#!python str`*                              | Name of the sublattice to wich the bond connects to change the onsite energy.                                                      |
                | `#!python hoppings`:*`#!python float` or `#!python array_like`* | The hopping energy between the different sublattices at the given positions, with the right shape to connect between the orbitals. |
    
    :   !!! declaration-function "<span id='structuraldisorder-add_local_onsite_disorder'>*function* `#!python add_local_onsite_disorder(relative_index, sub, value)`</span>"
            
            
        :   AInternal function to add one onsite disorder to chosen position.
            
             **Parameters**

            :  | Parameter                                                     | Description                                                                                    |
               |--------------------------------------------------------------|------------------------------------------------------------------------------------------------|
               | `#!python relative_index`:*`#!python array_like`*            | Relative index of the position to change the onsite energy.                                    |
               | `#!python sub`:*`#!python str`*                              | Name of the sublattice to change the onsite energy.                                            |
               | `#!python value`:*`#!python float` or `#!python array_like`* | The onsite energy of sublattice at the given positions, with the right shape for the orbitals. |
    
    :   !!! declaration-function "<span id='structuraldisorder-map_the_orbital'>*function* `#!python map_the_orbital(orb, nodes_map)`</span>"
            
            
        :   Internal function to map the orbitals to the NodesMap.
            
             **Parameters**

            :   | Parameter                              | Description                                                   |
                |----------------------------------------|---------------------------------------------------------------|
                | `#!python orb`:*`#!python str`*        | Name of the sublattice to give a unique value.                |
                | `#!python nodes_map`:*`#!python dict`* | The object that stores the unique values for the sublattices. |

## Disorder

!!! declaration-class "*class* `#!python kite.Disorder(lattice)`"
    
     
:   Class that introduces Disorder into the initially built [`#!python pb.Lattice`][lattice].

    The informations about the disorder are the *type*, *mean value*, and *standard deviation*.
    The function that you can use in the bulding of the [`#!python pb.Lattice`][lattice] is `#!python add_disorder()`.
    The class method takes care of the shape of the disorder chosen
    (it needs to be same as the number of orbitals at a given atom),
    and takes care of the conversion to the c++ orbital-only format.

:   **Parameters**
    : <span id='disorder-lattice'>`#!python lattice`: *[`#!python pb.Lattice`][lattice]*</span>
        : The lattice used to build the disorder.

:   **Methods**
    :   | Method                                                                                | Description                                              |
        | ------------------------------------------------------------------------------------- | -------------------------------------------------------- |
        | [`#!python add_disorder(sublattice [, ...])`][disorder-add_disorder]                  | Add the disorder to the lattice.                         |
        | [`#!python add_local_disorder(sublattice_name [, ...])`][disorder-add_local_disorder] | Internal function to add the disorder to the positions.  |
    
    :   !!! declaration-function "<span id='disorder-add_disorder'>*function* `#!python add_disorder(sublattice, dis_type, mean_value, standard_deviation=0.)`</span>"
            
            
        :   Add the disorder to the lattice.
            
             **Parameters**

             :  | Parameter                                                                  | Description                                                                                                                                                                                              |
                |----------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
                | `#!python sublattice`:*`#!python str` or `#!python list(str)`*             | Name of the sublattice to give a unique value.                                                                                                                                                           |
                | `#!python dis_type`:*`#!python int` or `#!python list(int)`*               | The type of disorder to apply, possible values are `#!python "Gaussian"`, `#!python "Uniform"`, `#!python "Deterministic"`,  `#!python "gaussian"`,  `#!python "uniform"` or `#!python "deterministic"`. |
                | `#!python dis_mean_valuetype`:*`#!python float` or `#!python list(float)`* | Mean value of the deformation.                                                                                                                                                                           |
                | `#!python standard_deviation`:*`#!python float` or `#!python list(float)`* | Standard deviation of the deformation.                                                                                                                                                                   |

    :   !!! declaration-function "<span id='disorder-add_local_disorder'>*function*`#!python add_local_disorder(sublattice_name, dis_type, mean_value, standard_deviation)`</span>"
            
            
        :   Internal function to add the disorder to the positions.
            
             **Parameters**

             :  | Parameter                                              | Description                                                                                                                                                                                              |
                |--------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
                | `#!python sublattice`:*`#!python list(str)`*           | Name of the sublattices to give a unique value.                                                                                                                                                          |
                | `#!python dis_type`:*`#!python list(int)`*             | The type of disorder to apply, possible values are `#!python "Gaussian"`, `#!python "Uniform"`, `#!python "Deterministic"`,  `#!python "gaussian"`,  `#!python "uniform"` or `#!python "deterministic"`. |
                | `#!python dis_mean_valuetype`:*`#!python list(float)`* | Mean value of the deformation.                                                                                                                                                                           |
                | `#!python standard_deviation`:*`#!python list(float)`* | Standard deviation of the deformation.                                                                                                                                                                   |

## Modification
!!! declaration-class "*class* `#!python kite.Modification(magnetic_field=None, flux=None)`"


:   Class that modifies the initially built [`#!python pb.Lattice`][lattice] with a [magnetic field][magnetic-field].

:   **Parameters**
    : <span id="modification-par-magnetic_field">`#!python magnetic_field`: *`#!python float`*</span>
        : Add the magnetic field to the lattice. The field will point along the second primitive lattice vector of the lattice. The magnetic field is in units of $Tesla$, if the [`#!python pb.Lattice`][lattice] is in units of $nm$. The magnetic field is rounded down to the nearest flux quantum.
    : <span id="modification-par-flux">`#!python flux`: *`#!python float`*</span>
        : Add the magnetic flux to the lattice.

:   **Attributes**
    :   | Attribute                                                                                      | Description                                                                                                                                      |
        | ---------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
        | <span id="modification-atr-magnetic_field">`#!python magnetic_field`:*`#!python float`*</span> | The added magnetic field to the lattice.                                                                                                         |
        | <span id="modification-atr-flux">`#!python flux`:*`#!python float`*</span>                     | The added magnetic flux to the lattice. *This is **not** the exact value used in the calculation, but the value added using the parameter above. |

## Configuration
!!! declaration-class "*class* `#!python kite.Configuration(divisions=(1, 1, 1), length=(1, 1, 1), boundaries=('open', 'open', 'open'), is_complex=False, precision=1, spectrum_range=None, angles=(0,0,0), custom_local=False, custom_local_print=False)`"
    
     
:   Define the basic parameters used in the calculation

:   **Parameters**
    : <span id="configuration-divisions">`#!python divisions`: *`#!python int` or `#!python tuple(int, int)` or `#!python tuple(int, int, int)`*</span>
        : Number of decomposition parts of the system.
    : <span id="configuration-length">`#!python length`: *`#!python int` or `#!python tuple(int, int)` or `#!python tuple(int, int, int)`*</span>
        : Number of unit cells in each direction.
    : <span id="configuration-boundaries">`#!python boundaries`: *`#!python str` or `#!python tuple(str, str)` or `#!python tuple(str, str, str)`*</span>
        : Periodic boundary conditions each direction. Possible values are `#!python "periodic"`, `#!python "open"`,
          `#!python "twisted"`*(this option needs the extra argument `#!python angles` and `#!python "random"` 
    : <span id="configuration-is_complex">`#!python is_complex`: *`#!python bool`*</span>
        : Boolean that reflects whether the type of Hamiltonian is complex or not.
    : <span id="configuration-precision">`#!python precision`: *`#!python int`*</span>
        : Integer which defines the precision of the number used in the calculation,
          `#!python float` (`#!python 0`), `#!python double` (`#!python 1`), `#!python long double` (`#!python 2`).
    : <span id="configuration-spectrum_range">`#!python spectrum_range`: *`#!python tuple(float, float)`*</span>
        : Energy scale which defines the scaling factor of all the energy related parameters.
          The scaling is done automatically in the background after this definition.
          If the term is not specified, a rough estimate of the bounds is found.
            
            !!! Warning

                Automatic scaling can lead to segmentation-errors due to an error in [pybinding].

    : <span id="configuration-angles">`#!python angles`: *`#!python float` or `#!python tuple(float, float)` or `#!python tuple(float, float, float)`*</span>
        : The angles used for the twisted boundary conditions when `#!python boundary="twist"` is selected.
          The values of `#!python angle`must be in the interval $[0, M \cdot 2 \pi]$.
    : <span id="configuration-custom_local">`#!python custom_local`: *`#!python bool`*</span>
        : Boolean that reflects whether the calculation should use the user-defined local potential.
    : <span id="configuration-custom_local_print">`#!python custom_local_print`: *`#!python bool`*</span>
        : Boolean that reflects whether the calculation should use output the values for the local potential for the various sublattices of the [`#!python pb.Lattice`][lattice].

:   **Attributes**

    :   | Attribute                                                                                                                                                 | Description                                                                                                                                                                                                                                                                                                                                                   |
        |-----------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
        | <span id="configuration-energy_scale">`#!python energy_scale`:*`#!python float`*</span>                                                                   | Returns the energy scale of the hopping parameters.                                                                                                                                                                                                                                                                                                           |
        | <span id="configuration-energy_shift">`#!python energy_shift`:*`#!python float`*</span>                                                                   | Returns the energy shift of the hopping parameters around which the spectrum is centered.                                                                                                                                                                                                                                                                     |
        | <span id="configuration-comp">`#!python comp`:*`#!python int`*</span>                                                                                     | Returns `#!python 0` if hamiltonian is real and `#!python 1` elsewise.                                                                                                                                                                                                                                                                                        |
        | <span id="configuration-prec">`#!python prec`:*`#!python int`*</span>                                                                                     | Returns `#!python 0`, `#!python 1`, `#!python 2` if precision if `#!python float`, `#!python double`, and `#!python long double` respectively.                                                                                                                                                                                                                |
        | <span id="configuration-div">`#!python div`:*`#!python int`*</span>                                                                                       | Returns the number of decomposed elements of matrix in $x$, $y$ and/or $z$ direction. Their product gives the total number of threads spawn.                                                                                                                                                                                                                  |
        | <span id="configuration-bound">`#!python bound`:*`#!python tuple(array_like, array_like)`*</span>                                                         | Returns the boundary conditions in each direction, the first argument describes the boundary conditions for the various dimensions with `#!python 0` - open boundary condtions, `#!python 1` - periodic or twisted boundary conditions, `#!python 2` - random boundary conditions, the second gives the angle used if `#!python boundary="twist"` was chosen. |
        | <span id="configuration-leng">`#!python leng`:*`#!python array_like`*</span>                                                                              | Return the number of unit cell repetitions in each direction.                                                                                                                                                                                                                                                                                                 |
        | <span id="configuration-type">`#!python type`:*`#!python np.float32` or `#!python np.float64` or `#!python np.float128` or `#!python np.float256`*</span> | Return the type of the Hamiltonian complex or real, and float, double or long double.                                                                                                                                                                                                                                                                         |
        | <span id="configuration-custom_pot">`#!python custom_pot`:*`#!python bool`*</span>                                                                        | Return custom potential flag.                                                                                                                                                                                                                                                                                                                                 |
        | <span id="configuration-print_custom_pot">`#!python print_custom_pot`:*`#!python bool`*</span>                                                            | Return print custom potential flag.                                                                                                                                                                                                                                                                                                                           |


:   **Methods**
    :   | Method                                         | Description                                                                     |
        | ----------------------------------------------- | ------------------------------------------------------------------------------- |
        | [`#!python set_type()`][configuration-set_type] | Internal function to determine the precision to be used during the calculation. |

    :   !!! declaration-function "<span id="configuration-set_type">*function* `#!python set_type()`</span>"
            
            
        :   Internal function to determine the precision to be used during the calculation.

## Calculation
!!! declaration-class "*class* `#!python kite.Calculation(configuration=None)`"
    
     
:   Class that contains the required target functions to calculate in later steps.

:   **Parameters**
    : <span id="calculation-configuration">`#!python configuration`: *[`#!python kite.Configuration`][configuration]*</span>
        : A [`#!python kite.Configuration`][configuration] that contains the settings for the calculation.

:   **Attributes**
    :   | Attribute                                                                                                                        | Description                                                                                                                 |
        |----------------------------------------------------------------------------------------------------------------------------------| --------------------------------------------------------------------------------------------------------------------------- |
        | <span id="calculation-get_dos">`#!python get_dos`:*`#!python dict`*</span>                                                       | Returns the requested DOS functions.                                                                                        |
        | <span id="calculation-get_ldos">`#!python get_ldos`:*`#!python dict`*</span>                                                     | Returns the requested LDOS functions.                                                                                       |
        | <span id="calculation-get_arpes">`#!python get_arpes`:*`#!python dict`*</span>                                                   | Returns the requested ARPES functions.                                                                                      |
        | <span id="calculation-get_gaussian_wave_packet">`#!python get_gaussian_wave_packet`:*`#!python dict`*</span>                     | Returns the requested wave packet time evolution function, with a gaussian wavepacket mutiplied with different plane waves. |
        | <span id="calculation-get_conductivity_dc">`#!python get_conductivity_dc`:*`#!python dict`*</span>                               | Returns the requested DC conductivity functions.                                                                            |
        | <span id="calculation-get_conductivity_optical">`#!python get_conductivity_optical`:*`#!python dict`*</span>                     | Returns the requested optical conductivity functions.                                                                       |
        | <span id="calculation-get_conductivity_optical_nonlinear">`#!python get_conductivity_optical_nonlinear`:*`#!python dict`*</span> | Returns the requested nonlinear optical conductivity functions.                                                             |
        | <span id="calculation-get_singleshot_conductivity_dc">`#!python get_singleshot_conductivity_dc`:*`#!python dict`*</span>         | Returns the requested singleshot DC conductivity functions.                                                                 |


:   **Methods**
    :   | Method                                                                                         | Description                                                                 |
        |------------------------------------------------------------------------------------------------| --------------------------------------------------------------------------- |
        | [`#!python dos(num_points, num_moments, [, ...]`)][calculation-dos]                            | Calculate the density of states as a function of energy.                    |
        | [`#!python ldos(energy, num_moments, [, ...])`][calculation-ldos]                              | Calculate the local density of states as a function of energy.              |
        | [`#!python arpes(k_vector, weight [, ...])`][calculation-arpes]                                | Calculate the spectral contribution for given k-points and weights.         |
        | [`#!python gaussian_wave_packet(num_points [, ...])`][calculation-gaussian_wave_packet]        | Calculate the time evolution function of a wave packet.                     |
        | [`#!python conductivity_dc(direction, [, ...])`][calculation-conductivity_dc]                  | Calculate the DC conductivity for a given direction.                        |
        | [`#!python conductivity_optical(direction, [, ...])`][calculation-conductivity_optical]        | Calculate optical conductivity for a given direction.                       |
        | [`#!python conductivity_optical_nonlinear([...])`][calculation-conductivity_optical_nonlinear] | Calculate nonlinear optical conductivity for a given direction.             |
        | [`#!python singleshot_conductivity_dc(energy, [...])`][calculation-singleshot_conductivity_dc] | Calculate the DC conductivity using KITEx for a given direction and energy. |

    :   !!! declaration-function "<span id="calculation-dos">*function* `#!python dos(num_points, num_moments, num_random, num_disorder=1)`</span>"
            
            
        :   Calculate the density of states as a function of energy.
            
            **Parameters**

            :   | Parameter                                | Description                                                                     |
                |------------------------------------------|---------------------------------------------------------------------------------|
                | `#!python num_points`:*`#!python int`*   | Number of energy point inside the spectrum at which the DOS will be calculated. |
                | `#!python num_moments`:*`#!python int`*  | Number of polynomials in the Chebyshev expansion.                               |
                | `#!python num_random`:*`#!python int`*   | Number of random vectors to use for the stochastic evaluation of trace.         |
                | `#!python num_disorder`:*`#!python int`* | Number of different disorder realisations.                                      |

    
    :   !!! declaration-function "<span id="calculation-ldos">*function*`#!python ldos(energy, num_moments, position, sublattice, num_disorder=1)`</span>"
            
            
        :   Calculate the local density of states as a function of energy. 
            
             **Parameters**

            :   | Parameter                                                   | Description                                                        |
                |-------------------------------------------------------------|--------------------------------------------------------------------|
                | `#!python energy`:*`#!python array_like`*                   | List of energy points at which the LDOS will be calculated.        |
                | `#!python num_moments`:*`#!python int`*                     | Number of polynomials in the Chebyshev expansion.                  |
                | `#!python position`:*`#!python int`*                        | Relative index of the unit cell where the LDOS will be calculated. |
                | `#!python sublattice`:*`#!python list`*                     | Name of the sublattice at which the LDOS will be calculated.       |
                | `#!python num_disorder`:*`#!python str` or `#!python list`* | Number of different disorder realisations.                         |
    
    :   !!! declaration-function "<span id="calculation-arpes">*function*`#!python arpes(k_vector, weight, num_moments, num_disorder=1)`</span>"
            
            
        :   Calculate the spectral contribution for given k-points and weights.
            
            **Parameters**

            :   | Parameter                                                   | Description                                                                                                   |
                |-------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------|
                | `#!python k_vector`:*`#!python array_like`*                 | List of K points with respect to reciprocal vectors b0 and b1 at which the band structure will be calculated. |
                | `#!python weight`:*`#!python array_like`*                   | List of orbital weights used for ARPES.                                                                       |
                | `#!python num_moments`:*`#!python int`*                     | Number of polynomials in the Chebyshev expansion.                                                             |
                | `#!python num_disorder`:*`#!python int`*                    | Number of different disorder realisations.                                                                    |

    
    :   !!! declaration-function "<span id="calculation-gaussian_wave_packet">*function*`#!python gaussian_wave_packet(num_points, num_moments, timestep, k_vector, spinor, width, mean_value, num_disorder=1, probing_point=0)`</span>"
            
            
        :   Calculate the time evolution function of a wave packet.
            
            **Parameters**

            :   | Parameter                                                          | Description                                                                                                     |
                |--------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------|
                | `#!python num_points`:*`#!python int`*                             | Number of time points for the time evolution.                                                                   |
                | `#!python num_moments`:*`#!python int`*                            | Number of polynomials in the Chebyshev expansion.                                                               |
                | `#!python timestep`:*`#!python float`*                             | Timestep for calculation of time evolution.                                                                     |
                | `#!python k_vector`:*`#!python array_like`*                        | Different wave vectors, components corresponding to vectors b0 and b1.                                          |
                | `#!python spinor`:*`#!python array_like`*                          | Spinors for each of the k vectors.                                                                              |
                | `#!python width`:*`#!python float`*                                | Width of the gaussian.                                                                                          |
                | `#!python mean_value`:*`#!python tuple(float, float)`*             | Mean value of the gaussian envelope.                                                                            |
                | `#!python num_disorder`:*`#!python int`*                           | Number of different disorder realisations.                                                                      |
                | `#!python probing_point`:*`#!python int` or `#!python array_like`* | Forward probing point, defined with x, y coordinate were the wavepacket will be checked at different timesteps. |

    
    :   !!! declaration-function "<span id="calculation-conductivity_dc">*function*`#!python conductivity_dc(direction, num_points, num_moments, num_random, num_disorder=1, temperature=0)`</span>"
            
            
        :   Calculate the DC conductivity for a given direction.
            
            **Parameters**

            :   | Parameter                                 | Description                                                                                                                                                                                                                                      |
                |-------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
                | `#!python direction`:*`#!python str`*     | Direction in $xyz$-coordinates along which the conductivity is calculated, supports `#!python "xx"`, `#!python "yy"`, `#!python "zz"`, `#!python "xy"`, `#!python "xz"`, `#!python "yx"`, `#!python "yz"`, `#!python "zx"`, `#!python "zy"`.     |
                | `#!python num_points`:*`#!python int`*    | Number of energy point inside the spectrum at which the DOS will be calculated.                                                                                                                                                                  |
                | `#!python num_moments`:*`#!python int`*   | Number of polynomials in the Chebyshev expansion.                                                                                                                                                                                                |
                | `#!python num_random`:*`#!python int`*    | Number of random vectors to use for the stochastic evaluation of trace.                                                                                                                                                                          |
                | `#!python num_disorder`:*`#!python int`*  | Number of different disorder realisations.                                                                                                                                                                                                       |
                | `#!python temperature`:*`#!python float`* | Value of the temperature at which we calculate the response. If $eV$ is used as unit for energy, then $k_B\cdot T$ is also in $eV$. To define the temperature in arbitraty units, specify the quantity $K_B \cdot T$, which has units of energy. |

    
    
    
    :   !!! declaration-function "<span id="calculation-conductivity_optical">*function*`#!python conductivity_optical(direction, num_points, num_moments, num_random, num_disorder=1, temperature=0)`</span>"
            
            
        :   Calculate optical conductivity for a given direction.
                        
            **Parameters**

            :   | Parameter                                 | Description                                                                                                                                                                                                                                      |
                |-------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
                | `#!python direction`:*`#!python str`*     | Direction in $xyz$-coordinates along which the conductivity is calculated, supports `#!python "xx"`, `#!python "yy"`, `#!python "zz"`, `#!python "xy"`, `#!python "xz"`, `#!python "yx"`, `#!python "yz"`, `#!python "zx"`, `#!python "zy"`.     |
                | `#!python num_points`:*`#!python int`*    | Number of energy point inside the spectrum at which the DOS will be calculated.                                                                                                                                                                  |
                | `#!python num_moments`:*`#!python int`*   | Number of polynomials in the Chebyshev expansion.                                                                                                                                                                                                |
                | `#!python num_random`:*`#!python int`*    | Number of random vectors to use for the stochastic evaluation of trace.                                                                                                                                                                          |
                | `#!python num_disorder`:*`#!python int`*  | Number of different disorder realisations.                                                                                                                                                                                                       |
                | `#!python temperature`:*`#!python float`* | Value of the temperature at which we calculate the response. If $eV$ is used as unit for energy, then $k_B\cdot T$ is also in $eV$. To define the temperature in arbitraty units, specify the quantity $K_B \cdot T$, which has units of energy. |

    
    :   !!! declaration-function "<span id="calculation-conductivity_optical_nonlinear">*function*`#!python conductivity_optical_nonlinear(direction, num_points, num_moments, num_random, num_disorder=1, temperature=0, special=0)`</span>"
            
            
        :   Calculate nonlinear optical conductivity for a given direction.
            
            **Parameters**

            :   | Parameter                                 | Description                                                                                                                                                                                                                                                                  |
                |-------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
                | `#!python direction`:*`#!python str`*     | Direction in $xyz$-coordinates along which the conductivity is calculated, supports all the combinations of the directions `#!python "x"`, `#!python "y"` and `#!python "z"` with length 3 like `#!python "xxx"`, `#!python "zzz"`, `#!python "xxy"`, `#!python "xxz"` etc.  |
                | `#!python num_points`:*`#!python int`*    | Number of energy point inside the spectrum at which the DOS will be calculated.                                                                                                                                                                                              |
                | `#!python num_moments`:*`#!python int`*   | Number of polynomials in the Chebyshev expansion.                                                                                                                                                                                                                            |
                | `#!python num_random`:*`#!python int`*    | Number of random vectors to use for the stochastic evaluation of trace.                                                                                                                                                                                                      |
                | `#!python num_disorder`:*`#!python int`*  | Number of different disorder realisations.                                                                                                                                                                                                                                   |
                | `#!python temperature`:*`#!python float`* | Value of the temperature at which we calculate the response. If $eV$ is used as unit for energy, then $k_B\cdot T$ is also in $eV$. To define the temperature in arbitraty units, specify the quantity $K_B \cdot T$, which has units of energy.                             |
                | `#!python special`:*`#!python int`*       | Optional, a parameter that can simplify the calculation for some materials.                                                                                                                                                                                                  |
    
    
    :   !!! declaration-function "<span id="calculation-singleshot_conductivity_dc">*function*`#!python singleshot_conductivity_dc(energy, direction, eta, num_moments, num_random, num_disorder=1, preserve_disorder=False)`</span>"
            
            
        :   Calculate the DC conductivity using KITEx for a given direction and energy.
            
            !!! Info "Processing the output of `#!python singleshot_conductivity_dc()`"

                `#!python singleshot_conductivity_dc()` works different from the other target-functions in that a
                single run with [KITEx][kitex] is sufficient. The results don't have to be processed by
                [KITE-tools][kitetools].
                As such, the results are already available in the [HDF5]-file. You can extract the results from the
                [HDF5]-file [as explained in the tutorial][tutorial-hdf5].

                There is also [a script][kitetools-output] that automates this process with, `#!python "output.h5"`
                the name of the [HDF5]-file processed by [KITEx][kitex]:
                
                ``` bash
                python3 process_single_shot.py output.h5
                ``` 
                
                this will result in a data-file `#!python "output.dat"`, as explained in the
                [API of KITE-tools][kitetools-output].                

            **Parameters**

            :   | Parameter                                                     | Description                                                                                                                               |
                |---------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------|
                | `#!python energy`:*`#!python array_like` or `#!python float`* | Array or a single value of energies at which `#!python singleshot_conductivity_dc` will be calculated.                                    |
                | `#!python direction`:*`#!python str`*                         | Direction in $xyz$-coordinates along which the conductivity is calculated, supports `#!python "xx"`, `#!python "yy"` and `#!python "zz"`. |
                | `#!python eta`:*`#!python int`*                               | Parameter that affects the broadening of the kernel function.                                                                             |
                | `#!python num_moments`:*`#!python int`*                       | Number of polynomials in the Chebyshev expansion.                                                                                         |
                | `#!python num_random`:*`#!python int`*                        | Number of random vectors to use for the stochastic evaluation of trace.                                                                   |
                | `#!python num_disorder`:*`#!python int`*                      | Number of different disorder realisations.                                                                                                |
                | `#!python preserve_disorder`:*`#!python bool`*                | Optional.                                                                                                                                 |

## make_pybinding_model

:   !!! declaration-function "*function* `#!python kite.make_pybinding_model(lattice, disorder=None, disorder_structural=None, shape=None)`"
            
            
    :   Build a Pybinding model with disorder used in Kite. Bond disorder or magnetic field are not currently supported.

        **Parameters**
    
        :   | Parameter                                                                                  | Description                                                                                                                                                                                                |
            |--------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
            | `#!python lattice`:*[`#!python pb.Lattice`][lattice]*                                      | Pybinding lattice object that carries the info about the unit cell vectors, unit cell cites, hopping terms and onsite energies.                                                                            |
            | `#!python disorder`:*[`#!python kite.Disorder`][disorder]*                                 | Class that introduces [`#!python kite.Disorder`][disorder] into the initially built lattice. For more info check the [`#!python kite.Disorder`][disorder] class.                                           |
            | `#!python disorder_structural`:*[`#!python kite.StructuralDisorder`][structural_disorder]* | Class that introduces [`#!python kite.StructuralDisorder`][structural_disorder] into the initially built lattice. For more info check the [`#!python kite.StructuralDisorder`][structural_disorder] class. |
            | `#!python shape`:*[`#!python pb.Shape`][shape]*                                            | Optional argument [`#!python pb.Shape`][shape].                                                                                                                                                            |

## estimate_bounds

:   !!! declaration-function "*function* `#!python kite.estimate_bounds(lattice, disorder=None, disorder_structural=None)`"
            
            
    :   Estimate the bounds for the energy, given the [`#!python pb.Lattice`][lattice] and/or [`#!python kite.Disorder`][disorder] and/or [`#!python kite.StructuralDisorder`][structural_disorder].

        **Parameters**
    
        :   | Parameter                                                                                  | Description                                                                                                                                                                                                |
            |--------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
            | `#!python lattice`:*[`#!python pb.Lattice`][lattice]*                                      | Pybinding lattice object that carries the info about the unit cell vectors, unit cell cites, hopping terms and onsite energies.                                                                            |
            | `#!python disorder`:*[`#!python kite.Disorder`][disorder]*                                 | Class that introduces [`#!python kite.Disorder`][disorder] into the initially built lattice. For more info check the [`#!python kite.Disorder`][disorder] class.                                           |
            | `#!python disorder_structural`:*[`#!python kite.StructuralDisorder`][structural_disorder]* | Class that introduces [`#!python kite.StructuralDisorder`][structural_disorder] into the initially built lattice. For more info check the [`#!python kite.StructuralDisorder`][structural_disorder] class. |

## config_system

:   !!! declaration-function "*function* `#!python kite.config_system(lattice, config, calculation, modification=None, filename="kite_config.h5", disorder=None, disorder_structural=None)`"
            
            
    :   Export the lattice and related parameters to the *.h5 file.

        **Parameters**
    
        :   | Parameter                                                                                  | Description                                                                                                                                                                                                |
            |--------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
            | `#!python lattice`:*[`#!python pb.Lattice`][lattice]*                                      | Pybinding lattice object that carries the info about the unit cell vectors, unit cell cites, hopping terms and onsite energies.                                                                            |
            | `#!python config`:*[`#!python kite.Configuration`][configuration]*                         | [`#!python kite.Configuration`][configuration] object, basic parameters defining size, precision, energy scale and number of decomposition parts in the calculation.                                       |
            | `#!python calculation`:*[`#!python kite.Calculation`][calculation]*                        | [`#!python kite.Calculation`][calculation] object that defines the requested functions for the calculation.                                                                                                |
            | `#!python modification`:*[`#!python kite.Modification`][modification]*                     | If specified [`#!python kite.Modification`][modification] object, has the magnetic field selection, either in terms of field, or in the number of flux quantum through the selected system.                |
            | `#!python filename`: *`#!python str`*                                                      | Filename for the output HDF5-file.                                                                                                                                                                         |
            | `#!python disorder`:*[`#!python kite.Disorder`][disorder]*                                 | Class that introduces [`#!python kite.Disorder`][disorder] into the initially built lattice. For more info check the [`#!python kite.Disorder`][disorder] class.                                           |
            | `#!python disorder_structural`:*[`#!python kite.StructuralDisorder`][structural_disorder]* | Class that introduces [`#!python kite.StructuralDisorder`][structural_disorder] into the initially built lattice. For more info check the [`#!python kite.StructuralDisorder`][structural_disorder] class. |

## LoudDeprecationWarning

:   Deprecationwarning.


[comment]: <> (Classes from Pybinding)
[lattice]: https://docs.pybinding.site/en/stable/_api/pybinding.Lattice.html
[pybinding]: https://docs.pybinding.site/
[shape]: https://docs.pybinding.site/en/stable/api.html#shapes

[comment]: <> (Class StructuralDisorder)
[structural_disorder]: #structuraldisorder
[comment]: <> (Class Parameters)
[structuraldisorder-lattice]: #structuraldisorder-lattice
[structuraldisorder-concentration]: #structuraldisorder-concentration
[structuraldisorder-position]: #structuraldisorder-position
[comment]: <> (Class Methods)
[structuraldisorder-add_vacancy]: #structuraldisorder-add_vacancy
[structuraldisorder-add_structural_disorder]: #structuraldisorder-add_structural_disorder
[structuraldisorder-add_local_vacancy_disorder]: #structuraldisorder-add_local_vacancy_disorder
[structuraldisorder-add_local_bond_disorder]: #structuraldisorder-add_local_bond_disorder
[structuraldisorder-add_local_onsite_disorder]: #structuraldisorder-add_local_onsite_disorder
[structuraldisorder-map_the_orbital]: #structuraldisorder-map_the_orbital

[comment]: <> (Class Disorder)
[disorder]: #disorder
[comment]: <> (Class Parameters)
[disorder-lattice]: #disorder-lattice
[comment]: <> (Class Methods)
[disorder-add_disorder]: #disorder-add_disorder
[disorder-add_local_disorder]: #disorder-add_local_disorder

[comment]: <> (Class Modification)
[modification]: #modification
[comment]: <> (Class Parameters)
[modification-par-magnetic_field]: #modification-par-magnetic_field
[modification-par-flux]: #modification-par-flux
[comment]: <> (Class Attributes)
[modification-atr-magnetic_field]: #modification-atr-magnetic_field
[modification-atr-flux]: #modification-atr-flux

[comment]: <> (Class Calculation)
[calculation]: #calculation
[comment]: <> (Class Parameters)
[calculation-configuration]: #calculation-configuration
[comment]: <> (Class Attributes)
[calculation-get_dos]: #calculation-get_dos
[calculation-get_ldos]: #calculation-get_ldos
[calculation-get_arpes]: #calculation-get_arpes
[calculation-get_gaussian_wave_packet]: #calculation-get_gaussian_wave_packet
[calculation-get_conductivity_dc]: #calculation-get_conductivity_dc
[calculation-get_conductivity_optical]: #calculation-get_conductivity_optical
[calculation-get_conductivity_optical_nonlinear]: #calculation-get_conductivity_optical_nonlinear
[calculation-get_singleshot_conductivity_dc]: #calculation-get_singleshot_conductivity_dc
[comment]: <> (Class Methods)
[calculation-dos]: #calculation-dos
[calculation-ldos]: #calculation-ldos
[calculation-arpes]: #calculation-arpes
[calculation-gaussian_wave_packet]: #calculation-gaussian_wave_packet
[calculation-conductivity_dc]: #calculation-conductivity_dc
[calculation-conductivity_optical]: #calculation-conductivity_optical
[calculation-conductivity_optical_nonlinear]: #calculation-conductivity_optical_nonlinear
[calculation-singleshot_conductivity_dc]: #calculation-singleshot_conductivity_dc

[comment]: <> (Class Configuration)
[configuration]: #configuration
[comment]: <> (Class Parameters)
[configuration-divisions]: #configuration-divisions
[configuration-length]: #configuration-length
[configuration-boundaries]: #configuration-boundaries
[configuration-is_complex]: #configuration-is_complex
[configuration-precision]: #configuration-precision
[configuration-spectrum_range]: #configuration-spectrum_range
[configuration-angles]: #configuration-angles
[configuration-custom_local]: #configuration-custom_local
[configuration-custom_local_print]: #configuration-custom_local_print
[comment]: <> (Class Attributes)
[configuration-energy_scale]: #configuration-energy_scale
[configuration-energy_shift]: #configuration-energy_shift
[configuration-comp]: #configuration-comp
[configuration-prec]: #configuration-prec
[configuration-div]: #configuration-div
[configuration-bound]: #configuration-bound
[configuration-leng]: #configuration-leng
[configuration-type]: #configuration-type
[configuration-custom_pot]: #configuration-custom_pot
[configuration-print_custom_pot]: #configuration-print_custom_pot
[comment]: <> (Class Methods)
[configuration-set_type]: #configuration-set_type

[make_pybinding_model]: #make_pybinding_model
[estimate_bounds]: #estimate_bounds
[config_system]: #config_system
[loud_deprecation_warning]: #louddeprecationwarning

[kitex]: kitex.md
[kitetools]: kite-tools.md
[kitetools-output]: kite-tools.md#output
[HDF5]: https://www.hdfgroup.org
[tutorial-hdf5]: ../documentation/editing_hdf_files.md

[magnetic-field]: ../documentation/magnetic.md