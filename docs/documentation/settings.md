KITE uses the classes [`#!python kite.Configuration`][configuration] and [`#!python kite.Calculation`][calculation] to define the calculation settings.

This is what a typical [`#!python kite.Configuration`][configuration] object (for a 2D lattice) looks like:

``` python linenums="1"
nx = ny = 2
lx = ly = 128
mode = 'periodic'
configuration = kite.Configuration(
  divisions=[nx, ny],
  length=[lx, ly],
  boundaries=[mode, mode],
  is_complex=False,
  precision=1 
)
```
Below, we explain each of the arguments in the [`#!python kite.Configuration`][configuration] object.

!!! Info

     The current syntax assumes that the user will simulate a 2D or 3D system as a default. However, 1D lattices can be easily constructed by simulating a strip geometry with open boundary conditions in the transverse (y) direction and `#!python ly, ny = 1`.

## [Divisions][configuration-divisions]
: The [`#!python divisions`][configuration-divisions] is an integer number that defines the number of decomposition parts in each spatial direction.
  [KITEx][kitex] implements an efficient domain decomposition technique to divide the lattice into various partitions that are computed in parallel.
[KITEx][kitex]  decomposition algorithms are optimized at the design level to deliver optimal multithreading scaling. Their usage for large-scale simulations is highly recommended.
  
: To activate this feature, set a number of decomposition parts larger than one, i.e. `#!python nx * ny > 1` (2D) or `#!python nx * ny * nz > 1` (3D).
    !!! Warning
    
        The product `#!python nx * ny (* nz)` equals the number of threads used by KITEx and thus **must not exceed** the number of **available cores** in the compute node.

## [Length][configuration-length]
: The [`#!python length`][configuration-length] is an integer number of unit cells along the direction of lattice vectors (for example, `#!python lx, ly, lz = 256, 256, 256`). 
  The lateral size of the decomposed parts are given by `#!python lx/nx`, `#!python ly/ny` and  `#!python lz/nz`.

    !!! Warning
    
        The lateral sizes `#!python lx/nx`, `#!python ly/ny`, `#!python lz/nz` **must be integers**.
          
: Note that when using a 2D lattice, only `#!python lx, ly, nx, ny ` are needed.

## [Boundaries][configuration-boundaries]
:  KITE's default boundary conditions (BCs) are: periodic, open, and twisted. Moreover, a "random BCs" option is available, whereby statistical averages over ensembles of random vectors (or disorder configurations) are accelerated via the use of random twist angles drawn from a uniform distribution. This special option is particularly useful to simulate the infinite-size “bulk", since it efficiently eliminates finite size effects.  


   
    
    
    The [`#!python boundaries`][configuration-boundaries] is a string. 
    Use `#!python 'periodic'` for *periodic* BCs, `#!python 'open'` for *open* BCs, `#!python 'twisted'` for *twisted* BCs and `#!python 'random'` for *random* BCs. 
    In all cases, the system has the geometry of the unit cell, which is replicated `#!python lx, ly, lz ` times in the directions of the unit vectors. 
    If *twisted* boundary conditions are used, the twist [`#!python angles`][configuration-angles] must be provided in radians.

    !!! Info

        Different BCs can be used along the different directions. 
        For example, impose open BCs along one spatial direction to build <b>ribbons</b> in 2D and <b>slabs</b> in 3D. 
 
    ### Twisted BC

    For twisted BCs, the twist phase angles need to be specified by the user. This is done by means of an extra argument `#!python ' angles=[phi_x,phi_y,phi_z]'` where `#!python 'phi_{x,y,z} \in [0, 2*M_PI]'`. The syntax is simple. For example, for a twist angle of pi/2.0 along both lattice directions in a 2D system, we can use: 
    
    ``` python linenums="1"
    nx = ny = 2
    lx = ly = 128
    mode = 'twisted'
    twsx = twsy = np.pi/2.0 
    
    configuration = kite.Configuration(
      divisions=[nx, ny],
      length=[lx, ly],
      boundaries=[mode, mode],
      angles = [twsx,twsy]
      is_complex=False,
      precision=1 
    )
    ```

    ### Random BC

    Random BCs are defined using `mode = 'random'`. No extra arguments are required, but this option implicitly assumes that many random vectors (and/or disorder configurations) will be used  (see Sec. [Calculation][calculation]). For a single system realization, this option is equivalent to a mere twisted-BC simulation with randomly chosen twist-angles along `#!python x, y` and `#!python z` lattice directions.


## [Complex][configuration-is_complex]
: The [`#!python is_complex`][configuration-is_complex] is a boolean value.
  For optimization purposes, [KITEx][kitex] only considers and stores complex data with the setting `#!python is_complex=True` activated.
  `#!python is_complex=False` should be used for real symmetric Hamiltonians.


## [Precision][configuration-precision]
: The [`#!python precision`][configuration-precision] is an integer identifier for the used data type.
  [KITEx][kitex] allows users to define the precision of the calculation.
  Use `#!python 0` for float, `#!python 1` for double, and `#!python 2` for long double.

## [Spectrum Range][configuration-spectrum_range]
: The optional [`#!python spectrum_range`][configuration-spectrum_range] is an array of real values.
  By default, KITE executes an automated rescaling of the Hamiltonian (Sec. [Documentation][documentation]).
  Advanced users are encouraged to override this feature and specify the energy interval manually using `#!python spectrum_range=[Emin,Emax]`, where `#!python Emin, Emax` are the minimum, maximum energy eigenvalues.
  

    _Lower/upper bounds should be used if exact energy eigenvalues of the lattice model are unknown (e.g. due to the presence of a disorder landscape); see Sec. [Disorder] for more information_. 

    To manually set the energy spectrum limits, it is necessary to add an extra parameter ([`#!python spectrum_range`][configuration-spectrum_range])
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

[HDF5]: https://www.hdfgroup.org
[pybinding]: https://docs.pybinding.site/en/stable
[lattice]: https://docs.pybinding.site/en/stable/_api/pybinding.Lattice.html
[documentation]: ../documentation/optimization.md
[tightbinding]: ../documentation/tight_binding.md

[lattice-tutorial]: tb_model.md

[kitepython]: ../api/kite.md
[kitex]: ../api/kitex.md
[kitetools]: ../api/kite-tools.md

[calculation]: calculation.md
[DOS]: index.md
[conductivity]: index.md
[modifications]: index.md
[Disorder]: disorder.md 
[Examples]: examples/graphene.md

[configuration]: ../api/kite.md#configuration
[configuration-divisions]: ../api/kite.md#configuration-divisions
[configuration-length]: ../api/kite.md#configuration-length
[configuration-boundaries]: ../api/kite.md#configuration-boundaries
[configuration-is_complex]: ../api/kite.md#configuration-is_complex
[configuration-precision]: ../api/kite.md#configuration-precision
[configuration-spectrum_range]: ../api/kite.md#configuration-spectrum_range
[configuration-angles]: ../api/kite.md#configuration-angles
[configuration-custom_local]: ../api/kite.md#configuration-custom_local
[configuration-custom_local_print]: ../api/kite.md#configuration-custom_local_print
[calculation]: calculation.md
