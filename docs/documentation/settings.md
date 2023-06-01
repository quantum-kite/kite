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

## [Divisions][configuration-divisions]
: The [`#!python divisions`][configuration-divisions] is an integer number that defines the number of decomposition parts in each spatial direction.
  [KITEx][kitex] implements a domain decomposition technique to divide the lattice into various partitions that are computed in parallel.
  The domain decomposition is optimized at the design level and allows a substantial speed up of multithreaded calculations, it's usage is recommended.
  
: To activate this feature, set a number of decomposition parts larger than one `#!python nx * ny * nz > 1`.
    !!! Warning
    
        The product `#!python nx * ny * nz` equals the number of threads used by KITEx and thus **must not exceed** the number of **available cores** in the computer.

## [Length][configuration-length]
: The [`#!python length`][configuration-length] is an integer number of unit cells along the direction of lattice vectors `#!python lx, ly, lz = 256, 256, 256`. 
  The lateral size of the decomposed parts are given by `#!python lx/nx` and `#!python ly/ny`.

    !!! Warning
    
        The laterial sizes `#!python lx/nx`, `#!python ly/ny`, `#!python lz/nz` **must be integers**.
          
: When using a 2D lattice, only `#!python lx, ly, nx, ny ` are needed.

## [Boundaries][configuration-boundaries]
:  KITE has 3 standard types of boundary conditions (BCs) implemented, namely: periodic, open, and twisted. Moreover, a "random BCs" option is available, whereby statistical averages over ensembles of random vectors (or disorder configurations) are done with the help of random twist angles drawn from a uniform distribution. This special option is particularly useful to simulate the infinite-size “bulk", since it efficiently eliminates finite size effects.  


    !!! Info
        
        It is possible to impose open BCs along one spatial direction to build ribbons in 2D and slabs in 3D. 
    
    
    The [`#!python boundaries`][configuration-boundaries] is a string, use `#!python 'periodic'` for *periodic* BCs, `#!python 'open'` for *open* BCs, `#!python 'twisted'` for *twisted* BCs and `#!python 'random'` for *random* BCs. In all cases, the system has the geometry of the unit cell, which is replicated `#!python lx, ly, lz ` times in the directions of the unit vectors. Different BCs can be used along the `#!python x, y` and `#!python z` axis. If *twisted* boundary conditions are used, the twist [`#!python angles`][configuration-angles] must be included in radians.
    
    ### Twisted BC

    For twisted BCs, the twist phase angles need to be specified by the user. This is done by means of an extra argument `#!python ' angles=[phi_1,..,phi_DIM]'` where `#!python ' phi_i \in [0, 2*M_PI]'`. The syntax is simple:
    
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

    Random BCs are defined using `mode = 'random'`. No extra arguments are required, but this option implicitely assumes that many random vectors (and/or disorder configurations) will be used. For a single system realization (Sec. [calculation] for calculation settings), this option is equivalent to a mere twisted-BC simulation with randomly chosen twist-angles along `#!python x, y` and `#!python z` axis.


    !!! Warning
    
        The usage of `#!python True` or `#!python False` for the boundaries is *deprecated*.

## [Complex][configuration-is_complex]
: The [`#!python is_complex`][configuration-is_complex] is a boolean value.
  For optimisation purposes, KITEx only considers and stores complex data with the setting `#!python is_complex=True`.
  `#!python False` should be used for real symmetric Hamiltonians.


## [Precision][configuration-precision]
: The [`#!python precision`][configuration-precision] is an integer identifier for the used data type.
  [KITEx][kitex] allows users to define the precision of the calculation.
  Use `#!python 0` for float, `#!python 1` for double, and `#!python 2` for long double.

## [Spectrum Range][configuration-spectrum_range]
: The optional [`#!python spectrum_range`][configuration-spectrum_range] is an array of reals.
  By default, [KITEx][kitex] executes an automated rescaling of the Hamiltonian, see the [Documentation][documentation].
  Advanced users should avoid the automated rescaling and override this feature using `#!python spectrum_range=[Emin,Emax]`, where `#!python Emin, Emax` are the minimum, maximum eigenvalues of the TB matrix. _Lower/upper bounds on smallest/largest energy eigenvalues should be used if exact eigenvalues are unknown_ _(often the case in systems with disorder); see Sec. [Disorder] for more information_. 


    To manually set the [`#!python spectrum_range`][configuration-spectrum_range], it is necessary to add an extra parameter
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
[documentation]: ../documentation/index.md
[tightbinding]: ../documentation/tight_binding.md

[lattice-tutorial]: tb_model.md

[kitepython]: ../api/kite.md
[kitex]: ../api/kitex.md
[kitetools]: ../api/kite-tools.md

[calculation]: index.md
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
[calculation]: ../api/kite.md#calculation
