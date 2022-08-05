KITE uses the classes [`#!python kite.Configuration`][configuration] and [`#!python kite.Calculation`][calculation] to define the calculation settings.

The class [`#!python kite.Configuration`][configuration] carries the following information:

[Divisions][configuration-divisions]
: The [`#!python divisions`][configuration-divisions] is an integer number that defines the number of decomposition parts in each spatial direction.
  [KITEx][kitex] implements a domain decomposition technique to divide the lattice into various partitions that are computed in parallel.
  The domain decomposition is optimized at the design level and allows a substantial speed up of multithreaded calculations, it's usage is recommended.
  
: To activate this feature, set a number of decomposition parts larger than one `#!python nx * ny *nz > 1`.
    !!! Warning
    
        The product `#!python nx * ny *nz` equals the number of threads used by KITEx and thus **must not exceed** the number of **available cores** in the computer.

[Length][configuration-length]
: The [`#!python length`][configuration-length] is an integer number of unit cells along the direction of lattice vectors `#!python lx, ly, lz = 256, 256, 256`. 
  The lateral size of the decomposed parts are given by `#!python lx/nx` and `#!python ly/ny`.

    !!! Warning
    
        The laterial sizes `#!python lx/nx`, `#!python ly/ny`, `#!python lz/nz` **must be integers**.
        
When using a 2D lattice, only `#!python lx, ly, nx, ny ` are needed

[Boundaries][configuration-boundaries]
: The [`#!python boundaries`][configuration-boundaries] is a string, use `#!python 'periodic'` for periodic boundary conditions and `#!python 'open'` for open boundary conditions.
  Additionally, *twisted* and *random twisted* boundary conditions can be implemented using `#!python 'twisted'` and `#!python 'random'` respectively.
  If *twisted* boundary conditions are used, the twist [`#!python angles`][configuration-angles] must be included in radians.

    !!! Info
    
        The usage of `#!python True` or `#!python False` for the boundaries is *deprecated*.

[Complex][configuration-is_complex]
: The [`#!python is_complex`][configuration-is_complex] is a boolean value.
  For optimisation purposes, KITEx only considers and stores complex data with the setting `#!python is_complex=True`.
  `#!python False` should be used for real symmetric Hamiltonians.


[Precision][configuration-precision]
: The [`#!python precision`][configuration-precision] is an integer identifier for the used data type.
  [KITEx][kitex] allows users to define the precision of the calculation.
  Use `#!python 0` for float, `#!python 1` for double, and `#!python 2` for long double.

[Spectrum Range][configuration-spectrum_range]
: The optional [`#!python spectrum_range`][configuration-spectrum_range] is an array of reals.
  By default, [KITEx][kitex] executes an automated rescaling of the Hamiltonian, see the [Documentation][documentation].
  Advanced users can override this feature using `#!python spectrum_range=[Emin,Emax]`, where `#!python Emin, Emax` are the minimum, maximum eigenvalues of the TB matrix.

The [`#!python kite.Configuration`][configuration] object is thus structured in the following way:

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
[disorder]: index.md
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
