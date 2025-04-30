## Structure

KITE has three different layers:

* [User interface (*Python*)][kitepython]
* [Main program (*C++*)][kitex]
* [Post-processing tool (*C++*)][kitetools]

 
The [tight-binding][tightbinding] (TB) model is defined on a [Python interface][kitepython] based on [Pybinding]. The TB parameterisation enjoys from a number of KITE-specific advanced features, including disorder patterns and magnetic-field modifications. The KITE model -- which includes the desired target-function calculations, such as DOS and conductivies -- is exported to a [HDF5]-file, together with the calculation settings (i.e. system size, boundary conditions, parallelization options, etc.). This file is then given as an input to the main program (*[KITEx][kitex]*). The *input* and *output* for the main program are written to the same [HDF5] file. The complete workflow is summarized in the figure below.

<div>
  <figure>
    <img src="../../assets/images/getting_started/schematic_kite.png" style="width: 40em" />
    <figcaption>The different components of KITE and its workflow.</figcaption>
  </figure>
</div>

## Steps

1. Start by building a [`#!python pb.Lattice`][lattice] that describes a regular tight-binding model (*Section 2*)
2. Add optional terms to the TB Hamiltonian, including disorder patterns and a magnetic field (covered in *Section 6 and 7*)
3. Specify the _calculations settings_ (*Section 3*) and the desired _target functions_ to be computed (*Section 4*)
4. Export your KITE model to the [HDF5] file and run [KITEx][kitex] (*Section 4*)
5. Post-process the data using the post-processing tools [KITE-tools][kitetools] and visualise the data (*Section 5*)

!!! Tip 
    
    It is possible to use a [simple python script][script] for the whole workflow.



[HDF5]: https://www.hdfgroup.org
[Pybinding]: https://docs.pybinding.site/en/stable
[lattice]: https://docs.pybinding.site/en/stable/_api/pybinding.Lattice.html
[script]: index.md
[tightbinding]: ../background/tight_binding.md

[kitepython]: ../api/kite.md
[kitex]: ../api/kitex.md
[kitetools]: ../api/kite-tools.md
