KITE is an open-source Python/C++ software suite for real-space tight-binding (TB) simulations of electronic structure and bulk quantum transport properties of disordered systems scalable to multi billions of atomic orbitals.

KITE v1.0 is shipped with the following capabilities (for clean and disordered systems in two and three spatial dimensions):

* Quantum Transport: longitudinal and transverse DC conductivities at zero and finite temperature;
* Optical Properties: AC longitudinal conductivity;
* Optical Properties: nonlinear optical susceptibility;
* Unitary time evolution: Gaussian wave-packet propagation (diffusion);
* Electronic Structure: high-resolution density of states (DoS);
* Electronic Structure: local DoS and ARPES;
* Special Features: automated magnetic field in 2D;
* Special Features: lattice perturbations (e.g., strain) and ‘disorder cell’ concept for implementation of multi-orbital/-bond impurities;
* Algorithms: ultra-high-resolution CPGF full spectral (DoS);
* Algorithms: high-resolution CPGF full spectral (finite temperature response functions);
* Algorithms: ultra-high-resolution CPGF single-shot (zero temperature DC conductivity)

For further details on the algorithms and the implementation please refer to S. M. Joao et al., [R. Soc. Open Sci. 7, 191809 (2020)](https://royalsocietypublishing.org/doi/full/10.1098/rsos.191809).

## Getting Started

For understanding the main functionalities, how to setup a TB model and disorder and evaluate the KITE target functions, we suggest you to check the tutorial [Getting Started](tutorial/index.md).
More advanced calculations are explained in [Examples](tutorial/examples). After downloading the repository, you can find a copy of the tutorial scripts under the examples folder.

### Prerequisites

Before installing the code, following prerequisites need to be satisfied:

* Eigen3
* Python (version 3.5 or newer)
* HDF5 (version 1.8.13 or newer)
* Pybinding ([see](https://github.com/dean0x7d/pybinding) the requirements)
* GCC compiler (version 4.8.1 or newer, for the wavepacket functionality gcc 8.0 is needed)
* CMake
* Make

### Installing

After meeting prerequisites and downloading the repository, you can compile the KITE code using the already available CMake files. Please check the CMake files first, and if required libraries/packages are installed elsewhere, edit them accordingly. For the full installation procedure, please refer to the [Installation](installation.md) section.

## Contact 

Please share your feedback and bug reports.
We'd love to hear your ideas for future improvements (contact email: support@quantum-kite.com).

## Contributors

If you would like to collaborate with us on the KITE project, send us a message to our email.

## Authors

Simão M. João, João V. Lopes (Universidade do Porto, Portugal), Tatiana G. Rappoport (Universidade Federal Rio de Janeiro, Brazil), Miša Anđelković, Lucian Covaci (University of Antwerp, Belgium) and Aires Ferreira (University of York, UK).

## License

This project is licensed under the GPL v3 License – see the LICENSE.md file for details

## Acknowledgments

KITE’s open source project was funded by The Royal Society through grant NA150043 *“Efficient Numerical Solver for Spin, Charge and Optical Conductivity”* (T. Rappoport and A. Ferreira). KITE’s team thanks the partial support from EPSRC (A. Ferreira), The Royal Society (UF130385, A. Ferreira) and FLAG-ERA (TRANS2DTMD, M. Anđelković and L. Covaci). We thank Miguel Dias Costa (HPC at NUS, Singapore), Killian Murphy (HPC at University of York, UK) and Julia Giannella (Web Design) for technical support.

[contact]: #contact
