<img src=https://user-images.githubusercontent.com/39924384/41094707-9e4ead6e-6a25-11e8-9e16-070a3236c8da.png width="100">

**KITE** is an open-source Python/C++ software suite for efficient real-space tight-binding (TB) simulations of electronic structure and bulk quantum transport properties of disordered systems scalable to multi billions of atomic orbitals.

KITE v1.1 (this version) is shipped with the following capabilities:

* Quantum Transport: longitudinal and transverse DC conductivities at zero and finite temperature;
* Optical Properties: AC longitudinal conductivity;
* Optical Properties: nonlinear optical susceptibility;
* Unitary time evolution: Gaussian wave-packet propagation (diffusion);
* Electronic Structure: high-resolution density of states (DoS);
* Electronic Structure: local DoS (LDoS) and ARPES;
* Special Features: automated magnetic field in 2D;
* Special Features: lattice perturbations (e.g., strain) and ‘disorder cell’ concept for implementation of multi-orbital/-bond impurities;
* Algorithms: high-resolution CPGF full spectral (DoS);
* Algorithms: high-resolution CPGF full spectral (finite temperature response functions);
* Algorithms: high-resolution CPGF single-shot (zero temperature DC conductivity)

New functionalities in v1.1:

* Algorithms: High-resolution CPGF single-shot algorithm extended to handle multiple Fermi energies simultaneously  
* Special features: twisted boundary conditions for arbitrary lattice models
* Special features: customized on-site potential landscapes
* Special features: automated magnetic field in 3D cubic systems

For further details on the algorithms and the implementation please refer to S. M. Joao et al., R. Soc. Open Sci. 7, 191809 (2020) [https://royalsocietypublishing.org/doi/full/10.1098/rsos.191809].

## Getting Started

For installation instructions and an overview of the main functionalities of **KITE**, please refer the documentation online [Getting Started](https://quantum-kite.com/category/getting-started/).
Advanced examples are explained in [Examples][kite-examples]. After downloading the repository, you can find a copy of the tutorial scripts under the examples folder.

### Prerequisites

Before installing the core components, KITEx and KITE-tools, the following prerequisites need to be satisfied:

* Eigen3

* Python (version 3.5 or newer)

* HDF5 (version 1.8.13 or newer)

* Pybinding ([see](https://github.com/dean0x7d/pybinding) the requirements)

* GCC compiler (version 4.8.1 or newer, for the wavepacket functionality gcc 8.0 is needed)

* CMake

* Make.

### Installation

After meeting prerequisites and downloading the repository, you can compile KITE using the already available Cmake files. Please check the Cmake files first, and make any necessary edits to the libraries/packages' paths. For step by step instructions, please refer to the [Installation](https://quantum-kite.com/installation/) section.

## License

This project is licensed under the GPL v3 License - see the [LICENSE.md](LICENSE.md) file for details

## Project founders

**Simão M. João**, **João V. Lopes** (Universidade do Porto, Portugal), **Tatiana G. Rappoport** (Universidade Federal Rio de Janeiro, Brazil), **Miša Anđelković**, **Lucian Covaci** (University of Antwerp, Belgium) and **Aires Ferreira** (University of York, UK).

## Contributors

* KITEx/KITE-tools development: **João P. Santos Pires** (Porto)
* Python interface / documentation: **Bert Jorissen** (Antwerp), **Emile Aerts** (Antwerp), **Robin Smeyers** (Antwerp), **David T. S. Perkins** (York)

If you would like to collaborate with us on the KITE project, please email a team member directly or send us a message through the [Contact](https://quantum-kite.com/contact/) form.

## Acknowledgments

The inception of KITE's open-source project was funded by **The Royal Society** through grant NA150043 "*Efficient Numerical Solver for Spin, Charge and Optical Conductivity*" (T. Rappoport and A. Ferreira). The KITE team is also grateful for the support received from EPSRC (A. Ferreira), The Royal Society (UF130385, A. Ferreira) and FLAG-ERA (TRANS2DTMD, M. Anđelković and L. Covaci). We thank Miguel Dias Costa (HPC at NUS, Singapore), Killian Murphy (HPC at University of York, UK) and Julia Giannella (Web Design) for the technical support provided in the early stages of the KITE project.

[kite-examples]: https://quantum-kite.com/documentation/examples/