
## Do I need to use the python script to generate the HDF file?

No. Advanced users can analyse the structure of the HDF file and construct their own hdf file as an input to run KITE. 

## Why do we need a post-processing tool? It could be easier to have a single code to directly calculate the quantities.

One of the main advantages of a Chebyshev expansion is the possibility to obtain new values of a given quantity without having to recalculate the Chebyshev expansion. This is the case, for example, of the DC and optical conductivities. One can use the same Chebyshev expansion to obtain these two quantities at different temperatures or the optical conductivity for a set of Fermi energies. KITE is structured to take advantage of this possibility. KITEx calculates the Chebyshev expansion while KITE-tools uses the expansion to obtain these quantities while allowing the user to change several parameters without the need to recalculate the terms of the expansion.

## What are the units of the magnetic field?  

KITE assumes the lattice constants of the tight-binding model are provided in nanometers (nm) and the magnetic field is given in Tesla. If the user consider other units for the lattice constant, the magnetic field is rescaled according to the definition of magnetic flux per unit cell. 


## What are the units of the temperature? 
For KITE, the temperature parameter is given in units of energy, considering in fact $k_BT$ and follows the definition of the hopping parameters. To convert to Kelvins, the user needs to devide KITE's temperature by $k_BT$.

##How to cite:

```
@ARTICLE{Joao2020,
  title     = "{KITE}: high-performance accurate modelling of electronic
               structure and response functions of large molecules, disordered
               crystals and heterostructures",
  author    = "Jo{\~a}o, Sim{\~a}o M and An{\dj}elkovi{\'c}, Mi{\v s}a and
               Covaci, Lucian and Rappoport, Tatiana G and Lopes, Jo{\~a}o M V
               P and Ferreira, Aires",
  publisher = "The Royal Society",
  volume    =  7,
  number    =  2,
  pages     = "191809",
  month     =  feb,
  year      =  2020,
  keywords  = "Chebyshev expansions; disorder; electronic structure; optical
               response; quantum transport; tight-binding simulations",
  language  = "en"
}
```

