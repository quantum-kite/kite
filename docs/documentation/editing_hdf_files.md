## What is an HDF file?

Some extracts from the HDF group [HDF](https://support.hdfgroup.org/documentation/index.html):

Hierarchical Data Format 5 (HDF5) is a unique open source technology suite for managing data collections of all sizes and complexity. HDF5 has features of other formats, but it can do much more.
HDF5 is similar to XML in that HDF5 files are self-describing and allow users to specify complex data relationships and dependencies.
In contrast to XML documents, HDF5 files can contain binary data (in many representations) and allow direct access to part of the file without first parsing the entire contents.

**HDF5 also allows hierarchical data objects to be expressed in a natural manner (similar to directories and files)**, in contrast to the tables in a relational database.
Whereas relational databases support tables, HDF5 supports n-dimensional datasets and each element in the dataset may itself be a complex object.
Relational databases offer excellent support for queries based on field matching,
but are not well-suited for sequentially processing all records in the database or for selecting a subset of the data based on coordinate-style lookup.

## Editing the file

Thanks to its underlying (Chebyshev) machinery, KITE can easily re-calculate a physical quantity for different choices of parameters at the post-processing level, i.e. without the need for re-calculating Chebyshev moments.
As explained in the [Post-processing tools documentation](postprocessing.md), this can be done by via the several options available in [KITE-tools][kitetools]. Here, we discuss an alternative (more advanced) approach, based upon the editing of the HDF file. 
Suppose we would like to change the post-processing parameters first specified when first creating the HDF file (e.g., the temperature or number of energy points of a conductivity calculation). For that purpose, we provide a simple python script that rewrites specific parts of our .h5 files.
As discussed above, the .h5 contains hierarchical data objects that are similar to the structure of directories and files.

When modifying a parameter like temperature, we need to locate in the .h5 file the quantity that is going to be calculated and modify its value.
The script describes how to list the parameters associated to each quantity and how to edit one parameter.

``` python linenums="1"
file_name = 'archive.h5'
f = h5py.File(file_name, 'r+')     # open the file

# List all groups
print('All groups')
for key in f.keys():  # Names of the groups in HDF5 file.
    print(key)
print()

# Get the HDF5 group
group = f['Calculation']

# Checkout what keys are inside that group.
print('Single group')
for key in group.keys():
    print(key)
print()
#if you want to modify other quantity, check de list and change the subgroup below
# Get the HDF5 subgroup
subgroup = group['conductivity_dc']

# Checkout what keys are inside that subgroup.
print('Subgroup')
for key in subgroup.keys():
    print(key)
print()

new_value = 70
data = subgroup['Temperature']  # load the data
data[...] = new_value  # assign new values to data
f.close()  # close the file

# To confirm the changes were properly made and saved:

f1 = h5py.File(file_name, 'r')
print(np.allclose(f1['Calculation/conductivity_dc/Temperature'].value, new_value))
```


[kitex]: ../api/kitex.md
[kitetools]: ../api/kite-tools.md