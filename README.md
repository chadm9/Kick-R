# Kick-R
---

## Description
Kick-R is a Python program for carrying out a stochastic search for any type of molecular structure. 
The general approach follows Saunders' seminal <a href='http://onlinelibrary.wiley.com/doi/10.1002/jcc.10407/abstract'>"Kick"</a> procedure, but adds "Restrictions" (-R) on 
what constitues an acceptable guess structure.  These restrictions are that no pair of atoms in a 
guess structure can exhibit distances less than the sum of their combined covalent radii, and that
each atom in a guess structure must be connected to every other atom by a network of bonds (where
bonds are defined based on geometric criteria).  Graph Theory is employed to check for adherence
to the latter restriction.  
<br>
<br>
Compared to Suanders' orginal Kick method, Kick-R delivers substantial improvements
in search efficiency due to its truncation of the molecular configuration space searched to
regions on the potential energy surface which are chemically relevant.  In the present implementation, Kick-R is designed to interface
with the Gaussian package for performing ab initio quantum chemistry optimizations.  However, if use
of another optimization package is desired, this code can be used simply
for guess structure geometry initialization.
<br>
<br>
Please note that this program randomly initiates geometries, and then filters the
resulting guess structures based on the two criteria outlined above.  Due to the random
nature of the geometry intialization, filtering times increase rapidly with system size.
For this reason, use of this code is not recommended for systems larger than 
~10 atoms.



## Author: Chad McKee
* <a href="https://github.com/chadm9">GitHub</a>
* <a href="http://wchadmckee.com/">Personal Website</a>
* <a href="https://www.linkedin.com/in/w-chad-mckee-88939163/">LinkedIn</a>

## Languages and Technologies used
* Python 2.7


## Dependencies and Plugins
None, but note that the present implementation of Kick-R is designed to 
interface with the Gaussian software suite for performing ab initio quantum
mechanical optimizations.

## Usage Example

The following assumes that Kick-R.py, Kick-R.in, and kicksum.py
are in your current working directory.
<br>
<br>
An example input file (Kick-R.in) is provided below:

```
kick_number 250
box_size 8.68x8.68x8.68
multiplicity 2
dft_type PBEPBE
basis_set LANL2DZ
charge 0
structure 7
Au        0.000000000      0.000000000      0.000000000
Au        0.000000000      0.000000000      0.000000000
Au        0.000000000      0.000000000      0.000000000
Au        0.000000000      0.000000000      0.000000000
Au        0.000000000      0.000000000      0.000000000
Au        0.000000000      0.000000000      0.000000000
Au        0.000000000      0.000000000      0.000000000
```

This input will produce 250 candidate structures inside of a
8.68x8.68x8.68 Angstrom^3 box.  The multiplicity, DFT method,
basis set and charge are then set as they would appear in a typical
Gaussian input file.  The structure keyword gives the number of atoms,
and the atoms and their initial coodinates are then listed below.
<br>
<br>
From the command line, Kick-R can be used to generate acceptable candidate structures
by executing the following command:
```
python Kick-R.py
```
This will produce a summary file (Kick-R.out) and 250 Gaussian inpupt files which are to be submitted for 
optimization.  An example of one such input file is:

```
%mem=4gb
%lindaworkers=LINDA
# UPBEPBE/LANL2DZ Opt=maxcycles=50

Title

0 2
Au  2.03907558319       -1.2628997452   1.18246662681
Au  -1.14620277447      2.47635533133   0.610008900135
Au  0.594469179038      1.09074156964   0.00957865796746
Au  -1.27735656662      -0.516589194542 -0.916367620249
Au  -0.435523996789     -1.34989418473  1.6629152073
Au  -2.01855389511      -1.94726269008  -2.8221232037
Au  -0.961636128879     -3.08273676614  -0.12877012753

END
```
Assuming Gaussian has been used to optimize the candidate structures,
the script 'kicksum.py' can then be called to collect all unique results
and order them according to their energy.  A usage example is:
```
python kicksum.py 250 > results
```
This will write the results to a file called 'results'.






## Code Snippet

### Checking for Connectivity in the Molecular Graph
Kick-R checks to ensure that in all candidate structures every atom is connected
to every other atom by a network of bonds.  This is common task in graph theory.
In Python, graphs can be easily represented via a dictionary, where the key values
 pairs represent atoms which are bonded.  Given a molecular graph, the follwoing
 function checks to see if two atoms are connected by a path consisting of 
 bonds between atoms
<br>
```Python
def isConnected(graph, start, end, path=[]):
    path = path + [start]
    if start == end:
        return True
    for node in graph[start]:
        if node not in path:
            newpath = isConnected(graph, node, end, path)
            if newpath: return True
    return False
```
Using the above function, it is now easy to check to see if all
atoms in a candidate structure are connected by a network of bonds.
```
def ConnectedCriteria(graph):
    passed, i, j = True, 1, 2
    while i < len(graph.keys()):
        while j <= len(graph.keys()):
            if not isConnected(graph, i, j):
                passed = False
                break
            j = j + 1
        if passed == False:
            break
        i = i + 1
        j = i + 1
    return passed
```
Please direct any questions regarding useage ofthis code to chadm_@yahoo.com.
