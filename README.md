# grasp2018_Rmatrix_interface

Main codes:

stg1d0_grasp2k.f90 - parses rwfn.out and rmcdhf.sum to get the radial orbitals and generalised
                     occupation numbers. these are then printed to TARGET.INP, which is used 
                     by the darc codes.

                     compile: gfortran stg1d0_grasp2k.f90 -o stg1d0_grasp2k.x
                     usage  : ./stg1d0_grasp2k.x (by symlink or otherwise) in the same directory 
                     as rwfn.out and rmcdhf.sum.


grasp2k_to_stg2.py - parses rcsfs.out to produce a (condensed) stg2 input deck, DSTG2.INP

                     usage: python3 path_to/grasp2k_to_stg2.py -f (rcsf.out file) -l (num levels in cc) -p (num partial waves, starts at either j=0 or 0.5 depending on parity of electron number.
                     right now it always starts at either of these, need to implement a higher starting point for split runs)


grasp2k_to_nr.py   - parses rcsfs.out and prints out the corresponding nr valence configutations

                     in a format akin to the input deck of grasp0 (http://connorb.freeshell.org)
                     usage: python3 path_to/grasp2k_to_nr.py -f (rcsf.out file)


Libraries:
read_grasp2k.py  - routines for reading grasp stuff
stg2_lib.py      - routines for writing DSTG2.INP