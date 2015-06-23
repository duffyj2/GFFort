Should put tips and tricks for running things in here

Compling lapack/blas
gfortran yourfile.f95 /usr/lib64/atlas/liblapack.so.3.0 /usr/lib64/libblas.so.3.4.2 

Realistically this directory is going to be entirely dedicated to your basic methods... NMS and MKL

Moving away from this idea and towards the idea of all SLATEC.
It's portable. That's a thing.
Also, i don't think NMS has a routine for infinite integrals.

Currently, NMS is beating the shit of of QAG and QAGS. 
Even if abs/rel error are taken into account. 
It's always possible that NMS is just better. But it might be throwing out a bunch of errors that you don't see. 
You could also optimise your routines for whatever integration you're doing. 
And you could just include numerical routines to do whatever type of integration you want.

