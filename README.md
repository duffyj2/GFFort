Should put tips and tricks for running things in here

Compling lapack/blas
gfortran yourfile.f95 /usr/lib64/atlas/liblapack.so.3.0 /usr/lib64/libblas.so.3.4.2 

Realistically this directory is going to be entirely dedicated to your basic methods... NMS and MKL

Moving away from this idea and towards the idea of all SLATEC.
It's portable. That's a thing.
Also, i don't think NMS has a routine for infinite integrals.

Currently, NMS is beating the shit of of QAG. 
But I think this might be because its relative and absolute error are the same thing.
Well, setting both to the same thing in QAG doesn't improve matters.
Maybe you should make use of the flowchart in the QUADPACK website and use a more appropriate routine?
That would be QAGS as far as we currently know, but it report a bunch of errors and is slow too.
It's always possible that NMS is just better. But it might be throwing out a bunch of errors that you don't see. 
You could also optimise your routines for whatever integration you're doing. 

