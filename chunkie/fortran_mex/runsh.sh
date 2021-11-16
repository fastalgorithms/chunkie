gfortran -c -O3 -march=native -std=legacy hank103.f -o hank103.o
gfortran -c -O3 -march=native -std=legacy hank103_wrap.f -o hank103_wrap.o
/Applications/MATLAB_R2021a.app/bin/mex -v hank103.c hank103.o hank103_wrap.o -largeArrayDims -DMWF77_UNDERSCORE1 -D_OPENMP -L/usr/local/lib -lfmm3d -lgomp -lgfortran -output hank103_jgh -L/usr/local/lib/gcc/10
cd ../+chnk/+helm2d
cp -rf besselh01_mex.m besselh01.m


