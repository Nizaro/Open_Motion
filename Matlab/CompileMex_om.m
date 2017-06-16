%% compile Mex SAD

OMPATH='/Users/scchia/Desktop/Open_Motion/'
cd ([OMPATH,'Matlab'])

mex -v -O CXXOPTIMFLAGS="\$CXXOPTIMFLAGS -mavx2 -Ofast -march=native" -lopenmotion -lm -D_BSD_SOURCE 'om_Matlab.c' 'utils_.c'