gfortran -std=f2003 -I. -Wall -O2 -o run parameters.f03 connectivity.f03 initialization.f03 element.f03 explicit.f03 BC.f03 extract.f03 main.f03 -L/share/apps/lib64 -lblas
