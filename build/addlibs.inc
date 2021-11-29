# here the following variables must be defined
# NCL netCDF lib location
# NCI netCDF include location
# ADDLIB any additional libs required by netCDF lib to compile
ADDLIB = 
NETCDF = no

ifeq ($(FC),gfortran)
	NCL = /home/ipolonsk/libs/lib
	NCI = /home/ipolonsk/libs/include
	INCLUDES:= -I. -I$(NCI) 
	NETCDF = yes
	ADDLIB = -L/home/ipolonsk/libs/hdf5/lib -lhdf5_hl -lhdf5 -ldl -lz
endif

