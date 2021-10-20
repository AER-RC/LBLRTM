# here the following variables must be defined
# NCL netCDF lib location
# NCI netCDF include location
# ADDLIB any additional libs required by netCDF lib to compile
ADDLIB = 
NETCDF = no

ifeq ($(FC),gfortran)
	NCL = /usr/lib
	NCI = /usr/include
	INCLUDES:= -I. -I$(NCI) 
	NETCDF = yes
endif

