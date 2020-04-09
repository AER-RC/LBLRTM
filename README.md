# LBLRTM

LBLRTM (Line-By-Line Radiative Transfer Model) is an accurate and efficient line-by-line radiative transfer model derived from the Fast Atmospheric Signature Code (FASCODE). LBLRTM has been, and continues to be, extensively validated against atmospheric radiance spectra from the ultraviolet to the sub-millimeter.

The [HITRAN database](http://cfa-www.harvard.edu/hitran) provides the basis for the line parameters used in LBLRTM. These line parameters, as well as additional line parameters from other sources, are extracted for use in LBLRTM by a line file creation program called LNFL. A line parameter database built from HITRAN and suitable for use with LNFL is available from the [AER RT web site](http://rtweb.aer.com).

[Add plantUML diagram]

# Cloning

Assuming the output directory should be `LBLRTM`:

```
% git clone --recursive git@github.com:AER-RC/LBLRTM.git
```

`--recursive` is important, because this repository is linked with our [common FORTRAN modules repository](https://github.com/AER-RC/aer_rt_utils) that are required in the model builds. If this keyword is forgotten, one can do:

```
git submodule init
git submodule update
```

in the `LBLRTM` directory.

# General LNFL/LBLRTM File Information

## Platforms on which LBLRTM can be run

It is recommended that LNFL and LBLRTM be compiled in Fortran 90. LBLRTM has previously been run on DEC alpha, Cray, MS-DOS, and HP platforms.

Some users have ported the code to the Windows/DOS environment. AER presently does not officially support this implementation; however, the following description of how LBLRTM was used in XP by a user (Christopher Rice, Air Force Institute of Technology [AFIT]):

* Obtain the newest Intel Fortran Compiler (v 9.1) and Visual Studio.net 2003.
* Install visual studio .net 2003.
*	Install Intel Fortran Compiler (Fortran compiler options will appear in MSVS 2003 after
the Intel Fortran Compiler is installed.
*	Decompress source codes for LBLRTM and LNFL into their appropriate folders.
*	Compile LNFL:
  * Open Visual Studio 2003.
	* Open a new project.
  * Select “Intel Fortran Projects” under “Project Types” and choose “Console Application” from “Templates”.
  * Name this project accordingly. (e.g., `LNFL_Exe`)
  * After the project is created a “Solution Explorer” will show the solution named as above. Right click on “Source Files” and add the following files:
      * `lnfl.f90`
      * `util_dos.f90`
  * Note that the `util_linux_intel` makefile can be used as a reference for the files to be included, replacing `util_linux.f90` with the file `util_dos.f90`.
  * Make sure the compiler is set to “release” NOT “debug”.
  * Build the Project.
  * On successful build find the `.exe` in the “release” folder where the project is saved
* Compile LBLRTM:
  * Open Visual Studio 2003.
  * Open a new project.
  * Select “Intel Fortran Projects” under “Project Types” and choose “Console Application” from “Templates”.
  * Name this project accordingly. (e.g., `LBLRTM_Exe`)
  * After the project is created a “Solution Explorer” will show a solution named as above. Right click on “Source Files” and add the following files:
    * `contnm.f90`
    *	`fftscn.f90`
    *	`lblatm.f90`
    *	`lbldum.f90`
    *	`lbllow.f90`
    *	`lblrtm.f90`
    *	`nonlte.f90`
    * `oprop.f90`
    *	`pltlbl.f90`
    *	`postsub.f90`
    *	`solar.f90`
    *	`testmm.f90`
    *	`util_dos.f90`
    *	`xmerge.f90`
  * Note that the `util_linux_intel` makefile can be used as a reference for the files to be included, replacing `util_linux.f90` with the file `util_dos.f90`.
  * Make sure the compiler is set to “release” NOT “debug” – LBLRTM will not operate correctly when compiled in “Debug” mode.
  * Build the Project.
  * On successful build find the exe in the “release” folder where the project is saved.
* Use LNFL with `TAPE5` to create `TAPE3` as described in the documentation,
* Use `TAPE3` with LBLRTM to satisfy your requirements as described in documentation.

## Issues relating to unformatted files on UNIX and LINUX systems

Unformatted files are often not compatible between systems due to differences in the way the bytes are written to the files (big-endian versus little-endian).  Note that the `byteswap` option available with most compilers will not work with most LBLRTM unformatted output files because of the mixing of real and integer data within records.

## LNFL/LBLRTM Naming Convention

Specific information on the input/output files from LNFL and LBLRTM is located in their respective input files, `lnfl_instructions` and `lblrtm_instructions`, and the examples provided in the code `tar` files.

Most file names are given as `TAPEx` where `x` is a one- or two-digit number.  The name is case-sensitive, and is uppercase.  Tape numbers may be same for LNFL and LBLRTM but do not represent identical files. For example, the primary LNFL input file is `TAPE5`, and the primary LBLRTM input file is `TAPE5`.  However, they have neither the same input information nor the same formatting. The instruction manual for each code details the input file information.

## LNFL/LBLRTM Input File (TAPE5) Format

The `TAPE5` input files are read as formatted FORTRAN. As a consequence of the formatted read, any blank space will be read as "zero".  Thus, one may leave blanks for most of the parameters and within the code they will default to an acceptable value.

Real numbers format input as either `E` or `F` format, with the entire number within the range specified in the input instructions.  Integers must be specified exactly in the integer format.  For example, the spectral bandwidth (_v<sub>1</sub>_ to _v<sub>2</sub>_) in LBLRTM `TAPE5` is input as 10 character real numbers. This means that the value can be written anywhere within these 10 characters, as long as there is a decimal point (e.g. "---600.000" or "-600.000--", where "-" is a blank space).

Integers are read in with the `I` format.  For example, the model atmosphere (`iatm`) in LBLRTM `TAPE5` is input as `I5`, so it must be "----2", and not "2----" as this will be read as 20000.

## LBLRTM Output File Format

The general structure of the files involves the use of panels, which are blocks of output usually containing 2400 points. Each panel contains a header to describe the starting and ending points of the panel (_v<sub>1</sub>_ and _v<sub>2</sub>_), the spectral spacing of the points (`dvp`), and the number of points in the panel (`npts`). The panel header is followed by either one or two (see below) blocks of output, consisting of `npts` points.

`TAPE12`: Radiances and transmittances
(1) file header
(2,i)-panel header
(3,i) radiances
(4,i) transmittances
Lines 2-4 repeat for i=1,N times to cover the entire spectral region.

`TAPE11`: Filtered radiance or transmittance (also applies to any user-designated output file which contains radiances, transmittances, or optical depths, such as the "ODint" file)
(1)file header
(2,i)-panel header
(3,i) radiances or transmittances
Lines 2-3 repeat for i=1,N times to cover the entire spectral region.

Note that a limited amount of spectral output information may also be put in the `TAPE6` using the `MPTS`/`NPTS` options of `TAPE5` record 1.2.

# Instructions and Tips for Running LNFL

LNFL is used to generate a unformatted file (`TAPE3`) of all the line parameters required by LBLRTM.

## Input files for LNFL

1. `TAPE1`: The line parameter database in ASCII format (also available on [RTWeb](www.rtweb.aer.com)).

2. `TAPE5`: LNFL input file.

## Output files for LNFL

1. `TAPE3`: Unformatted LNFL output file containing the line parameters for LBLRTM.
2. `TAPE6`: Informational output file.
3. `TAPE7`: Optional output file containing ASCII version of the parameters contained in `TAPE3`.

## Sequence for running LNFL

*	Download latest LNFL tar code (containing the source code) and the latest line parameter database from [RTWeb](rtweb.aer.com).
*	Compile LNFL using the makefiles found in the LNFL tar file.  Note: one needs to compile in the `build` directory.
*	Link the line parameter database to `TAPE1` in the LNFL working directory.
*	Remove `TAPE3` file from the LNFL working directory.
*	Edit necessary parameters in the `TAPE5` input file.  Note that the beginning and ending wavenumber (_v<sub>1</sub>_, _v<sub>2</sub>_) in `TAPE5` must extend at least 25 cm<sup>-1</sup> beyond each end of the desired spectral range for the LBLRTM calculations.
*	Run the LNFL executable.

# Instructions and Tips for Compiling and Running LBLRTM

LBLRTM is used to generate line-by-line upwelling and downwelling transmittances and radiances.

## Required input files for LBLRTM
1. `TAPE3`: Unformatted file containing line parameter information, generated by LNFL (see above). The `TAPE3` file should include lines from at least 25 cm<sup>-1</sup> on either end of the calculation region.
2. `TAPE5`: Input file required to run LBLRTM.

The spectral interval (_v<sub>1</sub>_, _v<sub>2</sub>_) for any LBLRTM run must not exceed 2000 cm<sup>-1</sup> (see instruction manual).

Other input files are required if you are using the solar source function, cross sections, surface emissivity, etc. See the LBLRTM instruction manual and provide example.

## Layer numbering scheme

The LBLRTM convention is that layer 1 is at the highest pressure level (lowest altitude).  The layer information for a given run may be found in `TAPE6`.

## Output files for LBLRTM
1. `TAPE6`: Informational output file
2. `TAPE11`: Unformatted file containing filtered output, if requested in `TAPE5`.
3. `TAPE12`: Unformatted file containing transmittances/radiances.

ASCII file of unformatted unformatted files can be requested in the LBLRTM TAPE5 (see `pltlbl` variable in Record 12).

Unformatted optical depth files can be requested in the LBLRTM using options specified in `TAPE5`.

## Sequence for running LBLRTM
* Download latest LBLRTM tar code (containing the source code) and the latest line parameter database from [RTWeb](http://rtweb.aer.com).
* Compile LBLRTM following makefiles in the LBLRTM tar file. Note, one needs to compile in the `build` directory
* Link the line parameter database (`TAPE3` from LNFL) to the LBLRTM working directory.
* Edit any parameters necessary in the input file `TAPE5`.
* Run the LBLRTM executable.

# General Questions

1. What is the difference between a line-by-line calculation and a band-model calculation?

Absorption/emission spectra are comprised of a complicated array of spectral lines.  The HITRAN 2008 Database (Version 13.0) contains over 2,713,000 lines for 39 different molecules. In order to resolve these individual lines, a nominal spectral sampling rate of less than the mean line half width must be utilized.  Such highly resolved radiative transfer calculations are called line-by-line (LBL) calculations. The computational time associated with calculating broadband fluxes from LBL calculations is formidable.  A band model aims to simplify radiative transfer calculations by using approximations to represent the line-by-line characteristics of a particular spectral interval.  Band models are appropriate for situations where the desired spectral resolution is much smaller than the Lorentz and Doppler widths of the spectral lines. Such approximations are also of use in general circulation models.
