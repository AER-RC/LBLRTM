# LBLRTM

---
**Contents**

1. [Introduction](#intro)
2. [Cloning the Latest Release](#cloning)
3. [LBLRTM and Docker](#docker)
4. [General LNFL/LBLRTM File Information](#general)
    1. [Platforms on which LBLRTM can be run](#platforms)
    2. [Issues relating to unformatted files on UNIX and LINUX systems](#unformatted)
	3. [LNFL/LBLRTM Naming Convention](#nomenclature)
	4. [LNFL/LBLRTM Input File (TAPE5) Format](#tape5)
	5. [LBLRTM Output File Format](#lblout)
5. [Instructions and Tips for Running LNFL](#runlnfl)
	1. [Input files for LNFL](#lnflin)
	2. [Output files for LNFL](#lnflout)
	3. [Sequence for running LNFL](#lnflseq)
6. [Instructions and Tips for Compiling and Running LBLRTM](#runlbl)
	1. [Required input files for LBLRTM](#lblin)
	2. [Layer numbering scheme](#laynum)
	3. [Output files for LBLRTM](#lblout)
	4. [Sequence for running LBLRTM](#lblseq)
7. [Tests](#tests)
8. [Frequently Asked Questions](#faq)

# Introduction <a name="intro"></a>

LBLRTM (Line-By-Line Radiative Transfer Model) is an accurate and efficient line-by-line radiative transfer model derived from the Fast Atmospheric Signature Code (FASCODE). LBLRTM has been, and continues to be, extensively validated against atmospheric radiance spectra from the ultraviolet to the sub-millimeter.

The [HITRAN database](http://cfa-www.harvard.edu/hitran) provides the basis for the line parameters used in LBLRTM. These line parameters, as well as additional line parameters from other sources, are extracted for use in LBLRTM by a line file creation program called LNFL. A line parameter database built from HITRAN and suitable for use with LNFL can be downloaded with the [AER Line File retrieval code](https://github.com/AER-RC/AER_Line_File) or directory from the [Zenodo repository](https://zenodo.org/record/4019086).

LBLRTM uses the line parameters and [MT_CKD continuum](https://github.com/AER-RC/MT_CKD) in its calculations. The models and data are thus linked. For the latest release, the relationships are:

| LBLRTM Release | MT_CKD Release | Line File |
| :---: | :---: | :---: |
| [v12.11](https://github.com/AER-RC/LBLRTM/releases/tag/v12.11) | [v3.5](https://github.com/AER-RC/MT_CKD/releases/tag/v3.5) | [v3.8](https://zenodo.org/record/4019086/files/aer_v_3.8.tar.gz?download=1) |

If any build or run issues occur, please [create an issue](https://github.com/AER-RC/LBLRTM/issues) or contact the [AER-RC Group](https://github.com/AER-RC).

[Add plantUML diagram]

For more, please see the [Wiki page](https://github.com/AER-RC/LBLRTM/wiki/)

# Cloning the Latest Release <a name="cloning"></a>

Assuming the output directory should be `LBLRTM`:

```
% git clone --recursive git@github.com:AER-RC/LBLRTM.git
```

`--recursive` is important, because this repository is linked with our [common FORTRAN modules repository](https://github.com/AER-RC/aer_rt_utils) that are required in the model builds. The [cross section database](https://github.com/AER-RC/cross-sections) is also added as a submodule (it is not required for all model runs). If this keyword is forgotten, one can do:

```
git submodule init
git submodule update
```

in the `LBLRTM` directory.

Currently, the latest release is LBLRTM v12.11, and it is recommended that this be the version that users clone and checkout (rather than the `master` branch). To do this, one needs to simply checkout the `v12.11` tag:

```
git checkout tags/v12.11
```

No releases before v12.9 are available via GitHub, but they can be requested by contacting the [AER-RC Group](https://github.com/AER-RC). For information on previous releases, please visit the [What's New Wiki page](https://github.com/AER-RC/LBLRTM/wiki/What's-New).

Instead of cloning, users can also download an LBLRTM [tarball](https://github.com/AER-RC/LBLRTM/archive/v12.11.zip) and unpack it:

```
tar xvf v12.11.tar.gz
mv LBLRTM-12.11/ lblrtm
```

Though not necessary, the move to `lblrtm` is for consistency with previous release packages and the associated documentation.

# LBLRTM and Docker <a name="docker"></a>

More doc to come, but see the [GitHub package page](https://github.com/AER-RC/LBLRTM/packages/200551) for Docker image pull directions. And to run:

```
docker pull docker.pkg.github.com/aer-rc/lblrtm/lblrtm:latest

docker tag docker.pkg.github.com/aer-rc/lblrtm/lblrtm:latest lblrtm

docker run -it --rm -v ~/Work/RC/LBLRTM/LBL_In:/LBLRTM/LBLRTM_In -v ~/Work/RC/LBLRTM/LBL_Out:/LBLRTM/LBLRTM_Out lblrtm
```

Volume mounts are necessary to provide LBLRTM inputs and for the user to have access to the outputs. Currently, the `TAPE3`, `TAPE5`, and cross section database are assumed to be in the LBLRTM input directory. `EMISSIVITY` and `REFLECTIVITY` could conceivably work with the correct volume mounts. Cross sections will be their own submodule at some point. The [LBLRTM input file naming convention](#lblin) is assumed.

# General LNFL/LBLRTM File Information <a name="general"></a>

## Platforms on which LBLRTM can be run <a name="platforms"></a>

It is recommended that LNFL and LBLRTM be compiled in Fortran 90. LBLRTM has previously been run on DEC alpha, Cray, MS-DOS, and HP platforms.

Some users have ported the code to the Windows/DOS environment. AER presently does not officially support this implementation; however, the following description of how LBLRTM was used in XP by a user (Christopher Rice, Air Force Institute of Technology [AFIT]):

> * Obtain the newest Intel Fortran Compiler (v 9.1) and Visual Studio.net 2003.
> * Install visual studio .net 2003.
> *	Install Intel Fortran Compiler (Fortran compiler options will appear in MSVS 2003 after
> the Intel Fortran Compiler is installed.
> *	Decompress source codes for LBLRTM and LNFL into their appropriate folders.
> *	Compile LNFL:
>   * Open Visual Studio 2003.
> 	* Open a new project.
>   * Select “Intel Fortran Projects” under “Project Types” and choose “Console Application” from “Templates”.
>   * Name this project accordingly. (e.g., `LNFL_Exe`)
>   * After the project is created a “Solution Explorer” will show the solution named as above. Right click on “Source Files” and add the following files:
>       * `lnfl.f90`
>       * `util_dos.f90`
>   * Note that the `util_linux_intel` makefile can be used as a reference for the files to be included, replacing `util_linux.f90` with the file `util_dos.f90`.
>   * Make sure the compiler is set to “release” NOT “debug”.
>   * Build the Project.
>   * On successful build find the `.exe` in the “release” folder where the project is saved
> * Compile LBLRTM:
>   * Open Visual Studio 2003.
>   * Open a new project.
>   * Select “Intel Fortran Projects” under “Project Types” and choose “Console Application” from “Templates”.
>   * Name this project accordingly. (e.g., `LBLRTM_Exe`)
>   * After the project is created a “Solution Explorer” will show a solution named as above. Right click on “Source Files” and add the following files:
>     * `contnm.f90`
>     *	`fftscn.f90`
>     *	`lblatm.f90`
>     *	`lbldum.f90`
>     *	`lbllow.f90`
>     *	`lblrtm.f90`
>     *	`nonlte.f90`
>     * `oprop.f90`
>     *	`pltlbl.f90`
>     *	`postsub.f90`
>     *	`solar.f90`
>     *	`testmm.f90`
>     *	`util_dos.f90`
>     *	`xmerge.f90`
>   * Note that the `util_linux_intel` makefile can be used as a reference for the files to be included, replacing `util_linux.f90` with the file `util_dos.f90`.
>   * Make sure the compiler is set to “release” NOT “debug” – LBLRTM will not operate correctly when compiled in “Debug” mode.
>   * Build the Project.
>   * On successful build find the exe in the “release” folder where the project is saved.
> * Use LNFL with `TAPE5` to create `TAPE3` as described in the documentation,
> * Use `TAPE3` with LBLRTM to satisfy your requirements as described in documentation.

## Issues relating to unformatted files on UNIX and LINUX systems <a name="unformatted"></a>

Unformatted files are often not compatible between systems due to differences in the way the bytes are written to the files (big-endian versus little-endian).  Note that the `byteswap` option available with most compilers will not work with most LBLRTM unformatted output files because of the mixing of real and integer data within records.

## LNFL/LBLRTM Naming Convention <a name="nomenclature"></a>

Specific information on the input/output files from LNFL and LBLRTM is located in their respective input files, `lnfl_instructions` and `lblrtm_instructions`, and the examples provided in the code `tar` files.

Most file names are given as `TAPEx` where `x` is a one- or two-digit number.  The name is case-sensitive, and is uppercase.  Tape numbers may be same for LNFL and LBLRTM but do not represent identical files. For example, the primary LNFL input file is `TAPE5`, and the primary LBLRTM input file is `TAPE5`.  However, they have neither the same input information nor the same formatting. The instruction manual for each code details the input file information.

## LNFL/LBLRTM Input File (TAPE5) Format <a name="tape5"></a>

The `TAPE5` input files are read as formatted FORTRAN. As a consequence of the formatted read, any blank space will be read as "zero".  Thus, one may leave blanks for most of the parameters and within the code they will default to an acceptable value.

Real numbers format input as either `E` or `F` format, with the entire number within the range specified in the input instructions.  Integers must be specified exactly in the integer format.  For example, the spectral bandwidth (_v<sub>1</sub>_ to _v<sub>2</sub>_) in LBLRTM `TAPE5` is input as 10 character real numbers. This means that the value can be written anywhere within these 10 characters, as long as there is a decimal point (e.g. "---600.000" or "-600.000--", where "-" is a blank space).

Integers are read in with the `I` format.  For example, the model atmosphere (`iatm`) in LBLRTM `TAPE5` is input as `I5`, so it must be "----2", and not "2----" as this will be read as 20000.

## LBLRTM Output File Format <a name="lblout"></a>

The general structure of the files involves the use of panels, which are blocks of output usually containing 2400 points. Each panel contains a header to describe the starting and ending points of the panel (_v<sub>1</sub>_ and _v<sub>2</sub>_), the spectral spacing of the points (`dvp`), and the number of points in the panel (`npts`). The panel header is followed by either one or two (see below) blocks of output, consisting of `npts` points.

[this needs work]

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

# Instructions and Tips for Running LNFL <a name="runlnfl"></a>

LNFL is used to generate a unformatted file (`TAPE3`) of all the line parameters required by LBLRTM.

## Input files for LNFL <a name="lnflin"></a>

1. `TAPE1`: The line parameter database in ASCII format (downloaded with the [AER Line File repository](https://github.com/AER-RC/AER_Line_File) or from [Zenodo](https://zenodo.org/record/4019086)).

2. `TAPE5`: LNFL input file.

## Output files for LNFL <a name="lnflout"></a>

1. `TAPE3`: Unformatted LNFL output file containing the line parameters for LBLRTM.
2. `TAPE6`: Informational output file.
3. `TAPE7`: Optional output file containing ASCII version of the parameters contained in `TAPE3`.

## Sequence for running LNFL <a name="lnflseq"></a>

*	[Clone the latest LNFL code](https://github.com/AER-RC/LNFL#cloning-)  and download the latest line parameter database with the [AER Line File repository](https://github.com/AER-RC/AER_Line_File) or from [Zenodo](https://zenodo.org/record/4019086).
*	Compile LNFL using the makefiles found in the LNFL tar file.  Note: one needs to compile in the `build` directory.
*	Link the line parameter database to `TAPE1` in the LNFL working directory.
*	Remove `TAPE3` file from the LNFL working directory.
*	Edit necessary parameters in the `TAPE5` input file.  Note that the beginning and ending wavenumber (_v<sub>1</sub>_, _v<sub>2</sub>_) in `TAPE5` must extend at least 25 cm<sup>-1</sup> beyond each end of the desired spectral range for the LBLRTM calculations.
*	Run the LNFL executable.

# Instructions and Tips for Compiling and Running LBLRTM <a name="runlbl"></a>

LBLRTM is used to generate line-by-line upwelling and downwelling transmittances and radiances.

## Required input files for LBLRTM <a name="lblin"></a>
1. `TAPE3`: Unformatted file containing line parameter information, generated by LNFL (see above). The `TAPE3` file should include lines from at least 25 cm<sup>-1</sup> on either end of the calculation region.
2. `TAPE5`: Input file required to run LBLRTM.

The spectral interval (_v<sub>1</sub>_, _v<sub>2</sub>_) for any LBLRTM run must not exceed 2000 cm<sup>-1</sup> (see instruction manual).

Other input files are required if you are using the solar source function, cross sections, surface emissivity, etc. See the LBLRTM instruction manual and provide example.

## Layer numbering scheme <a name="laynum"></a>

The LBLRTM convention is that layer 1 is at the highest pressure level (lowest altitude).  The layer information for a given run may be found in `TAPE6`.

## Output files for LBLRTM <a name="lblout"></a>
1. `TAPE6`: Informational output file
2. `TAPE11`: Unformatted file containing filtered output, if requested in `TAPE5`.
3. `TAPE12`: Unformatted file containing transmittances/radiances.

ASCII file of unformatted unformatted files can be requested in the LBLRTM TAPE5 (see `pltlbl` variable in Record 12).

Unformatted optical depth files can be requested in the LBLRTM using options specified in `TAPE5`.

## Sequence for running LBLRTM <a name="lblseq"></a>
* [Clone the latest LBLRTM code](https://github.com/AER-RC/LBLRTM#cloning-the-latest-release-) and download the latest line parameter database with the [AER Line File repository](https://github.com/AER-RC/AER_Line_File) or from [Zenodo](https://zenodo.org/record/4019086).
* Compile LBLRTM following makefiles in the LBLRTM tar file. Note, one needs to compile in the `build` directory
* Link the line parameter database (`TAPE3` from LNFL) to the LBLRTM working directory.
* Edit any parameters necessary in the input file `TAPE5`.
* Run the LBLRTM executable.

# Tests <a name="tests"></a>

As of LBLRTM v12.10, a [run example package](https://github.com/AER-RC/LBLRTM/releases/download/v12.11/lblrtm_v12.11.examples.tar) is provided separately from the code repository. It can be used to validate building and running of the model for select atmospheric specifications and model configurations. See `README.setup` in top level of the package for further direction.

# Frequently Asked Questions <a name="faq"></a>

1. **What is the difference between a line-by-line calculation and a band-model calculation?**

Absorption/emission spectra are comprised of a complicated array of spectral lines.  The HITRAN 2008 Database (Version 13.0) contains over 2,713,000 lines for 39 different molecules. In order to resolve these individual lines, a nominal spectral sampling rate of less than the mean line half width must be utilized.  Such highly resolved radiative transfer calculations are called line-by-line (LBL) calculations. The computational time associated with calculating broadband fluxes from LBL calculations is formidable.  A band model aims to simplify radiative transfer calculations by using approximations to represent the line-by-line characteristics of a particular spectral interval.  Band models are appropriate for situations where the desired spectral resolution is much smaller than the Lorentz and Doppler widths of the spectral lines. Such approximations are also of use in general circulation models.

2. **What are the standard units used in LBLRTM calculations?**

| Output Variable | Units |
| :---: | :---: |
| Wavenumber | cm<sup>-1</sup> |
| Radiance | W cm<sup>-2</sup> sr<sup>-1</sup> /  cm-1 |
| Brightness Temperature | K |
| Analytic Jacobians (dR/dx) | <ul><li>molecules: W cm<sup>-2</sup> sr<sup>-1</sup> / cm<sup>-1</sup> / log(VMR)</li><li>temperature: W cm<sup>-2</sup> sr<sup>-1</sup> / cm-1  / K |

3. **Radiance Derivatives (Jacobians)**

This section describes the Analytical Jacobian capability in LBLRTM, Version 10.0 and later. The results from earlier versions are not reliable and should NOT be used.
The implementation of the analytic Jacobians in LBLRTM has been designed to require a minimal amount of setup on the part of the user while exploiting pre-existing LBLRTM calculation options. There are three steps required to obtain layer and level Analytic Jacobians:

  1. Create the `ODint_lll` files using `IMRG=1` and `IOD =3`; "lll" is the layer number.  Following LBLRTM convention, layer 1 is at the highest pressure level (lowest altitude). The spectral grid for the monochromatic calculation is determined by the DV of the highest layer (`DVSET` is automatically set by the program).

  2. Create the `RDDN_lll` files which provide the radiances from the top of the profile to level "lll" using `IMRG=40` and `IOD=3`.  Note that the `RDDNlayer_001` file is the downwelling radiance at the surface. The `ODint` files from Step 1 are used for this calculation.

  3. Create the layer and level Analytic Jacobian files in directory `AJ`; `RDderivUPW_xx_lll` and `LEV_RDderivUPW_xx_lll` are the layer and level derivatives files of the upwelling radiance at the upper boundary of the profile (`IMRG=41` and `IOD=3`) taken with respect to state parameter xx. Similarly, the files `RDderivDNW_xx_lll` and `LEV_RDderivDNW_xx_lll` are obtained for the downwelling radiance at the lower boundary (`IMRG=40` and `IOD=3`) with the Jacobian taken with respect to state vector
type xx where:

  * xx = -1 surface parameters
  * xx = 0 temperature
  * xx = mol molecule number (1-38)

Multiple runs of Step 3 may be performed once Step 1 and Step 2 have been run, the radiometric representation for the Jacobians having been established,

All the unformatted files resulting from these operations are fully consistent with LBLRTM files. The layer and level Jacobian files have the Radiance Jacobians in the normal radiance panel and the transmittance from the lowest boundary of the problem to the bottom level of the layer lll in the transmittance panel or to the designated level in the case of level Jacobians. These transmittances are included only to retain consistency with the LBLRTM file structures. Consequently, postprocessing including the application of selected instrument functions is accomplished in a manner identical to that for radiances. Note that selecting brightness temperature in the postprocessing will not provide a meaningful result.  In general, the input parameters for AJ calculations is described in the LBLRTM instructions where the required parameter is described, particularly note RECORD 1.2.b.  The `scanmrg` option (`IMRG=42,43` in this case) has not been tested and should be used with extreme caution.

Finally, a script has been included with ALL necessary files to run a sample set of Jacobian calculations. A second script has been provided to perform symmetric finite difference calculations to check the AJ results.  A PowerPoint document with plots showing the results from the two scripts has been included. The file `lblrtm_AJ_readme.txt` explains the two scripts.

4. **Is it possible to scale the profile of one or more species?**

Yes. The instructions for using this capability are provided in the lblrtm instructions. See records 1.3, 1.3.a and 1.3.b.

5. **Does LBLRTM include heavy molecule parameters (cross-sectional species)?**

Heavy molecules (such as CCL4, F11, and others listed in Table II of the lblrtm_instruction manual) can be included in LBLRTM calculations by setting the `IXSECT` input variable to 1 and adding Record 2.2 or Record 3.7 to the LBLRTM `TAPE5`  An additional file (`FSCDXS`) and directory (`xs`) are required for these calculations and can be obtained from the [cross sections repository](https://github.com/AER-RC/cross-sections) (which is cloned with LBLRTM if the directions in the [Cloning][#cloning] section are followed) or from the LBLRTM example tar file (available in the [LBLRTM v12.11 Release](https://github.com/AER-RC/LBLRTM/releases/download/v12.11/lblrtm_v12.11.examples.tar)).

6. **Format of external surface emissivity/reflectivity files**

Sea surface spectral emissivity and reflectivity files are provided with the example (available in the [LBLRTM v12.11 Release](https://github.com/AER-RC/LBLRTM/releases/download/v12.11/lblrtm_v12.11.examples.tar)). The files must have the file names of `EMISSIVITY` and `REFLECTIVITY`. The format is as follows:

| Parameter | Format | Description |
| :--- | :---: | :--- |
| `V1EMIS` | `E10.3` | Initial emissivity/reflectivity frequency value [cm<sup>-1</sup>] |
| `V2EMIS` | `E10.3` | Finial emissivity/reflectivity frequency value [cm<sup>-1</sup>] |
| `DVEMIS` | `E10.3` | Frequency Increment [cm<sup>-1</sup>] |
| `NLIMEM` | `I5` | Number of spectral emissivity/reflectivity points in the file |
| `ZEMIS` | `E15.7` | Emissivity at each spectral point |

**NOTE**: It is assumed that the spectral emissivity/reflectivity points are equally spaced and there is a maximum number of points (see instructions).

7. **Absorption due to clouds/aerosols and LOWTRAN5 routines**

Absorption due to clouds and aerosols can be computed in LBLRTM by setting the `IAERSL` flag in the input `TAPE5` file (refer to instructions).  This flag allows for LBLRTM to utilize the aerosol capabilities of LOWTRAN5.

8. **Solar Radiance**

Solar radiance calculations can be performed by utilizing LBLRTM input options such as `IEMIT=2` and a particular solar source function file `SOLAR.RAD`. A `SOLAR.RAD` file can be generated with program `extract_solar` available on the AER RT web site and the Kurucz solar source function.  The Kurucz solar source function has been used in AER’s research in shortwave radiation and is based on theoretical radiative transfer calculations for the solar atmosphere. The solar source function is available is at a high spectral resolution (i.e. for monochromatic calculations) and 1 cm<sup>-1</sup> resolution.

9. **Line coupling/mixing**

Line coupling parameters are utilized in LBLRTM for O<sub>2</sub>, CO<sub>2</sub> and CH<sub>4</sub>. The line coupling parameters are provided in the AER line parameter database (available in the [AER Line File repository](https://github.com/AER-RC/AER_Line_File) or on [Zenodo](https://zenodo.org/record/4019086)) and are written to the line parameter input file (`TAPE3`) by LNFL.

10. **What is the appropriate reference for LBLRTM calculations in journal articles and presentations?**

Clough SA, Shephard MW, Mlawer EJ, Delamere JS, Iacono MJ, Cady-Pereira K, Boukabara S, Brown PD.Atmospheric radiative transfer modeling: a summary of the AER codes • SHORT COMMUNICATION• _J Quant. Spectrosc. and Radiat Transfer_, **91**, 233-244 (2005).

Also, please refer to the [References Wiki Page](https://github.com/AER-RC/LBLRTM/wiki/References-and-Publications) for the complete list of references.

11. **How do you calculate fluxes?**

Source code and instructions available: [RADSUM](http://rtweb.aer.com/radsum_frame.html)
