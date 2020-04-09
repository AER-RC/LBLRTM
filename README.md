# LBLRTM
Line-By-Line Radiatve Transfer Model by Atmospheric and Environmental Research

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

in the `LBLRTM` directory
