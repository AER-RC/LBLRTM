# Local Version

The `master` branch of hosted in GitHub *is* the local version of LBLRTM. Releases are in "Releases" and have tags. At AER, we have a single clone accessible to everyone for the purpose of having builds to which we can all point (`/nas/project/rc/rc1/lblrtm_local_version`, with PGI and Intel builds), but this clone IS NOT FOR DEVELOPMENT! For developing the model, please refer to the [Development section](#Dev).

# Development <a href="Dev"></a>

"Incoming" is no more. If you have contrbutions to make, branch off of master, make your edits, push to your branch, then create a merge request as soon as your commit has been pushed. Owners and/or maintainers will then execute the merge if it is acceptable. The development steps are (new users must generate an [SSH Key](https://docs.github.com/en/github/authenticating-to-github/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent#generating-a-new-ssh-key) and [add it to their GitHub account](https://docs.github.com/en/github/authenticating-to-github/adding-a-new-ssh-key-to-your-github-account)):

```
git clone --recursive git@github.com:AER-RC/LBLRTM.git
cd LBLRTM
git branch your_branch
git checkout your_branch
... (file modifications)
git commit -a -m 'made some changes to LBLRTM' # be more specific than this
git push origin your_branch
```

Merging and pulling should be done in the [web interface](https://github.com/AER-RC/LBLRTM/branches).

# Documentation

Document this model at will. `FAQ_LBLRTM.doc` has been removed from version control (and its PDF cousin will be gone, too) and converted to Markdown format and is now the top-level `README.md`. This README is the first place we should consider placing general information. Things for specific projects or special cases should have their own page in the [Wiki](https://github.com/AER-RC/LBLRTM/wiki).

# Release Procedure

A few of our repositories are dependent on each other. Contributors are encouraged to peruse the Wiki pages of [LBLRTM](https://github.com/AER-RC/LBLRTM/wiki/LBLRTM-New-Release-Procedure), [aer_rt_utils](https://github.com/AER-RC/aer_rt_utils/wiki/New-Release-Procedure), and [cross-section](https://github.com/AER-RC/cross-sections/wiki/New-Release-Procedure) for guidance on new releases of each repository and its associated submodules. Essentially, **each release of code base should be linked to specific releases of its submodules, and not the `master` branches** (assuming we know the submodule releases work with the given repository).

Releases should be done in the following order (submodules in parentheses):

1. [aer_rt_utils](https://github.com/AER-RC/aer_rt_utils/wiki/New-Release-Procedure)
2. [cross-sections](https://github.com/AER-RC/cross-sections/wiki/New-Release-Procedure)
3. [LNFL](https://github.com/AER-RC/LNFL/wiki/LNFL-New-Release-Procedure) (aer_rt_utils)
4. [LBLRTM](https://github.com/AER-RC/LBLRTM/wiki/LBLRTM-New-Release-Procedure) (aer_rt_utils, cross-sections)
5. [MT_CKD](https://github.com/AER-RC/MT_CKD/wiki/MT_CKD-New-Release-Procedure) (aer_rt_utils, LBLRTM)

Since LNFL and LBLRTM are not dependent on each other, the order of these two releases does not matter.
