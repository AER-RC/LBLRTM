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
