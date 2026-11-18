# Publishing SprayNPray through Bioconda

The V2/cohort tools currently run from the repository checkout. The supplied
`environment.yml` installs their dependencies. This document describes the
remaining work to publish SprayNPray itself as a package named `spraynpray`;
adding an environment file or merging this PR does not publish that package.

## 1. Prepare an installable release

- Add Python packaging metadata (`pyproject.toml`) with the project name,
  a consistent distribution version, supported Python version, authors,
  homepage and license metadata matching `LICENSE`.
- Install the V2, comparison, wrapper and helper modules together. The wrapper
  currently locates sibling script files through `__file__` and launches them
  using the active Python interpreter, so packaging must preserve that layout
  or replace the subprocess calls with installed module entry points.
- Add executable entry points. Proposed command names are `spraynpray` for the
  genome-folder wrapper, `spraynpray-v2`, `spraynpray-compare` and
  `spraynpray-helper`. These are proposed packaging interfaces, not commands
  installed by the current checkout/environment file.
- Keep the original workflow available. If the original `spray-and-pray.py`
  is included, also install its auxiliary scripts/data and declare its own
  additional dependencies; test those resources from outside the checkout.
- Test a built distribution in a fresh environment and a different working
  directory: help commands, a small local-reference workflow, output creation
  and a second run that reuses completed alignments/ANI.

The V2/cohort runtime dependencies are Python, Prodigal, DIAMOND, FastANI,
NCBI Datasets CLI, NumPy, SciPy and Matplotlib. Record supported constraints
based on testing, distinguishing Python package requirements from external
executables. Do not bundle the large protein database, genome references,
protein-to-genome map or taxdump; users configure these as analysis inputs.

After the code and installation checks pass, merge the release changes and
publish a versioned GitHub release/tag. Bioconda can build directly from a
GitHub release source archive; a PyPI publication is not required.

## 2. Add a Bioconda recipe

Check package-name availability, then fork
[bioconda/bioconda-recipes](https://github.com/bioconda/bioconda-recipes)
and add `recipes/spraynpray/meta.yaml`. Include:

- Package name `spraynpray` and the released version.
- The versioned source archive URL and its actual SHA-256 checksum. Do not
  build from a moving `master` URL or leave a placeholder checksum.
- Build instructions using the Python packaging metadata, build number `0`
  for a new version, and appropriate build/host dependencies.
- Runtime requirements, including the external tool packages listed above.
- Entry-point/import checks and a small offline functional test. Large public
  reference downloads should not be required for recipe tests.
- Homepage, summary, description, redistribution-compatible license metadata,
  `license_file: LICENSE` and recipe maintainers.

For the pure Python V2/cohort distribution, `noarch: python` is normally
appropriate. Availability on a particular platform still depends on all its
external runtime packages. Including the original auxiliary resources may
require additional packaging work and dependencies.

Use the current [Bioconda contributor workflow](https://bioconda.github.io/contributor/index.html)
and [recipe guidelines](https://bioconda.github.io/contributor/guidelines.html)
to lint, build and test the recipe, then submit a pull request to the recipes
repository. Fix CI/review findings. After acceptance and successful publication,
the package becomes available from the Bioconda channel. Subsequent releases
need recipe version/checksum updates and validation.

## 3. User installation after publication

With conda-forge before Bioconda and strict channel priority already configured:

```bash
mamba install spraynpray
```

An explicit-channel installation into a new environment would be:

```bash
mamba create -n spraynpray --override-channels \
  -c conda-forge -c bioconda --strict-channel-priority spraynpray
conda activate spraynpray
```

The channel order follows [Bioconda's installation guidance](https://bioconda.github.io/).
Until the package has been published, use the source checkout and
`environment.yml` installation in [README.md](README.md).
