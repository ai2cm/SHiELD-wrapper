# SHiELD-wrapper

This repository contains code for a Python-wrapped version of the public
version of GFDL's SHiELD model.  It is based on [AI2's existing Python-wrapped
version of
FV3GFS](https://github.com/ai2cm/fv3gfs-fortran/tree/master/FV3/wrapper), which
was described in [McGibbon et
al. (2021)](https://gmd.copernicus.org/articles/14/4401/2021/gmd-14-4401-2021.html),
and used in several published hybrid machine learning projects.  Much of the
code is identical or nearly identical to this original wrapper; the main
difference is that the wrapper in this repository wraps the latest public
version of GFDL's SHiELD model instead of an old fork of FV3GFS.

## Development

Development and testing of SHiELD-wrapper is currently only supported within a
Docker container.  Similar to how things are done in AI2's fv3net repository,
tests are set up to be run from within the Docker image defined in the
Dockerfile within this repository.  To build the docker image locally, make
sure to first initialize all the submodules, and then run `make build`:

```
$ make update_submodules
$ make build
```

Once the Docker image has been built, you can enter it with the wrapper testing
code editable from outside the container with:

```
$ make enter
```

To run the tests, you can use the following:

```
$ make test
```

To reset the regression test checksums, use:

```
$ make test_regtest_reset
```

The Python dependencies used in building and testing the wrapper are defined in
the `requirements.in` file.  The full set of pinned dependencies these resolve
into are defined in `requirements.txt`, which is produced using `pip-compile`.
