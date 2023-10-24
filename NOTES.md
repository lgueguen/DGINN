# TODO

- allow to start the pipeline from intermediary step by providing input file. What about `Data` ?
- the possel/aln example has two input files
- the `all` rule should only have the last step files as output ?


# Running the pipeline


## Run with a local conda environment

To run locally, first create and activate a conda/mamba environment from environment.yml with:

```shell
mamba env create --file=environment.yml
mamba activate dginn
```

Then you can run the pipeline with:

```shell
snakemake --cores 1 --configfile config.yaml
```


## Run with docker

To run with docker:

```shell
docker build . -t dginn
docker run --rm -u $(id -u $USER) -v $(pwd):/opt/local dginn --cores 1 --configfile config.yaml
```


## Run with Singularity/Apptainer

To run with Singularity/Apptainer:

```shell
apptainer build dginn.sif Apptainer
apptainer run dginn.sif --cores 1 --configfile config.yaml
```

To export conda/mamba environment:

```shell
mamba env export > environment.yml
```


# Bio++ packages

The [Bio++](https://github.com/BioPP) packages in bioconda date back to the latest official release in 2018. To be able to used the more recent (if not yet unreleased) v3 version, it is needed to build conda packages locally.

Instructions follow those described in the bioconda documentation:

https://bioconda.github.io/contributor/building-locally.html


First, create a specific build environment:

```shell
mamba create -n dginn_build -c conda-forge -c bioconda bioconda-utils 
```

To build Bio++ conda packages locally, first verify the commit values in the recipes (files in `packages/recipes/bpp-*/meta.yaml`). Then build with:

```shell
cd packages
mamba activate dginn_build
# build all biopp packages
bioconda-utils build --docker --mulled-test recipes/
# build only biopp-core package
bioconda-utils build --docker --mulled-test --packages biopp-core recipes/
mamba deactivate
```

Verify that builds are ok:
```shell
mamba activate dginn_build
mamba search -c ${CONDA_PREFIX}/conda-bld bppsuite
mamba deactivate
```

To install them in your own local environment, and assuming that your `dginn` and `dginn_build` environements are in the same folder, you can use something like:
```shell
mamba activate dginn
mamba install -c ${CONDA_PREFIX}/../dginn_build/conda-bld bppsuite
```