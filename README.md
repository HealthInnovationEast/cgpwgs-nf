# cgpwgs-nf

Nextflow version of [dockstore-cgpwgs][ds-cgpwgs] wrapper, tailored for multi pair execution.

[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit)](https://github.com/pre-commit/pre-commit)

|                Main                |               Develop               |
| :--------------------------------: | :---------------------------------: |
| [![Main][gha-main]][gha-main-view] | [![Develop][gha-dev]][gha-dev-view] |

## How is this different

The original `dockstore-cgpwgs` codebase was designed to be run as a monolith process, all jobs on a single host managed
by a bash script.  This nextflow implementation allows for a more flexible execution across multiple hosts.  Additionally
the individual tools use their respective docker images (where they exist), rather than relying on a single image being
up to date for each tool.

In some cases additional optimisation has been possible over the original script.

## Reference files

### GRCh37

```
mkdir -p GRCh37/archives
cd GRCh37/archives
wget ftp.sanger.ac.uk/pub/cancer/dockstore/human/{core_ref_GRCh37d5.tar.gz,VAGrENT_ref_GRCh37d5_ensembl_75.tar.gz,CNV_SV_ref_GRCh37d5_brass6+.tar.gz
,SNV_INDEL_ref_GRCh37d5-fragment.tar.gz,qcGenotype_GRCh37d5.tar.gz}
```

### GRCh38

```
wget ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/{core_ref_GRCh38_hla_decoy_ebv.tar.gz,VAGrENT_ref_GRCh38_hla_decoy_ebv_ensembl_91.tar.gz,CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+.tar.gz,qcGenotype_GRCh38_hla_decoy_ebv.tar.gz,SNV_INDEL_ref_GRCh38_hla_decoy_ebv-fragment.tar.gz}
```

## Test data

A **GRCh37** mapped paired cell line is available from the original repository owners:

```
# case
wget http://ngs.sanger.ac.uk/production/cancer/dockstore/cgpwgs/sampled/COLO-829.{bam,bam.bai,bam.bas}
# control
wget http://ngs.sanger.ac.uk/production/cancer/dockstore/cgpwgs/sampled/COLO-829-BL.{bam,bam.bai,bam.bas}
```

## Execution

To execute with test data:

1. Download the *GRCh37* reference sets

1. Download the test data files

1. Copy `data/test.csv` and replace instances of `TEST_DATA_PATH` within the file with the path to the data files.

1. Run the following command providing suitable values for `PROFILES`, `PATH_TO_REF` and `PATH_TO_UPDATED_CSV`:

   ```bash
   nextflow run -r main https://github.com/cynapse-ccri/cgpwgs-nf \
    -profile $PROFILES \
    --core_ref $PATH_TO_REF/core_ref_GRCh37d5.tar.gz \
    --snv_indel $PATH_TO_REF/SNV_INDEL_ref_GRCh37d5-fragment.tar.gz \
    --cnv_sv $PATH_TO_REF/CNV_SV_ref_GRCh37d5_brass6+.tar.gz \
    --annot $PATH_TO_REF/VAGrENT_ref_GRCh37d5_ensembl_75.tar.gz \
    --qc_genotype $PATH_TO_REFs/qcGenotype_GRCh37d5.tar.gz \
    --pairs $PATH_TO_UPDATED_CSV/test.csv
   ```

## Development testing

Use the Nextflow `-stub-run` option with the `nextflow.stubRub.config`.

```
nextflow -c nextflow.stubRun.config run main.nf \
            -profile test -stub-run \
            ...
```

Once CYNAPSE is able to fully support the `-stub-run` option the additional config file will not be necessary.

### Profiles

Ensure you set appropriate values for `-profile`.  For example, to use the `test` data on a `slurm` compute farm using
`singularity` containers the variable should be set as:

```
-profile test,slurm,singularity
```

You can see all available profiles under `./conf/`.

NOTE: do not use the `test` profile for anything other than testing with GRCh37 as it will mask the majority of the genome.

#### CYNAPSE

On CYNAPSE you need to select 2 profiles:

1. `awsbatch`
1. `cynapse-pro-admin` or `cynapse-pro-wrkspc`
   - depending on use of admin or standard workspace - different queues

### Versions

|    Workflow     |                                                                                                                                               Images |
| :-------------: | ---------------------------------------------------------------------------------------------------------------------------------------------------: |
| 0.1.0 - current | ascatngs:4.5.0<br>brass:v6.3.4<br>cgpcavemanwrapper:1.18.2<br>cgppindel:3.10.0<br>dockstore-cgpwgs:2.1.1<br>dockstore-cgpwgs:2.1.1<br>vagrent:v3.7.0 |

<!-- refs -->

[ds-cgpwgs]: https://github.com/cancerit/dockstore-cgpwgs
[gha-dev]: https://github.com/cynapse-ccri/cgpwgs-nf/actions/workflows/build.yaml/badge.svg?branch=develop
[gha-dev-view]: https://github.com/cynapse-ccri/cgpwgs-nf/actions?query=branch%3Adevelop
[gha-main]: https://github.com/cynapse-ccri/cgpwgs-nf/actions/workflows/build.yaml/badge.svg?branch=main
[gha-main-view]: https://github.com/cynapse-ccri/cgpwgs-nf/actions?query=branch%3Amain
