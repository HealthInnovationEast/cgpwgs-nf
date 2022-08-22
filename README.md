# cgpwgs-nf

Nextflow version of dockstore-cgpwgs wrapper.

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
   nextflow run cgpwgs-nf/main.nf
    -profile $PROFILES
    --core_ref $PATH_TO_REF/core_ref_GRCh37d5.tar.gz
    --snv_indel $PATH_TO_REF/SNV_INDEL_ref_GRCh37d5-fragment.tar.gz
    --cvn_sv $PATH_TO_REF/CNV_SV_ref_GRCh37d5_brass6+.tar.gz
    --annot $PATH_TO_REF/VAGrENT_ref_GRCh37d5_ensembl_75.tar.gz
    --qc_genotype $PATH_TO_REFs/qcGenotype_GRCh37d5.tar.gz
    --pairs $PATH_TO_UPDATED_CSV/test.csv
   ```

### Profiles

Ensure you set appropriate values for `-profile`.  For example, to use the `test` data on a `slurm` compute farm using
`singularity` containers the variable should be set as:

```
-profile test,slurm,singularity
```

You can see all available profiles under `./conf/`.

NOTE: do not use the `test` profile for anything other than testing with GRCh37 as it will mask the majority of the genome.
