# cgpwgs-nf

Nextflow version of dockstore-cgpwgs wrapper.

## How is this different

The original `dockstore-cgpwgs` codebase was designed to be run as a monolith process, all jobs on a single host managed
by a bas script.  This nextflow implementation allows each data type to be executed in isolation.  Additionally the individual
tools use their respective docker images, rather than relying on one image being up to date for each tool.

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
