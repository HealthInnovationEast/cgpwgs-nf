# cgpwgs-nf

Nextflow version of dockstore-cgpwgs wrapper.

## How is this different

The original `dockstore-cgpwgs` codebase was designed to be run as a monolith process, all jobs on a single host managed
by a bas script.  This nextflow implementation allows each data type to be executed in isolation.  Additionally the individual
tools use their respective docker images, rather than relying on one image being up to date for each tool.

In some cases additional optimisation has been possible over the original script.