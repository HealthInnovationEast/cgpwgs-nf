#!/usr/bin/env nextflow
nextflow.enable.dsl=2
/*
========================================================================================
                         lifebit-ai/xxxx
========================================================================================
lifebit-ai/xxxx
 #### Homepage / Documentation
https://github.com/xxxx
----------------------------------------------------------------------------------------
*/

// Help message

def helpMessage() {
    // TODO
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --bams sample.bam [Options]
    
    Inputs Options:
    --input         Input file
    Resource Options:
    --max_cpus      Maximum number of CPUs (int)
                    (default: $params.max_cpus)  
    --max_memory    Maximum memory (memory unit)
                    (default: $params.max_memory)
    --max_time      Maximum time (time unit)
                    (default: $params.max_time)
    See here for more info: https://github.com/lifebit-ai/hla/blob/master/docs/usage.md
    """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

/*--------------------------------------------------------
  Defining and showing header with all params information 
----------------------------------------------------------*/

// Header log info

def summary = [:]

if (workflow.revision) summary['Pipeline Release'] = workflow.revision

summary['Output dir']                                  = params.outdir
summary['Launch dir']                                  = workflow.launchDir
summary['Working dir']                                 = workflow.workDir
summary['Script dir']                                  = workflow.projectDir
summary['User']                                        = workflow.userName
// then arguments
summary['pairs']                                       = params.pairs
summary['genomefa']                                    = params.genomefa
summary['snpgc']                                       = params.snpgc
summary['sexloci']                                     = params.sexloci

log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Importantly, in order to successfully introspect:
// - This needs to be done first `main.nf`, before any (non-head) nodes are launched. 
// - All variables to be put into channels in order for them to be available later in `main.nf`.

ch_repository         = Channel.of(workflow.manifest.homePage)
ch_commitId           = Channel.of(workflow.commitId ?: "Not available is this execution mode. Please run 'nextflow run ${workflow.manifest.homePage} [...]' instead of 'nextflow run main.nf [...]'")
ch_revision           = Channel.of(workflow.manifest.version)

ch_scriptName         = Channel.of(workflow.scriptName)
ch_scriptFile         = Channel.of(workflow.scriptFile)
ch_projectDir         = Channel.of(workflow.projectDir)
ch_launchDir          = Channel.of(workflow.launchDir)
ch_workDir            = Channel.of(workflow.workDir)
ch_userName           = Channel.of(workflow.userName)
ch_commandLine        = Channel.of(workflow.commandLine)
ch_configFiles        = Channel.of(workflow.configFiles)
ch_profile            = Channel.of(workflow.profile)
ch_container          = Channel.of(workflow.container)
ch_containerEngine    = Channel.of(workflow.containerEngine)

/*----------------------------------------------------------------
  Setting up additional variables used for documentation purposes  
-------------------------------------------------------------------*/

Channel
    .of(params.raci_owner)
    .set { ch_raci_owner } 

Channel
    .of(params.domain_keywords)
    .set { ch_domain_keywords }

/*----------------------
  Setting up input data  
-------------------------*/

// Define Channels from input
// only if not in dsl2

/*-----------
  Processes  
--------------*/

// Do not delete this process
// Create introspection report

process obtain_pipeline_metadata {
    publishDir "${params.tracedir}", mode: "copy"

    input:
    val repository from ch_repository
    val commit from ch_commitId
    val revision from ch_revision
    val script_name from ch_scriptName
    val script_file from ch_scriptFile
    val project_dir from ch_projectDir
    val launch_dir from ch_launchDir
    val work_dir from ch_workDir
    val user_name from ch_userName
    val command_line from ch_commandLine
    val config_files from ch_configFiles
    val profile from ch_profile
    val container from ch_container
    val container_engine from ch_containerEngine
    val raci_owner from ch_raci_owner
    val domain_keywords from ch_domain_keywords

    output:
    file("pipeline_metadata_report.tsv") into ch_pipeline_metadata_report
    
    shell:
    '''
    echo "Repository\t!{repository}"                  > temp_report.tsv
    echo "Commit\t!{commit}"                         >> temp_report.tsv
    echo "Revision\t!{revision}"                     >> temp_report.tsv
    echo "Script name\t!{script_name}"               >> temp_report.tsv
    echo "Script file\t!{script_file}"               >> temp_report.tsv
    echo "Project directory\t!{project_dir}"         >> temp_report.tsv
    echo "Launch directory\t!{launch_dir}"           >> temp_report.tsv
    echo "Work directory\t!{work_dir}"               >> temp_report.tsv
    echo "User name\t!{user_name}"                   >> temp_report.tsv
    echo "Command line\t!{command_line}"             >> temp_report.tsv
    echo "Configuration file(s)\t!{config_files}"    >> temp_report.tsv
    echo "Profile\t!{profile}"                       >> temp_report.tsv
    echo "Container\t!{container}"                   >> temp_report.tsv
    echo "Container engine\t!{container_engine}"     >> temp_report.tsv
    echo "RACI owner\t!{raci_owner}"                 >> temp_report.tsv
    echo "Domain keywords\t!{domain_keywords}"       >> temp_report.tsv
    awk 'BEGIN{print "Metadata_variable\tValue"}{print}' OFS="\t" temp_report.tsv > pipeline_metadata_report.tsv
    '''
}

process ascat_counts {
    input:
        tuple path('genome.fa'), path('genome.fa.fai'), path('genome.fa.dict')
        path('snp.gc')
        path('sex.loci')
        tuple val(groupId), val(type), val(sampleId), val(protocol), val(platform), file(htsfile), file(htsidx), file(htsStats)

    output:
        tuple path("${sampleId}.count.gz"), path("${sampleId}.count.gz.tbi")
        path("${sampleId}.is_male.txt")
        tuple val(groupId), val(type), val(sampleId), val(protocol), val(platform), path("${sampleId}.count.gz"), path("${sampleId}.count.gz.tbi"), path("${sampleId}.is_male.txt"), emit: to_ascat

    // makes sure pipelines fail properly, plus errors and undef values
    shell = ['/bin/bash', '-euo', 'pipefail']

    script:
        """
        # remove logs for sucessful jobs
        export PCAP_THREADED_REM_LOGS=1
        ascatCounts.pl -o . \
            -b $htsfile \
            -r genome.fa \
            -sg snp.gc \
            -l sex.loci \
            -c $task.cpus
        """
}

process ascat {
    input:
        tuple path('genome.fa'), path('genome.fa.fai'), path('genome.fa.dict')
        path('snp.gc')
        tuple val(groupId), val(types), val(sampleIds), val(protocol), val(platform), path(counts), path(indexes), path(ismale)

    output:
        tuple path('*.copynumber.caveman.vcf.gz'), path('*.copynumber.caveman.vcf.gz.tbi')
        path('*.png')
        tuple val(groupId), path('*.copynumber.caveman.csv'), path('*.samplestatistics.txt'), emit: ascat_for_caveman
        path('*.copynumber.txt.gz')
        tuple path("*.count.gz", includeInputs: true), path("*.count.gz.tbi", includeInputs: true)
    
    publishDir {
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        "${params.outdir}/${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}/ascat"
    }, mode: 'copy'

    // makes sure pipelines fail properly, plus errors and undef values
    shell = ['/bin/bash', '-euo', 'pipefail']

    script:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        """
        # remove logs for sucessful jobs
        export PCAP_THREADED_REM_LOGS=1
        SPECIES=`head -n 2 genome.fa.dict | tail -n 1 | perl -ne 'm/SP:([^\t]+)/;print \$1;'`
        ASSEMBLY=`head -n 2 genome.fa.dict | tail -n 1 | perl -ne 'm/AS:([^\t]+)/;print \$1;'`
        ascat.pl -nb -f -o . \
            -ra "\$ASSEMBLY" -rs "\$SPECIES" \
            -pr "${protocol[ctrl_idx]}" -pl "${platform[ctrl_idx]}" \
            -r genome.fa \
            -sg snp.gc \
            -g ${ismale[ctrl_idx]} \
            -t ${counts[case_idx]} -tn ${sampleIds[case_idx]} \
            -n ${counts[ctrl_idx]} -nn ${sampleIds[ctrl_idx]} \
            -c $task.cpus
        """
}

process pindel {
    input:
        tuple path('genome.fa'), path('genome.fa.fai'), path('genome.fa.dict')
        tuple path('badloci.bed.gz'), path('badloci.bed.gz.tbi')
        tuple val(groupId), val(types), val(sampleIds), val(protocols), val(platforms), file(htsfiles), file(htsindexes), file(htsStats)
        // optional
        val exclude
        path(exclude_file)

    output:
        tuple val(groupId), path('*.vcf.gz'), path('*.vcf.gz.tbi'), emit: vcf
        path('*_mt.*')
        path('*_wt.*')
    
    publishDir {
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        "${params.outdir}/${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}/pindel"
    }, mode: 'copy'

    shell = ['/bin/bash', '-euo', 'pipefail']

    script:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        def apply_exclude = exclude != false ? "-e '$exclude'" : ''
        def apply_exclude_file = exclude_file.name != 'NO_FILE' ? "-ef $exclude_file" : ''
        """
        # remove logs for sucessful jobs
        export PCAP_THREADED_REM_LOGS=1
        SPECIES=`head -n 2 genome.fa.dict | tail -n 1 | perl -ne 'm/SP:([^\t]+)/;print \$1;'`
        ASSEMBLY=`head -n 2 genome.fa.dict | tail -n 1 | perl -ne 'm/AS:([^\t]+)/;print \$1;'`
        pindel.pl -noflag -o result \
        -sp "\$SPECIES" \
        -as "\$ASSEMBLY" \
        ${apply_exclude} \
        ${apply_exclude_file} \
        -r genome.fa \
        -t ${htsfiles[case_idx]} \
        -n ${htsfiles[ctrl_idx]} \
        -b badloci.bed.gz \
        -st ${protocols[case_idx]} \
        -c ${task.cpus}
        # easier to link the files than use "publishDir saveAs:"
        ln -f result/*.vcf.gz* .
        ln -f result/*_mt.* .
        ln -f result/*_wt.* .
        """
}

process pindel_flag {

    input:
        path(filter)
        val unmatched_ext
        tuple path('codingexon_regions.indel.bed.gz'), path('codingexon_regions.indel.bed.gz.tbi')
        tuple path("unmatched.${unmatched_ext}.gz"), path("unmatched.${unmatched_ext}.gz.tbi")
        tuple path('simplerepeats.bed.gz'), path('simplerepeats.bed.gz.tbi')
        tuple val(groupId), path('input.vcf.gz'), path('input.vcf.gz.tbi'), val(types), val(sampleIds)
        // optional
        path(softfil)
        val skipgerm
    
    output:
        tuple val(groupId), path('*.pindel.flagged.vcf.gz'), path('*.pindel.flagged.vcf.gz.tbi'), emit: flagged
        tuple val(groupId), path('*.pindel.germline.bed.gz'), path('*.pindel.germline.bed.gz.tbi'), emit: germline

    publishDir {
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        "${params.outdir}/${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}/pindel"
    }, mode: 'copy'

    shell = ['/bin/bash', '-euo', 'pipefail']

    script:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        def vs_filename = "${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}"
        def apply_soft = softfil.name != 'NO_FILE' ? "-sr $softfil" : ''
        def apply_skipgerm = skipgerm ? '-sg' : ''
        """
        FlagVcf.pl -r ${filter} -a codingexon_regions.indel.bed.gz -u unmatched.${unmatched_ext}.gz -s simplerepeats.bed.gz -i input.vcf.gz -o ${vs_filename}.pindel.flagged.vcf ${apply_soft} ${apply_skipgerm}
        bgzip -c ${vs_filename}.pindel.flagged.vcf >${vs_filename}.pindel.flagged.vcf.gz
        tabix -p vcf ${vs_filename}.pindel.flagged.vcf.gz
        
        GERM_FLAG=\$(grep F012 ${filter})
        GERMLINE_BED=${vs_filename}.pindel.germline.bed
        pindel_germ_bed.pl \
            -f \$GERM_FLAG \
            -i ${vs_filename}.pindel.flagged.vcf \
            -o \$GERMLINE_BED
        sort -k1,1 -k2,2n -k3,3n \$GERMLINE_BED | bgzip -c > \$GERMLINE_BED.gz
        tabix -p bed \$GERMLINE_BED.gz
        rm -f \$GERMLINE_BED
        """
}

process vagrent {
    input:
        tuple val(groupId), path(in_vcf), path(in_tbi), val(types), val(sampleIds)
        tuple path('vagrent.cache.gz'), path('vagrent.cache.gz.tbi'), path('vagrent.fa'), path('vagrent.fa.fai')

    output:
        path('*.annotated.vcf.gz*')

    publishDir {
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        "${params.outdir}/${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}/annotated"
    }, mode: 'copy'

    shell = ['/bin/bash', '-euo', 'pipefail']

    script:
        def annot = in_vcf.toString().minus('vcf.gz') + 'annotated.vcf'
        """
        AnnotateVcf.pl -t -i ${in_vcf} -o ${annot} -c vagrent.cache.gz
        """
}

process caveman {
    input:
        tuple path('genome.fa'), path('genome.fa.fai'), path('genome.fa.dict')
        tuple val(groupId), val(types), val(sampleIds), val(protocols), val(platforms), file(htsfiles), file(htsindexes), file(htsStats), path('copynumber.caveman.csv'), path('samplestatistics.txt')
        path('highDepth.tsv')
        val cavereads
        val exclude

    output:
        tuple val(groupId), path('*.muts.ids.vcf.gz'), path('*.muts.ids.vcf.gz.tbi'), emit: to_flag
        tuple path('*.snps.ids.vcf.gz'), path('*.snps.ids.vcf.gz.tbi')
        path('*.no_analysis.bed')
    
    publishDir {
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        "${params.outdir}/${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}/caveman"
    }, mode: 'copy'

    shell = ['/bin/bash', '-euo', 'pipefail']

    script:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        def apply_exclude = exclude != false ? "-x '$exclude'" : ''
        // -b $REF_BASE/caveman/flagging
        // -ab $REF_BASE/vagrent
        // -c $SNVFLAG \
        // -f $REF_BASE/caveman/flagging/flag.to.vcf.convert.ini \
        """
        # Get the species and assembly from the dict file
        SPECIES=`head -n 2 genome.fa.dict | tail -n 1 | perl -ne 'm/SP:([^\t]+)/;print \$1;'`
        ASSEMBLY=`head -n 2 genome.fa.dict | tail -n 1 | perl -ne 'm/AS:([^\t]+)/;print \$1;'`

        # prep ascat outputs
        NORM_CONTAM=`perl -ne 'if(m/^rho\s(.+)\n/){print 1-\$1;}' samplestatistics.txt`
        perl -ne '@F=(split q{,}, \$_)[1,2,3,4]; \$F[1]-1; print join("\t",@F)."\n";' < copynumber.caveman.csv > norm.cn.bed
        perl -ne '@F=(split q{,}, \$_)[1,2,3,6]; \$F[1]-1; print join("\t",@F)."\n";' < copynumber.caveman.csv > tum.cn.bed

        # remove logs for sucessful jobs under PCAP::Threaded module
        export PCAP_THREADED_REM_LOGS=1
        caveman.pl \
        -r genome.fa.fai \
        -ig highDepth.tsv \
        -u panel \
        -s "\$SPECIES" \
        -sa "\$ASSEMBLY" \
        -t ${task.cpus} \
        -st ${protocols[case_idx]} \
        -tc tum.cn.bed \
        -nc norm.cn.bed \
        -td 5 -nd 2 \
        -tb ${htsfiles[case_idx]} \
        -nb ${htsfiles[ctrl_idx]} \
        -e ${cavereads} \
        -o ./ \
        ${apply_exclude} \
        -k \$NORM_CONTAM \
        -no-flagging
        """
}

process caveman_vcf_split {
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        tuple val(groupId), path('input.vcf.gz'), path('input.vcf.gz.tbi')
        val cavevcfsplit
    
    output:
        tuple val(groupId), path('split.*'), emit: split_vcf
    
    script:
        """
        cgpVCFSplit.pl -i input.vcf.gz -o split -l ${cavevcfsplit}
        """
}

process caveman_flag {
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        tuple path('genome.fa'), path('genome.fa.fai'), path('genome.fa.dict')
        tuple val(groupId), val(types), val(sampleIds), val(protocols), val(platforms), file(htsfiles), file(htsindexes), file(htsStats), path('pindel.germline.bed.gz'), path('pindel.germline.bed.gz.tbi'), path(splitvcf)
        tuple path('panel/unmatchedNormal.bed.gz'), path('panel/unmatchedNormal.bed.gz.tbi')
        tuple path('bedfiles/centromeric_repeats.bed.gz'), path('bedfiles/centromeric_repeats.bed.gz.tbi')
        tuple path('bedfiles/hi_seq_depth.bed.gz'), path('bedfiles/hi_seq_depth.bed.gz.tbi')
        tuple path('bedfiles/simple_repeats.bed.gz'), path('bedfiles/simple_repeats.bed.gz.tbi')
        tuple path('bedfiles/snps.bed.gz'), path('bedfiles/snps.bed.gz.tbi')
        tuple path('flag.to.vcf.convert.ini'), path('flag.vcf.config.ini')
        tuple path('vagrent/codingexon_regions.sub.bed.gz'), path('vagrent/codingexon_regions.sub.bed.gz.tbi')
        tuple path('vagrent/gene_regions.bed.gz'), path('vagrent/gene_regions.bed.gz.tbi')
    
    output:
        tuple val(groupId), path('flagged.vcf'), emit: flagged
    
    script:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        """
        # Get the species and assembly from the dict file
        SPECIES=`head -n 2 genome.fa.dict | tail -n 1 | perl -ne 'm/SP:([^\t]+)/;print \$1;'`
        ASSEMBLY=`head -n 2 genome.fa.dict | tail -n 1 | perl -ne 'm/AS:([^\t]+)/;print \$1;'`
    
        cgpFlagCaVEMan.pl \
            -i $splitvcf \
            -o flagged.vcf \
            -s "\$SPECIES" \
            -sa "\$ASSEMBLY" \
            -m ${htsfiles[case_idx]} \
            -n ${htsfiles[ctrl_idx]} \
            -ref genome.fa.fai \
            -t ${protocols[case_idx]} \
            -b ./bedfiles \
            -g pindel.germline.bed.gz \
            -umv ./panel \
            -ab ./vagrent \
            -v flag.to.vcf.convert.ini \
            -c flag.vcf.config.ini
        """
}

process caveman_flag_merge {
    input:
        tuple val(groupId), path('?.vcf'), val(types), val(sampleIds)
    
    output:
        tuple val(groupId), path('*.snvs.flagged.vcf.gz'), path('*.snvs.flagged.vcf.gz.tbi'), emit: flagged
    
    publishDir {
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        "${params.outdir}/${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}/caveman"
    }, mode: 'copy'

    shell = ['/bin/bash', '-euo', 'pipefail']

    script:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        """
        FLAG_VCF=${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.snvs.flagged.vcf

        vcf-concat *.vcf | vcf-sort > \$FLAG_VCF
        bgzip -c \$FLAG_VCF > \$FLAG_VCF.gz
        tabix -p vcf \$FLAG_VCF.gz
        """
}

process brass {
    input:
        tuple path('genome.fa'), path('genome.fa.fai'), path('genome.fa.dict')
        tuple val(groupId), val(types), val(sampleIds), val(protocols), val(platforms), file(htsfiles), file(htsindexes), file(htsStats), file('ascat.samplestatistics.txt')
        tuple path('vagrent.cache.gz'), path('vagrent.cache.gz.tbi'), path('vagrent.fa'), path('vagrent.fa.fai')
        path('HiDepth.bed.gz')
        tuple path('brass_np.groups.gz'), path('brass_np.groups.gz.tbi')
        path('viral.genomic.fa.2bit')
        path('all_ncbi_bacteria.*')
        path('500bp_windows.gc.bed.gz')
        path('CentTelo.tsv')
        path('cytoband.txt')

    output:
        path('*.brm.bam')
        path('*.brm.bam.bai')
        path('*.vcf.gz')
        path('*.vcf.gz.tbi')
        path('*.bedpe.gz')
        path('*.bedpe.gz.tbi')
        path('*.ngscn.abs_cn.bw')
        path('*.intermediates.tar.gz')

    publishDir {
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        "${params.outdir}/${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}/brass"
    }, mode: 'copy'

    shell = ['/bin/bash', '-euo', 'pipefail']
    
    script:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        """
        # Get the species and assembly from the dict file
        SPECIES=`head -n 2 genome.fa.dict | tail -n 1 | perl -ne 'm/SP:([^\t]+)/;print \$1;'`
        ASSEMBLY=`head -n 2 genome.fa.dict | tail -n 1 | perl -ne 'm/AS:([^\t]+)/;print \$1;'`

        brass.pl -j 4 -k 4 -pl ILLUMINA \
            -c ${task.cpus} \
            -d HiDepth.bed.gz \
            -f brass_np.groups.gz \
            -g genome.fa \
            -s "\$SPECIES" -as "\$ASSEMBLY" \
            -pr ${protocols[case_idx]} \
            -g_cache vagrent.cache.gz \
            -vi viral.genomic.fa.2bit \
            -mi all_ncbi_bacteria \
            -b 500bp_windows.gc.bed.gz \
            -ct CentTelo.tsv \
            -cb cytoband.txt \
            -t ${htsfiles[case_idx]} \
            -n ${htsfiles[ctrl_idx]} \
            -ss ascat.samplestatistics.txt \
            -o ./
        """
}

process genotypes {
    input:
        tuple val(groupId), val(types), val(sampleIds), file(htsfiles), file(htsindexes)
        file('general.tsv')
        file('sex.tsv')
    
    output:
        file('*.tsv.gz')
        file('*.genotype.json')
        file('*.genotype.txt')
    
    
    publishDir {
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        "${params.outdir}/${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}/genotype"
    }, mode: 'copy'

    shell = ['/bin/bash', '-euo', 'pipefail']

    script:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        """
        compareBamGenotypes.pl \
            -o ./ \
            -j ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.genotype.json \
            -tb ${htsfiles[case_idx]} \
            -nb ${htsfiles[ctrl_idx]} \
            -s general.tsv \
            -g sex.tsv
        gzip ${sampleIds[case_idx]}.*.tsv
        gzip ${sampleIds[ctrl_idx]}.*.tsv
        """
}


workflow {
    genome = tuple file(params.genomefa), file("${params.genomefa}.fai"), file("${params.genomefa}.dict")
    badloci = tuple file(params.badloci), file("${params.badloci}.tbi")
    pairs = Channel.fromPath(params.pairs)
    case_control_map = pairs.splitCsv(header: true).map { row -> tuple(row.groupId, row.type, row.sampleId, row.protocol, row.platform, file(row.reads), file(row.readIdx), file(row.readStats)) }

    exclude_file = file(params.exfile)
    filter = file(params.filter)
    softfil = file(params.softfil)
    exons_indel = tuple file(params.exons_indel), file("${params.exons_indel}.tbi")
    simrep = tuple file(params.simrep), file("${params.simrep}.tbi")
    unmatched = tuple file(params.unmatched), file("${params.unmatched}.tbi")
    unmatched_ext = params.unmatched.contains('gff3') ? 'gff3' : 'bed'
    vagrentset = tuple file(params.vcache), file("${params.vcache}.tbi"), file(params.vfasta), file("${params.vfasta}.fai")
    softfil = file(params.softfil)

    snpgc = file(params.snpgc)
    sexloci = file(params.sexloci)
    genotypeloci = file(params.genotypeloci)

    // caveman specific
    cavehigh = file(params.cavehigh)

    // caveman flagging
    cavepanel     = tuple file(params.cavepanel), file("${params.cavepanel}.tbi")
    cave_centrep  = tuple file(params.cave_centrep), file("${params.cave_centrep}.tbi")
    cave_hibed    = tuple file(params.cave_hibed), file("${params.cave_hibed}.tbi")
    cave_simrep   = tuple file(params.cave_simrep), file("${params.cave_simrep}.tbi")
    cave_snps     = tuple file(params.cave_snps), file("${params.cave_snps}.tbi")
    cave_flag_ini = tuple file(params.cave_flag_to), file("${params.cave_flagvcf}")
    exons_snv     = tuple file(params.exons_snv), file("${params.exons_snv}.tbi")
    genes_snv     = tuple file(params.genes_snv), file("${params.genes_snv}.tbi")

    // brass
    viral      = file(params.viral)
    brass_np   = tuple file(params.brass_np), file("${params.brass_np}.tbi")
    bacterial  = file("${params.bacterial}.*") // should be 4
    gc_bins    = file(params.gc_bins)
    centtelo   = file(params.centtelo)
    cytoband   = file(params.cytoband)
    brass_high = file(params.brass_high)

    main:

        ascat_counts(
            genome,
            snpgc,
            sexloci,
            case_control_map
        )
        
        ascat(
            genome,
            snpgc,
            ascat_counts.out.to_ascat.groupTuple()
        )

        grouped_align = case_control_map.groupTuple()

        pindel(
            genome,
            badloci,
            grouped_align,
            params.exclude,
            exclude_file
        )

        type_samp = grouped_align.map { idx, type, samp, prot, plat, reads, ridx, bas -> [idx, type, samp] }

        pindel_flag(
            filter,
            unmatched_ext,
            exons_indel,
            unmatched,
            simrep,
            pindel.out.vcf.combine(type_samp, by: 0),
            softfil,
            params.skipgerm
        )

        align_and_ascat = grouped_align.combine(ascat.out.ascat_for_caveman, by: 0)

        caveman(
            genome,
            align_and_ascat,
            cavehigh,
            params.cavereads,
            params.exclude
        )

        caveman_vcf_split(
            caveman.out.to_flag,
            params.cavevcfsplit
        )

        // handle cases where split file list is a single element and currently unable to force to a list.
        cleaned_split = caveman_vcf_split.out.split_vcf.map {
            idx, maybe_list -> [idx, maybe_list instanceof Collection ? maybe_list : [maybe_list]]
        }
        flag_set = grouped_align.combine(pindel_flag.out.germline, by:0)
                    .combine(cleaned_split, by: 0)
                    .transpose(by: 10)

        caveman_flag(
            genome,
            flag_set,
            cavepanel,
            cave_centrep,
            cave_hibed,
            cave_simrep,
            cave_snps,
            cave_flag_ini,
            exons_snv,
            genes_snv
        )

        caveman_flag_merge(
            caveman_flag.out.flagged.groupTuple(by: 0).combine(type_samp, by: 0)
        )

        to_annotate = pindel_flag.out.flagged.concat(caveman_flag_merge.out.flagged)

        vagrent(
            to_annotate.combine(type_samp, by: 0),
            vagrentset
        )

        // just the bits we need staging for this process
        sample_stats = ascat.out.ascat_for_caveman.map { idx, counts, samp_stats -> [idx, samp_stats] }
        align_and_ss = grouped_align.combine(sample_stats, by: 0)

        brass(
            genome,
            align_and_ss,
            vagrentset,
            brass_high,
            brass_np,
            viral,
            bacterial,
            gc_bins,
            centtelo,
            cytoband
        )

        genotypes(
            grouped_align.map { idx, type, samp, prot, plat, reads, ridx, bas -> [idx, type, samp, reads, ridx] },
            genotypeloci,
            sexloci
        )



    // # 1 pair
    // rm -rf results/* .nextflow* work && nextflow -bg run /home/kr525/git/cynapse-ccri/cgpwgs-nf/main.nf --pairs /home/kr525/git/cynapse-ccri/cgpwgs-nf/test.csv --outdir results --genomefa data/cgpwgs_ref/GRCh37/genome.fa --snpgc data/cgpwgs_ref/GRCh37/ascat/SnpGcCorrections.tsv --sexloci data/cgpwgs_ref/GRCh37/gender.tsv --genotypeloci  data/cgpwgs_ref/GRCh37/general.tsv --badloci data/cgpwgs_ref/GRCh37/pindel/HiDepth.bed.gz --exfile data/cgpwgs_ref/GRCh37/pindel/exclude.lst --filter data/cgpwgs_ref/GRCh37/pindel/WGS_Rules.lst --exfile data/cgpwgs_ref/GRCh37/pindel/exclude.lst --simrep data/cgpwgs_ref/GRCh37/pindel/simpleRepeats.bed.gz --unmatched data/cgpwgs_ref/GRCh37/pindel/pindel_np.gff3.gz --softfil data/cgpwgs_ref/GRCh37/pindel/softRules.lst --vagrentroot data/cgpwgs_ref/GRCh37/vagrent --cavepanel data/cgpwgs_ref/GRCh37/caveman/unmatchedNormal.bed.gz --cavehigh data/cgpwgs_ref/GRCh37/caveman/HiDepth.tsv --exclude 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,X,Y,MT,NC_007605,hs37d5,GL% --cavebedroot data/cgpwgs_ref/GRCh37/caveman/flagging --brassroot data/cgpwgs_ref/GRCh37/brass -profile singularity >& nf.bg.log& 
    // # resume (don't delete everything)
    // nextflow -bg run /home/kr525/git/cynapse-ccri/cgpwgs-nf/main.nf --pairs /home/kr525/git/cynapse-ccri/cgpwgs-nf/test.csv --outdir results --genomefa data/cgpwgs_ref/GRCh37/genome.fa --snpgc data/cgpwgs_ref/GRCh37/ascat/SnpGcCorrections.tsv --sexloci data/cgpwgs_ref/GRCh37/gender.tsv --genotypeloci  data/cgpwgs_ref/GRCh37/general.tsv --badloci data/cgpwgs_ref/GRCh37/pindel/HiDepth.bed.gz --exfile data/cgpwgs_ref/GRCh37/pindel/exclude.lst --filter data/cgpwgs_ref/GRCh37/pindel/WGS_Rules.lst --exfile data/cgpwgs_ref/GRCh37/pindel/exclude.lst --simrep data/cgpwgs_ref/GRCh37/pindel/simpleRepeats.bed.gz --unmatched data/cgpwgs_ref/GRCh37/pindel/pindel_np.gff3.gz --softfil data/cgpwgs_ref/GRCh37/pindel/softRules.lst --vagrentroot data/cgpwgs_ref/GRCh37/vagrent --cavepanel data/cgpwgs_ref/GRCh37/caveman/unmatchedNormal.bed.gz --cavehigh data/cgpwgs_ref/GRCh37/caveman/HiDepth.tsv  --exclude 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,X,Y,MT,NC_007605,hs37d5,GL% --cavebedroot data/cgpwgs_ref/GRCh37/caveman/flagging --brassroot data/cgpwgs_ref/GRCh37/brass -profile singularity -resume >& nf.bg.log&
    // # no background/resume
    // nextflow run /home/kr525/git/cynapse-ccri/cgpwgs-nf/main.nf --pairs /home/kr525/git/cynapse-ccri/cgpwgs-nf/test.csv --outdir results --genomefa data/cgpwgs_ref/GRCh37/genome.fa --snpgc data/cgpwgs_ref/GRCh37/ascat/SnpGcCorrections.tsv --sexloci data/cgpwgs_ref/GRCh37/gender.tsv --genotypeloci  data/cgpwgs_ref/GRCh37/general.tsv --badloci data/cgpwgs_ref/GRCh37/pindel/HiDepth.bed.gz --exfile data/cgpwgs_ref/GRCh37/pindel/exclude.lst --filter data/cgpwgs_ref/GRCh37/pindel/WGS_Rules.lst --exfile data/cgpwgs_ref/GRCh37/pindel/exclude.lst --simrep data/cgpwgs_ref/GRCh37/pindel/simpleRepeats.bed.gz --unmatched data/cgpwgs_ref/GRCh37/pindel/pindel_np.gff3.gz --softfil data/cgpwgs_ref/GRCh37/pindel/softRules.lst --vagrentroot data/cgpwgs_ref/GRCh37/vagrent --cavepanel data/cgpwgs_ref/GRCh37/caveman/unmatchedNormal.bed.gz --cavehigh data/cgpwgs_ref/GRCh37/caveman/HiDepth.tsv  --exclude 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,X,Y,MT,NC_007605,hs37d5,GL% --cavebedroot data/cgpwgs_ref/GRCh37/caveman/flagging --brassroot data/cgpwgs_ref/GRCh37/brass -profile singularity -resume

    // # 2 pair
    // rm -rf results/* .nextflow* work && nextflow -bg run /home/kr525/git/cynapse-ccri/cgpwgs-nf/main.nf --pairs /home/kr525/git/cynapse-ccri/cgpwgs-nf/test_2.csv --outdir results --genomefa data/cgpwgs_ref/GRCh37/genome.fa --snpgc data/cgpwgs_ref/GRCh37/ascat/SnpGcCorrections.tsv --sexloci data/cgpwgs_ref/GRCh37/gender.tsv --genotypeloci  data/cgpwgs_ref/GRCh37/general.tsv --badloci data/cgpwgs_ref/GRCh37/pindel/HiDepth.bed.gz --exfile data/cgpwgs_ref/GRCh37/pindel/exclude.lst --filter data/cgpwgs_ref/GRCh37/pindel/WGS_Rules.lst --exfile data/cgpwgs_ref/GRCh37/pindel/exclude.lst --simrep data/cgpwgs_ref/GRCh37/pindel/simpleRepeats.bed.gz --unmatched data/cgpwgs_ref/GRCh37/pindel/pindel_np.gff3.gz --softfil data/cgpwgs_ref/GRCh37/pindel/softRules.lst --vagrentroot data/cgpwgs_ref/GRCh37/vagrent --cavepanel data/cgpwgs_ref/GRCh37/caveman/unmatchedNormal.bed.gz --cavehigh data/cgpwgs_ref/GRCh37/caveman/HiDepth.tsv --exclude 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,MT,NC_007605,hs37d5,GL% --cavebedroot data/cgpwgs_ref/GRCh37/caveman/flagging --brassroot data/cgpwgs_ref/GRCh37/brass -profile singularity >& nf.bg.log&
    // resume (don't delete everything)
    // nextflow -bg run /home/kr525/git/cynapse-ccri/cgpwgs-nf/main.nf --pairs /home/kr525/git/cynapse-ccri/cgpwgs-nf/test_2.csv --outdir results --genomefa data/cgpwgs_ref/GRCh37/genome.fa --snpgc data/cgpwgs_ref/GRCh37/ascat/SnpGcCorrections.tsv --sexloci data/cgpwgs_ref/GRCh37/gender.tsv --genotypeloci  data/cgpwgs_ref/GRCh37/general.tsv --badloci data/cgpwgs_ref/GRCh37/pindel/HiDepth.bed.gz --exfile data/cgpwgs_ref/GRCh37/pindel/exclude.lst --filter data/cgpwgs_ref/GRCh37/pindel/WGS_Rules.lst --exfile data/cgpwgs_ref/GRCh37/pindel/exclude.lst --simrep data/cgpwgs_ref/GRCh37/pindel/simpleRepeats.bed.gz --unmatched data/cgpwgs_ref/GRCh37/pindel/pindel_np.gff3.gz --softfil data/cgpwgs_ref/GRCh37/pindel/softRules.lst --vagrentroot data/cgpwgs_ref/GRCh37/vagrent --cavepanel data/cgpwgs_ref/GRCh37/caveman/unmatchedNormal.bed.gz --cavehigh data/cgpwgs_ref/GRCh37/caveman/HiDepth.tsv  --exclude 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,MT,NC_007605,hs37d5,GL% --cavebedroot data/cgpwgs_ref/GRCh37/caveman/flagging --brassroot data/cgpwgs_ref/GRCh37/brass -profile singularity -resume >& nf.bg.log&
    // # no background/resume
    // nextflow run /home/kr525/git/cynapse-ccri/cgpwgs-nf/main.nf --pairs /home/kr525/git/cynapse-ccri/cgpwgs-nf/test_2.csv --outdir results --genomefa data/cgpwgs_ref/GRCh37/genome.fa --snpgc data/cgpwgs_ref/GRCh37/ascat/SnpGcCorrections.tsv --sexloci data/cgpwgs_ref/GRCh37/gender.tsv --genotypeloci  data/cgpwgs_ref/GRCh37/general.tsv --badloci data/cgpwgs_ref/GRCh37/pindel/HiDepth.bed.gz --exfile data/cgpwgs_ref/GRCh37/pindel/exclude.lst --filter data/cgpwgs_ref/GRCh37/pindel/WGS_Rules.lst --exfile data/cgpwgs_ref/GRCh37/pindel/exclude.lst --simrep data/cgpwgs_ref/GRCh37/pindel/simpleRepeats.bed.gz --unmatched data/cgpwgs_ref/GRCh37/pindel/pindel_np.gff3.gz --softfil data/cgpwgs_ref/GRCh37/pindel/softRules.lst --vagrentroot data/cgpwgs_ref/GRCh37/vagrent --cavepanel data/cgpwgs_ref/GRCh37/caveman/unmatchedNormal.bed.gz --cavehigh data/cgpwgs_ref/GRCh37/caveman/HiDepth.tsv  --exclude 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,MT,NC_007605,hs37d5,GL% --cavebedroot data/cgpwgs_ref/GRCh37/caveman/flagging --brassroot data/cgpwgs_ref/GRCh37/brass -profile singularity -resume

}


// process report {
//     publishDir "${params.outdir}/MultiQC", mode: 'copy'
//     // anything written to MultiQC/multiqc_report.html is shown in the UI
//     // it doesn't have to be anything to do with MultiQC

//     input:
//     file(report_dir) from ch_report_dir
//     file(aggregate_output_dir) from ch_aggregate_output
    
//     output:
//     file "multiqc_report.html" into ch_multiqc_report

//     script:
//     """
//     cp -r ${report_dir}/* .
//     # convert from pdf to png
//     for f in \$(ls $aggregate_output_dir/*.pdf); do
//        pdftoppm \$f -png > $aggregate_output_dir/\$(echo \$(basename \$f | cut -d. -f1)).png
//     done
//     Rscript -e "rmarkdown::render('report.Rmd',params = list(aggregate_output_dir='$aggregate_output_dir'))"
//     mv report.html multiqc_report.html
//     """
// }


//bam_ch = Channel.fromPath(params.bam)
// ch_bam = Channel.value(file(params.bam))
// idx = params.bam
// idx = idx.take(idx.lastIndexOf('.')) + '.bai'
// ch_bai = Channel.value(file(idx))
// // combine into tuples
// ch_bam_bai = ch_bam.combine(ch_bai)


// range_lst = params.ranges.split(',').collect()
// ch_range = Channel.fromList(range_lst)

// process count_reads {
//     publishDir "${params.outdir}", mode: 'copy'

//     // always include for safety
//     shell = ['/bin/bash', '-euo', 'pipefail']
    
//     input:
//     tuple val(sampleId), path("in.${params.htstype}"), path("in.${params.htsidx}")



//     output:
//     file "${sampleId}.count" into ch_result

//     script:
//     """
//     samtools view -c in.bam > ${sampleId}.count
//     """
// }