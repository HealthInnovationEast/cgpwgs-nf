#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def helpMessage() {
    // TODO
    log.info """
    Please see here for usage information: https://github.com/cynapse-ccri/cgpwgs-nf/blob/main/README.md#execution
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
summary['core_ref']                                    = params.core_ref
summary['snv_indel']                                   = params.snv_indel
summary['cnv_sv']                                      = params.cnv_sv
summary['annot']                                       = params.annot
summary['qc_genotype']                                 = params.qc_genotype
summary['exclude']                                     = params.exclude
summary['exfile']                                      = params.exfile
summary['skipgerm']                                    = params.skipgerm
summary['cavereads']                                   = params.cavereads
summary['cavevcfsplit']                                = params.cavevcfsplit
summary['cpus_caveman']                                = params.cpus_caveman
summary['cpus_pindel']                                 = params.cpus_pindel
summary['cpus_brass']                                  = params.cpus_brass
summary['cpus_counts']                                 = params.cpus_counts

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
      val(repository)
      val(commit)
      val(revision)
      val(script_name)
      val(script_file)
      val(project_dir)
      val(launch_dir)
      val(work_dir)
      val(user_name)
      val(command_line)
      val(config_files)
      val(profile)
      val(container)
      val(container_engine)
      val(raci_owner)
      val(domain_keywords)

    output:
      path("pipeline_metadata_report.tsv"), emit: pipeline_metadata_report

    // same as script except ! instead of $ for variables
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

    stub:
      '''
      touch pipeline_metadata_report.tsv
      '''
}

process prep_ref {
    input:
        file(core_ref)
        file(snv_indel)
        file(cnv_sv)
        file(annot)
        file(qc_genotype)

    output:
        path 'ref', type: 'dir', emit: ref
        path 'caveman', type: 'dir', emit: cave_flag
        path 'caveman_HiDepth.tsv', emit: cave_hidepth
        path 'pindel', type: 'dir', emit: pindel_flag
        tuple path('pindel_HiDepth.bed.gz'), path('pindel_HiDepth.bed.gz.tbi'), emit: pindel_hidepth
        path 'brass', type: 'dir', emit: brass
        path 'ascat/SnpGcCorrections.tsv', emit: snps_gc
        path 'general.tsv', emit: snps_genotype
        path 'gender.tsv', emit: snps_sex
        path 'verifyBamID_snps.vcf.gz' , emit: snps_verifyBamID
        path'vagrent', type: 'dir', emit: vagrent
        path'ref_cache', type: 'dir', emit: ref_cache

    shell = ['/bin/bash', '-euo', 'pipefail']

    stub:
        """
        mkdir -p ref
        mkdir -p caveman
        touch caveman_HiDepth.tsv
        mkdir -p pindel
        touch pindel_HiDepth.bed.gz
        touch pindel_HiDepth.bed.gz.tbi
        mkdir -p brass
        mkdir -p ascat
        touch ascat/SnpGcCorrections.tsv
        touch general.tsv
        touch gender.tsv
        touch verifyBamID_snps.vcf.gz
        mkdir -p vagrent
        mkdir -p ref_cache
        """

    script:
        """
        mkdir ref
        tar --strip-components 1 -C ref -zxvf $core_ref
        tar --strip-components 1 -zxvf $snv_indel
        tar --strip-components 1 -zxvf $cnv_sv
        tar --strip-components 1 -zxvf $annot
        tar --strip-components 1 -zxvf $qc_genotype

        # organise some elements further
        mv caveman/HiDepth.tsv caveman_HiDepth.tsv
        mv pindel/HiDepth.bed.gz pindel_HiDepth.bed.gz
        mv pindel/HiDepth.bed.gz.tbi pindel_HiDepth.bed.gz.tbi

        # build ref-cache
        seq_cache_populate.pl -subdirs 2 -root ./ref_cache ref/genome.fa
        """
}

process genotypes {
    input:
        tuple val(groupId), val(types), val(sampleIds), file(htsfiles), file(htsindexes)
        file('general.tsv')
        file('sex.tsv')
        path('ref_cache')

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

    stub:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        """
        touch ${sampleIds[case_idx]}.tsv.gz
        touch ${sampleIds[ctrl_idx]}.tsv.gz
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.genotype.json
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.genotype.txt
        """

    script:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        """
        export REF_CACHE=\$PWD/ref_cache/%2s/%2s/%s
        export REF_PATH=\$REF_CACHE

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

process ascat_counts {
    input:
        path('ref')
        path('snp.gc')
        path('sex.loci')
        tuple val(groupId), val(type), val(sampleId), val(protocol), val(platform), file(htsfile), file(htsidx), file(htsStats)

    output:
        tuple path("${sampleId}.count.gz"), path("${sampleId}.count.gz.tbi")
        path("${sampleId}.is_male.txt")
        tuple val(groupId), val(type), val(sampleId), val(protocol), val(platform), path("${sampleId}.count.gz"), path("${sampleId}.count.gz.tbi"), path("${sampleId}.is_male.txt"), emit: to_ascat

    // makes sure pipelines fail properly, plus errors and undef values
    shell = ['/bin/bash', '-euo', 'pipefail']

    stub:
        """
        touch ${sampleId}.count.gz
        touch ${sampleId}.count.gz.tbi
        touch ${sampleId}.is_male.txt
        """

    script:
        """
        # remove logs for sucessful jobs
        export PCAP_THREADED_REM_LOGS=1
        ascatCounts.pl -o . \
            -b $htsfile \
            -r ref/genome.fa \
            -sg snp.gc \
            -l sex.loci \
            -c $task.cpus
        """
}

process ascat {
    input:
        path('ref')
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

    stub:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        """
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.copynumber.caveman.vcf.gz
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.copynumber.caveman.vcf.gz.tbi
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.png
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.copynumber.caveman.csv
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.samplestatistics.txt
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.copynumber.txt.gz
        # others are from inputs
        """

    script:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        """
        # remove logs for sucessful jobs
        export PCAP_THREADED_REM_LOGS=1
        SPECIES=`head -n 2 ref/genome.fa.dict | tail -n 1 | perl -ne 'm/SP:([^\t]+)/;print \$1;'`
        ASSEMBLY=`head -n 2 ref/genome.fa.dict | tail -n 1 | perl -ne 'm/AS:([^\t]+)/;print \$1;'`
        ascat.pl -nb -f -o . \
            -ra "\$ASSEMBLY" -rs "\$SPECIES" \
            -pr "${protocol[ctrl_idx]}" -pl "${platform[ctrl_idx]}" \
            -r ref/genome.fa \
            -sg snp.gc \
            -g ${ismale[ctrl_idx]} \
            -t ${counts[case_idx]} -tn ${sampleIds[case_idx]} \
            -n ${counts[ctrl_idx]} -nn ${sampleIds[ctrl_idx]} \
            -c $task.cpus
        """
}

process pindel {
    input:
        path('ref')
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

    stub:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        """
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.pindel.vcf.gz
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.pindel.vcf.gz.tbi
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}_mt.pindel.bam
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}_mt.pindel.bam.bai
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}_mt.pindel.bw
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}_wt.pindel.bam
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}_wt.pindel.bam.bai
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}_wt.pindel.bw
        """

    script:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        def apply_exclude = exclude != false ? "-e '$exclude'" : ''
        def apply_exclude_file = exclude_file.name != 'NO_FILE' ? "-ef $exclude_file" : ''
        """
        # remove logs for sucessful jobs
        export PCAP_THREADED_REM_LOGS=1
        SPECIES=`head -n 2 ref/genome.fa.dict | tail -n 1 | perl -ne 'm/SP:([^\t]+)/;print \$1;'`
        ASSEMBLY=`head -n 2 ref/genome.fa.dict | tail -n 1 | perl -ne 'm/AS:([^\t]+)/;print \$1;'`
        pindel.pl -noflag -o result \
        -sp "\$SPECIES" \
        -as "\$ASSEMBLY" \
        ${apply_exclude} \
        ${apply_exclude_file} \
        -r ref/genome.fa \
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
        path('pindel')
        path('vagrent')
        tuple val(groupId), path('input.vcf.gz'), path('input.vcf.gz.tbi'), val(types), val(sampleIds), val(protocols)
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

    stub:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        """
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.pindel.flagged.vcf.gz
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.pindel.flagged.vcf.gz.tbi
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.pindel.germline.bed.gz
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.pindel.germline.bed.gz.tbi
        """

    script:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        def vs_filename = "${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}"
        def apply_skipgerm = skipgerm ? '-sg' : ''

        def softfile = 'pindel/softRules.lst'
        def filter = "pindel/${protocols[case_idx]}_Rules.lst"

        """
        SOFTRULES=''
        if [ -e ${softfile} ]; then
            SOFTRULES="-sr ${softfile}"
        fi

        FlagVcf.pl -r ${filter} \
            -a vagrent/codingexon_regions.indel.bed.gz \
            -u pindel/pindel_np.*.gz \
            -s pindel/simpleRepeats.bed.gz \
            -i input.vcf.gz \
            -o ${vs_filename}.pindel.flagged.vcf \
            \$SOFTRULES ${apply_skipgerm}
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

process caveman {
    input:
        path('ref')
        tuple val(groupId), val(types), val(sampleIds), val(protocols), val(platforms), file(htsfiles), file(htsindexes), file(htsStats), path('copynumber.caveman.csv'), path('samplestatistics.txt')
        path('highDepth.tsv')
        val cavereads
        val exclude
        val cavevcfsplit

    output:
        tuple val(groupId), path('*.muts.ids.vcf.gz'), path('*.muts.ids.vcf.gz.tbi')
        tuple path('*.snps.ids.vcf.gz'), path('*.snps.ids.vcf.gz.tbi')
        path('*.no_analysis.bed')
        tuple val(groupId), path('split.*'), emit: split_vcf

    publishDir {
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        "${params.outdir}/${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}/caveman"
    }, mode: 'copy', pattern: '*.{vcf.gz,vcf.gz.tbi,bed}'

    shell = ['/bin/bash', '-euo', 'pipefail']

    stub:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        """
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.muts.ids.vcf.gz
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.muts.ids.vcf.gz.tbi
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.snps.ids.vcf.gz
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.snps.ids.vcf.gz.tbi
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.no_analysis.bed
        touch split.1 split.2
        """

    script:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        def apply_exclude = exclude != false ? "-x '$exclude'" : ''
        """
        # Get the species and assembly from the dict file
        SPECIES=`head -n 2 ref/genome.fa.dict | tail -n 1 | perl -ne 'm/SP:([^\t]+)/;print \$1;'`
        ASSEMBLY=`head -n 2 ref/genome.fa.dict | tail -n 1 | perl -ne 'm/AS:([^\t]+)/;print \$1;'`

        # prep ascat outputs
        NORM_CONTAM=`perl -ne 'if(m/^rho\s(.+)\n/){print 1-\$1;}' samplestatistics.txt`
        perl -ne '@F=(split q{,}, \$_)[1,2,3,4]; \$F[1]-1; print join("\t",@F)."\n";' < copynumber.caveman.csv > norm.cn.bed
        perl -ne '@F=(split q{,}, \$_)[1,2,3,6]; \$F[1]-1; print join("\t",@F)."\n";' < copynumber.caveman.csv > tum.cn.bed

        # remove logs for sucessful jobs under PCAP::Threaded module
        export PCAP_THREADED_REM_LOGS=1
        caveman.pl \
        -r ref/genome.fa.fai \
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

        # split ready for flagging
        cgpVCFSplit.pl -i ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.muts.ids.vcf.gz -o split -l ${cavevcfsplit}
        """
}

process caveman_flag {
    input:
        path('ref')
        tuple val(groupId), val(types), val(sampleIds), val(protocols), val(platforms), file(htsfiles), file(htsindexes), file(htsStats), path('pindel.germline.bed.gz'), path('pindel.germline.bed.gz.tbi'), path(splitvcf)
        path('caveman')
        path('vagrent')

    output:
        tuple val(groupId), path('flagged.vcf'), emit: flagged

    shell = ['/bin/bash', '-euo', 'pipefail']

    stub:
        """
        touch flagged.vcf
        """

    script:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        """
        # Get the species and assembly from the dict file
        SPECIES=`head -n 2 ref/genome.fa.dict | tail -n 1 | perl -ne 'm/SP:([^\t]+)/;print \$1;'`
        ASSEMBLY=`head -n 2 ref/genome.fa.dict | tail -n 1 | perl -ne 'm/AS:([^\t]+)/;print \$1;'`

        cgpFlagCaVEMan.pl \
            -i $splitvcf \
            -o flagged.vcf \
            -s "\$SPECIES" \
            -sa "\$ASSEMBLY" \
            -m ${htsfiles[case_idx]} \
            -n ${htsfiles[ctrl_idx]} \
            -ref ref/genome.fa.fai \
            -t ${protocols[case_idx]} \
            -b ./caveman/flagging \
            -g pindel.germline.bed.gz \
            -umv ./caveman \
            -ab ./vagrent \
            -v ./caveman/flagging/flag.to.vcf.convert.ini \
            -c ./caveman/flag.vcf.config.${protocols[case_idx]}.ini
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

    stub:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        """
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.snvs.flagged.vcf.gz
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.snvs.flagged.vcf.gz.tbi
        """

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

process vagrent {
    input:
        tuple val(groupId), path(in_vcf), path(in_tbi), val(types), val(sampleIds)
        path('vagrent')

    output:
        path('*.annotated.vcf.gz*')

    publishDir {
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        "${params.outdir}/${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}/annotated"
    }, mode: 'copy'

    shell = ['/bin/bash', '-euo', 'pipefail']

    stub:
        def annot = in_vcf.toString().minus('vcf.gz') + 'annotated.vcf'
        """
        # include full input path name to prevent clash
        touch ${annot}.gz
        touch ${annot}.gz.tbi
        """

    script:
        def annot = in_vcf.toString().minus('vcf.gz') + 'annotated.vcf'
        """
        AnnotateVcf.pl -t -i ${in_vcf} -o ${annot} -c vagrent/vagrent.cache.gz
        """
}

process brass {
    input:
        path('ref')
        tuple val(groupId), val(types), val(sampleIds), val(protocols), val(platforms), file(htsfiles), file(htsindexes), file(htsStats), file('ascat.samplestatistics.txt')
        path('vagrent')
        path('brass')
        path('ref_cache')

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

    stub:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        """
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.brm.bam
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.brm.bam.bai
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.vcf.gz
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.vcf.gz.tbi
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.bedpe.gz
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.bedpe.gz.tbi
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.ngscn.abs_cn.bw
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.intermediates.tar.gz
        """

    script:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        """
        export REF_CACHE=\$PWD/ref_cache/%2s/%2s/%s
        export REF_PATH=\$REF_CACHE

        # Get the species and assembly from the dict file
        SPECIES=`head -n 2 ref/genome.fa.dict | tail -n 1 | perl -ne 'm/SP:([^\t]+)/;print \$1;'`
        ASSEMBLY=`head -n 2 ref/genome.fa.dict | tail -n 1 | perl -ne 'm/AS:([^\t]+)/;print \$1;'`

        brass.pl -j 4 -k 4 -pl ILLUMINA \
            -c ${task.cpus} \
            -d brass/HiDepth.bed.gz \
            -f brass/brass_np.groups.gz \
            -g ref/genome.fa \
            -s "\$SPECIES" -as "\$ASSEMBLY" \
            -pr ${protocols[case_idx]} \
            -g_cache vagrent/vagrent.cache.gz \
            -vi brass/viral.genomic.fa.2bit \
            -mi brass/all_ncbi_bacteria \
            -b brass/500bp_windows.gc.bed.gz \
            -ct brass/CentTelo.tsv \
            -cb brass/cytoband.txt \
            -t ${htsfiles[case_idx]} \
            -n ${htsfiles[ctrl_idx]} \
            -ss ascat.samplestatistics.txt \
            -o ./
        """
}

process verifybamid {
    input:
    //idx, types, samp, htsread, htsidx, cn
        tuple val(groupId), val(types), val(sampleIds), file(htsfiles), file(htsindexes), path('copynumber.caveman.csv')
        path('verifyBamID_snps.vcf.gz')
        path('ref_cache')

    output:
        path('case')
        path('ctrl')
        path('*.contamination.result.json')

    publishDir {
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        "${params.outdir}/${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}/verifyBamID"
    }, mode: 'copy'

    shell = ['/bin/bash', '-euo', 'pipefail']

    stub:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        """
        mkdir -p case
        mkdir -p ctrl
        touch ${sampleIds[case_idx]}_vs_${sampleIds[ctrl_idx]}.contamination.result.json
        """

    script:
        def case_idx = types.indexOf('case')
        def ctrl_idx = types.indexOf('control')
        """
        export REF_CACHE=\$PWD/ref_cache/%2s/%2s/%s
        export REF_PATH=\$REF_CACHE

        # control
        verifyBamHomChk.pl -d 25 \
        -o ./ctrl \
        -b ${htsfiles[ctrl_idx]} \
        -t ${task.cpus} \
        -j ${sampleIds[ctrl_idx]}.contamination.result.json \
        -s verifyBamID_snps.vcf.gz

        # case
        verifyBamHomChk.pl -d 25 \
        -o ./case \
        -b ${htsfiles[case_idx]} \
        -t ${task.cpus} \
        -j ${sampleIds[case_idx]}.contamination.result.json \
        -s verifyBamID_snps.vcf.gz \
        -a copynumber.caveman.csv
        """

}

workflow {
    core_ref     = file(params.core_ref)
    snv_indel    = file(params.snv_indel)
    cnv_sv       = file(params.cnv_sv)
    annot        = file(params.annot)
    qc_genotype  = file(params.qc_genotype)

    pairs = Channel.fromPath(params.pairs)
    exclude_file = file(params.exfile)

    case_control_map = pairs.splitCsv(header: true).map { row -> tuple(row.groupId, row.type, row.sampleId, row.protocol, row.platform, file(row.reads), file(row.readIdx), file(row.readStats)) }

    main:
        obtain_pipeline_metadata(
            ch_repository,
            ch_commitId,
            ch_revision,
            ch_scriptName,
            ch_scriptFile,
            ch_projectDir,
            ch_launchDir,
            ch_workDir,
            ch_userName,
            ch_commandLine,
            ch_configFiles,
            ch_profile,
            ch_container,
            ch_containerEngine,
            ch_raci_owner,
            ch_domain_keywords
        )
        prep_ref(
            core_ref,
            snv_indel,
            cnv_sv,
            annot,
            qc_genotype
        )

        grouped_align = case_control_map.groupTuple()

        genotypes(
            grouped_align.map { idx, type, samp, prot, plat, reads, ridx, bas -> [idx, type, samp, reads, ridx] },
            prep_ref.out.snps_genotype,
            prep_ref.out.snps_sex,
            prep_ref.out.ref_cache
        )

        ascat_counts(
            prep_ref.out.ref,
            prep_ref.out.snps_gc,
            prep_ref.out.snps_sex,
            case_control_map
        )
        ascat(
            prep_ref.out.ref,
            prep_ref.out.snps_gc,
            ascat_counts.out.to_ascat.groupTuple()
        )

        pindel(
            prep_ref.out.ref,
            prep_ref.out.pindel_hidepth,
            grouped_align,
            params.exclude,
            exclude_file
        )

        type_samp_prot = grouped_align.map { idx, type, samp, prot, plat, reads, ridx, bas -> [idx, type, samp, prot] }

        pindel_flag(
            prep_ref.out.pindel_flag,
            prep_ref.out.vagrent,
            pindel.out.vcf.combine(type_samp_prot, by: 0),
            params.skipgerm
        )

        align_and_ascat = grouped_align.combine(ascat.out.ascat_for_caveman, by: 0)

        caveman(
            prep_ref.out.ref,
            align_and_ascat,
            prep_ref.out.cave_hidepth,
            params.cavereads,
            params.exclude,
            params.cavevcfsplit
        )

        // handle cases where split file list is a single element and currently unable to force to a list.
        cleaned_split = caveman.out.split_vcf.map {
            idx, maybe_list -> [idx, maybe_list instanceof Collection ? maybe_list : [maybe_list]]
        }
        flag_set = grouped_align.combine(pindel_flag.out.germline, by:0)
                    .combine(cleaned_split, by: 0)
                    .transpose(by: 10)

        caveman_flag(
            prep_ref.out.ref,
            flag_set,
            prep_ref.out.cave_flag,
            prep_ref.out.vagrent
        )

        type_samp = grouped_align.map { idx, type, samp, prot, plat, reads, ridx, bas -> [idx, type, samp] }

        caveman_flag_merge(
            caveman_flag.out.flagged.groupTuple(by: 0).combine(type_samp, by: 0)
        )

        to_annotate = pindel_flag.out.flagged.concat(caveman_flag_merge.out.flagged)

        vagrent(
            to_annotate.combine(type_samp, by: 0),
            prep_ref.out.vagrent
        )

        align_cn_verify = align_and_ascat.map {
            idx, types, samp, prot, plat, htsread, htsidx, htsstats, cn, ss -> [idx, types, samp, htsread, htsidx, cn]
        }

        verifybamid(
            align_cn_verify,
            prep_ref.out.snps_verifyBamID,
            prep_ref.out.ref_cache
        )

        // just the bits we need staging for this process
        sample_stats = ascat.out.ascat_for_caveman.map { idx, counts, samp_stats -> [idx, samp_stats] }
        align_and_ss = grouped_align.combine(sample_stats, by: 0)

        brass(
            prep_ref.out.ref,
            align_and_ss,
            prep_ref.out.vagrent,
            prep_ref.out.brass,
            prep_ref.out.ref_cache
        )
}
