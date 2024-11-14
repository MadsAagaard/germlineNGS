#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


date=new Date().format( 'yyMMdd' )
user="$USER"
runID="${date}.${user}"
server=

//Unset parameters
params.help                     =false
params.panel                    =null
params.samplesheet              =null
params.preprocessOnly           =null
params.keepwork                 =null
params.nomail                   =null
params.hg38v1                   =null
params.hg38v2                   =null
params.cram                     =null
params.fastq                    =null
params.spring                   =null

params.skipJointGenotyping      =null
params.fastqInput               =null

params.skipSV                   =null
params.skipVariants             =null
params.skipQC                   =null
params.skipSTR                  =null
params.skipSMN                  =null
params.subdirs                  =null
params.gatk                     =null
params.copyCram                 =null
params.single                   =null
//Preset parameters:

params.server                   = "lnx01"
params.genome                   = "hg38"
params.outdir                   = "${launchDir.baseName}.Results"
params.rundir                   = "${launchDir.baseName}"


def helpMessage() {
    log.info"""

    Usage:

    KG Vejle Germline script (WGS, WES or panels)

    PANEL ANALYSIS:
    See https://github.com/KGVejle/germlineNGS for a description of the most common usecases

    WGS ANALYSIS:

    Example samplesheet for standard trio:
    johnDoe 123456789012    index   affected
    johnDoe 234567890123    mater   normal
    johnDoe 345678901234    pater   normal

    The above information can usually be extracted directly from the sample overview excel file

    If the inputdata (FastQ or CRAM) have been transferred to the data archive (which it is by default), the script will automatically find the relevant inputdata  and create symlinks for them in the output (results) directory.

    The script will automatically look for FastQ or CRAM files in subfolders at /lnx01_data2/shared/dataArchive/. This location contains read-only access to the data archive, containing all FastQ and CRAM files. There's no need to copy or move any data.

    The user can point to a specific folder containing raw data (FastQ) using the --fastq option  or alignment data (CRAM) using the --cram option
    This is only needed if input data (FastQ or CRAM) exists outside the data archive (e.g. if data are in personal folders), or if the script is run without samplesheet.

    If the script is run without samplesheet, the user MUST point to a folder containing inputdata with either the --fastq or --cram option.

    Main options:
      --help            Print this help message
      
      --genome          hg19 or hg38
                            Default: hg38 v3 (masked + decoys)

      --hg38v1          Use primary (full) hg38 assembly (UCSC primary).

      --hg38v2          Use hg38 v2 (ucsc.hg38.NGS.analysisSet.fa).

      --gatk            "danak" (v.4.1.9) or "new" (v.4.4.0.0)
                            Default: danak  
      
      --samplesheet     Path to samplesheet for samples to be analyzed (Only required for WGS analysis)
      
      --fastq            Path to folder with wgs fastq files
                            Default: /lnx01_data2/shared/dataArchive/{all subfolders}

      --fastqInput      Use fastq as input (i.e. perform trimming, and alignment)
                            Default: Not set - use CRAM as input.
      
      --cram             Path to folder with wgs CRAM files
                            Default: /lnx01_data2/shared/dataArchive/{all subfolders}
      
      --outdir          Manually set output directory
                            Default: {current_dir}.Results

      --keepwork        keep the workfolder generated by the nextflow script.
                            Default: Not set - removes the Work folder

      --nomail          Does not send a mail-message when completing a script
                            Default: Not set - sends mail message if the user is mmaj or raspau and only if the script has been running longer than 20 minutes.

      --panel           Type of paneldata to analyze. Currently supports AV1, CV5 and MV1
                            Default: Not set - assumes WGS data by dfault


    WGS Analysis: Select or modify analysis steps:

      --skipVariants    Do not call SNPs and INDELs at all
                            Default: Call SNPs and INDELs using GATK HaplotypeCaller

      --skipSV          Do not call Structural Variants (SV incl. CNVs) at all
                            Default: Call SVs using Manta, Lumpy, CNVNator and CNVKit

      --skipSTR         Do not call repeat expansions.
                            Default: Calls repeat expansions using Stripy and ExpansionHunter

      --skipQC          Do not run QC module (e.g. Picard Metrics, samtools, multiQC etc.)
                            Default: Run QC module

      --skipSMN         Do not call SMN1 and SMN2 variants
                            Default: Call SMN variants with SMNCopyNumberCaller

    """.stripIndent()
}
if (params.help) exit 0, helpMessage()

def errorMessage1() {

    log.info"""

    USER INPUT ERROR: If no samplesheet is selected, the user needs to point to a folder containing relevant fastq, CRAM or SPRING files... 
    Run the script with the --help parameter to see available options
    
    """.stripIndent()
}

if (!params.samplesheet && !params.fastq && !params.cram && !params.spring) exit 0, errorMessage1()

def FastqCRAM_error() {
    log.info"""
    USER INPUT ERROR: The user should point to either FastQ (--fastq parameter) or CRAM (--cram parameter) as input - not both! 
    """.stripIndent()
}

if (params.cram && params.fastq) exit 0, FastqCRAM_error()


switch (params.gatk) {

    case 'danak':
    gatk_image="gatk419.sif";
    break;
    case 'new':
    gatk_image="gatk4400.sif";
    break;
    case 'v45':
    gatk_image="gatk4500.sif";
    default:
        if (params.panel=="AV1" || params.panel=="GV3" || params.panel=="CV5"){
            gatk_image="gatk419.sif";
        }
        else {
            gatk_image="gatk4400.sif";
        }
    break;
}



switch (params.server) {
    case 'lnx02':
        s_bind="/data/:/data/,/lnx01_data2/:/lnx01_data2/,/fast/:/fast/,/lnx01_data3/:/lnx01_data3/";
        simgpath="/data/shared/programmer/simg";
        tmpDIR="/fast/TMP/TMP.${user}/";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/${gatk_image} gatk";
        multiqc_config="/data/shared/programmer/configfiles/multiqc_config.yaml"
        tank_storage="/home/mmaj/tank.kga/data/data.storage.archive/";

        dataStorage="/lnx01_data3/storage/";
        params.intervals_list="/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/hg38v3/hg38v3_scatter20_BWI/*.interval_list";
        dataArchive="/lnx01_data2/shared/dataArchive";
        refFilesDir="/fast/shared/genomes";
    break;

    case 'lnx01':
        s_bind="/data/:/data/,/lnx01_data2/:/lnx01_data2/";
        simgpath="/data/shared/programmer/simg";
        tmpDIR="/data/TMP/TMP.${user}/";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/${gatk_image} gatk";
        multiqc_config="/data/shared/programmer/configfiles/multiqc_config.yaml"
        tank_storage="/home/mmaj/tank.kga2/data/data.storage.archive/";
        dataStorage="/lnx01_data3/storage/";
        //params.intervals_list="/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/hg38v3/hg38v3_scatter10_IntervalSubdiv/*.interval_list";
        params.intervals_list="/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/hg38v3/hg38v3_scatter20_BWI/*.interval_list";
        modules_dir="/home/mmaj/scripts_lnx01/nextflow_lnx01/dsl2/modules";
        subworkflow_dir="/home/mmaj/scripts_lnx01/nextflow_lnx01/dsl2/subworkflows";
        dataArchive="/lnx01_data2/shared/dataArchive";
        refFilesDir="/data/shared/genomes";
    break;

    case 'rgi01':
        s_bind="/data/:/data/,/lnx01_data2/:/lnx01_data2/,/fast/:/fast/,/lnx01_data3/:/lnx01_data3/";
        simgpath="/data/shared/programmer/simg";
        tmpDIR="/fast/TMP/TMP.${user}/";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/${gatk_image} gatk";
        multiqc_config="/data/shared/programmer/configfiles/multiqc_config.yaml"
        tank_storage="/home/mmaj/tank.kga/data/data.storage.archive/";

        dataStorage="/lnx01_data3/storage/";
        params.intervals_list="/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/hg38v3/hg38v3_scatter20_BWI/*.interval_list";
        dataArchive="/lnx01_data2/shared/dataArchive";
        refFilesDir="/fast/shared/genomes";
    break;

    case 'kga01':
        modules_dir="/home/mmaj/LNX01_mmaj/scripts_lnx01/nextflow_lnx01/dsl2/modules";
        subworkflow_dir="/home/mmaj/LNX01_mmaj/scripts_lnx01/nextflow_lnx01/dsl2/subworkflows";
        dataArchive="/data/shared/dataArchive";

    break;
}



switch (params.panel) {

    case "AV1":
        reads_pattern_cram="*{.,-,_}{AV1}{.,-,_}*.cram";
        reads_pattern_crai="*{.,-,_}{AV1}{.,-,_}*.crai";
        reads_pattern_fastq="*{.,-,_}{AV1}{.,-,_}*R{1,2}*{fq,fastq}.gz";
        reads_pattern_spring="*AV1*.spring";
        panelID="AV1"
    break;

    case "CV5":
        reads_pattern_cram="*{.,-,_}{CV5}{.,-,_}*.cram";
        reads_pattern_crai="*{.,-,_}{CV5}{.,-,_}*.crai";
        reads_pattern_fastq="*{.,-,_}{CV5}{.,-,_}*R{1,2}*{fq,fastq}.gz";
        reads_pattern_spring="*CV5*.spring";
        panelID="CV5"
    break;

    case "GV3":
        reads_pattern_cram="*{GV1,GV2,GV3}*.cram";
        reads_pattern_crai="*{GV1,GV2,GV3}*.crai";
        reads_pattern_fastq="*{GV1,GV2,GV3}*R{1,2}*{fq,fastq}.gz";
        reads_pattern_spring="*{GV1,GV2,GV3}*.spring";
        panelID="GV3"
    break;

    case "GV_TEST":
        reads_pattern_cram="*.cram";
        reads_pattern_crai="*.crai";
        reads_pattern_fastq="*R{1,2}*{fq,fastq}.gz";
        panelID="GV_TEST"
    break;

    case "MV1":
        reads_pattern_cram="*{.,-,_}{MV1}{.,-,_}*.cram";
        reads_pattern_crai="*{.,-,_}{MV1}{.,-,_}*.crai";
        reads_pattern_fastq="*{.,-,_}{MV1}{.,-,_}*R{1,2}*{fq,fastq}.gz";
        reads_pattern_spring="*MV1*.spring";
        panelID="MV1"
    break;

    case "WES_2":
        reads_pattern_cram="*{-,.,_}{EV8,EV7,EV6}{-,.,_}*.cram";
        reads_pattern_crai="*{-,.,_}{EV8,EV7,EV6}{-,.,_}*.crai";
        reads_pattern_fastq="*{-,.,_}{EV8,EV7,EV6}{-,.,_}*R{1,2}*{fq,fastq}.gz";
        reads_pattern_spring="*{EV8,EV7,EV6}*.spring";
        panelID="WES"
    break;

    case "WES":
        reads_pattern_cram="*{-,.,_}{EV8_ALM,EV8_ONK}{-,.,_}*.cram";
        reads_pattern_crai="*{-,.,_}{EV8_ALM,EV8_ONK}{-,.,_}*.crai";
        reads_pattern_fastq="*{-,.,_}{EV8_ALM,EV8_ONK}{-,.,_}*R{1,2}*{fq,fastq}.gz";
        reads_pattern_spring="*{EV8_ALM,EV8_ONK}*.spring";
        panelID="WES_subpanel"
    break;

    case "WGS_CNV":
        reads_pattern_cram="*{-,.,_}{WG4_CNV}{-,.,_}*.cram";
        reads_pattern_crai="*{-,.,_}{WG4_CNV}{-,.,_}*.crai";
        reads_pattern_fastq="*{-,.,_}{WG4_CNV}{-,.,_}*R{1,2}*{fq,fastq}.gz";
        reads_pattern_spring="*{WG4_CNV}*.spring";
        panelID="WGS"
    break;

    case "NGC":
        reads_pattern_cram="*{-,.,_}{WG4_NGC,NGCWGS}{-,.,_}*.cram";
        reads_pattern_crai="*{-,.,_}{WG4_NGC,NGCWGS}{-,.,_}*.crai";
        reads_pattern_fastq="*{-,.,_}{WG4_NGC,NGCWGS}{-,.,_}*R{1,2}*{fq,fastq}.gz";
        reads_pattern_spring="*{WG4_NGC,NGCWGS}*.spring";
        panelID="WGS"
    break;

    default: 
        reads_pattern_cram="*.cram";
        reads_pattern_crai="*.crai";
        reads_pattern_fastq="*R{1,2}*{fq,fastq}.gz";
        reads_pattern_spring="*.spring";
        panelID="ALL"
    break;
}


switch (params.genome) {
    case 'hg19':
        assembly="hg19"
        // Genome assembly files:
        genome_fasta = "/data/shared/genomes/hg19/human_g1k_v37.fasta"
        genome_fasta_fai = "/data/shared/genomes/hg19/human_g1k_v37.fasta.fai"
        genome_fasta_dict = "/data/shared/genomes/hg19/human_g1k_v37.dict"
        genome_version="V1"
        break;


    case 'hg38':
        assembly="hg38"
        spliceai_assembly="grch38"
        smncaller_assembly="38"
        svdb_databases="/data/shared/genomes/hg38/inhouse_DBs/hg38v3/svdb_AF"
        // Genome assembly files:
        if (params.hg38v1) {
        genome_fasta = "${refFilesDir}/hg38/GRCh38.primary.fa"
        genome_fasta_fai = "${refFilesDir}/hg38/GRCh38.primary.fa.fai"
        genome_fasta_dict = "${refFilesDir}/hg38/GRCh38.primary.dict"
        genome_version="hg38v1"
        cnvkit_germline_reference_PON="/data/shared/genomes/hg38/inhouse_DBs/hg38v1_primary/cnvkit/wgs_germline_PON/jgmr_45samples.reference.cnn"
        cnvkit_inhouse_cnn_dir="/data/shared/genomes/hg38/inhouse_DBs/hg38v1_primary/cnvkit/wgs_persample_cnn/"
        inhouse_SV="/data/shared/genomes/hg38/inhouse_DBs/hg38v1_primary/"
        }
        
        if (params.hg38v2){
        genome_fasta = "${refFilesDir}/hg38/ucsc.hg38.NGS.analysisSet.fa"
        genome_fasta_fai = "${refFilesDir}/hg38/ucsc.hg38.NGS.analysisSet.fa.fai"
        genome_fasta_dict = "${refFilesDir}/hg38/ucsc.hg38.NGS.analysisSet.dict"
        genome_version="hg38v2"
        }

        // Current hg38 version (v3): NGC with masks and decoys.
        if (!params.hg38v2 && !params.hg38v1){
        genome_fasta = "${refFilesDir}/hg38/GRCh38_masked_v2_decoy_exclude.fa"
        genome_fasta_fai = "${refFilesDir}/hg38/GRCh38_masked_v2_decoy_exclude.fa.fai"
        genome_fasta_dict = "${refFilesDir}/hg38/GRCh38_masked_v2_decoy_exclude.dict"
        genome_version="hg38v3"
        cnvkit_germline_reference_PON="/data/shared/genomes/hg38/inhouse_DBs/hg38v3_primary/cnvkit/hg38v3_109samples.cnvkit.reference.cnn"
        cnvkit_inhouse_cnn_dir="/data/shared/genomes/hg38/inhouse_DBs/hg38v3_primary/cnvkit/wgs_persample_cnn/"
        inhouse_SV="/data/shared/genomes/hg38/inhouse_DBs/hg38v3/"
        }

        // Gene and transcript annotation files:

        gencode_gtf = "${refFilesDir}/hg38/gene.annotations/gencode.v36.annotation.gtf"
        gencode_gff3 = "${refFilesDir}/hg38/gene.annotations/gencode.v36.annotation.gff3"
     
        //Program  files:
        msisensor_list="${refFilesDir}/hg38/program_DBs/msisensor/hg38_msisensor_scan.txt"
        
      
        //Structural variants
        delly_exclude="/data/shared/genomes/hg38/program_DBs/delly/human.hg38.excl.tsv"
        
        smoove_exclude="/data/shared/genomes/hg38/interval.files/smoove/smoove.hg38.excluderegions.bed"
        smoove_gff="/data/shared/genomes/hg38/gene.annotations/GRCh38_latest_genomic.gff.gz"


        //inhouse SV AF databases: 
        mantaSVDB="${svdb_databases}/mantaSVDB315.db"
        lumpySVDB="${svdb_databases}/lumpySVDB218.db"
        cnvkitSVDB="${svdb_databases}/cnvkitSVDB313.db"
        //tidditSVDB="${svdb_databases}/tidditSVDB.db"
        dellySVDB="${svdb_databases}/dellySVDB112.db"


        //Repeat Expansions:
        expansionhunter_catalog="/data/shared/genomes/hg38/program_DBs/expansionHunter/expansionHunter_hg38_stripy.variant_catalog.json"
        hipSTR_bed="/data/shared/genomes/hg38/interval.files/STRs/GRCh38.hipstr_reference.bed"

        // Somatic calling files (GATK Mutect2 pipeline):
        gatk_wgs_pon="/data/shared/genomes/hg38/program_DBs/GATK/somatic/somatic-hg38_1000g_pon.hg38.vcf.gz"
        mutect_gnomad="/data/shared/genomes/hg38/program_DBs/GATK/somatic/somatic-hg38_af-only-gnomad.hg38.vcf.gz"
        gatk_contamination_ref="/data/shared/genomes/hg38/program_DBs/GATK/somatic/somatic-hg38_small_exac_common_3.hg38.vcf.gz"

        // Program indexes:
        pcgr_assembly="grch38"
        sequenza_cg50_wig="/data/shared/genomes/hg38/program_DBs/sequenza/GRCh38.primary.cg50.sequenza.wig.gz"


        // Regions & variants:
        qualimap_ROI="/data/shared/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.20bp.SM.6col.bed"
        gencode_exons_ROI="/data/shared/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.SM.bed"

        ROI="/data/shared/genomes/hg38/interval.files/exome.ROIs/211130.hg38.refseq.gencode.fullexons.50bp.SM.bed"
        
        //ROI="/data/shared/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.20bp.SM.bed"

        callable_regions="/data/shared/genomes/hg38/interval.files/GATK.hg38.callable.regions.bed"
        manta_callable_regions="/data/shared/genomes/hg38/interval.files/manta/GATK.hg38.callable.regions.bed.gz"

        dbsnp="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
        KGindels="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz"
        KGindels_idx="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"

        KGmills="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
        KGmills_idx="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
        KG_p1_High_snps="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz"

        hapmap="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz"
        omni="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz"
        AV1_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/av1.hg38.ROI.v2.bed"
        CV1_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/cv3.hg38.ROI.bed"
        CV2_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/cv3.hg38.ROI.bed"
        CV3_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/cv3.hg38.ROI.bed"
        CV4_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/cv4.hg38.ROI.bed"
        CV5_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/cv5.hg38.ROI.bed"
        GV3_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/gv3.hg38.ROI.v2.bed"
        NV1_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/nv1.hg38.ROI.bed"
        WES_ROI="/data/shared/genomes/hg38/interval.files/exome.ROIs/211130.hg38.refseq.gencode.fullexons.50bp.SM.bed"
        MV1_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/muc1.hg38.coordinates.bed"
        break;
}

channel
    .fromPath(params.intervals_list)
    .map { it -> tuple(it.baseName,it)}
    .set { haplotypecallerIntervalList }
////////////////////////////////////////////////////
////// INPUT DATA: FASTQ AS INPUT //////////////////
////////////////////////////////////////////////////

if (!params.fastq && params.fastqInput) {

    inputFastq="${dataArchive}/{lnx01,lnx02}/**/${reads_pattern_fastq}"
}
if (params.fastq) {
    inputFastq="${params.fastq}/${reads_pattern_fastq}"
}

if (params.fastq || params.fastqInput) {
    Channel.fromFilePairs("${inputFastq}", checkIfExists: true)
    | map { id, reads -> 
        (sample, ngstype)   = reads[0].baseName.tokenize("-")
        (panel,subpanel)    = ngstype.tokenize("_")
        meta = [id:sample+"_"+ngstype, npn:sample, fullpanel:ngstype,panel:panel, subpanel:subpanel]
        [meta, reads]
        |view
    }

    | branch {meta, reads ->
        WGS: (meta.panel=~/WG/ || meta.panel=~/NGC/)
            return [meta + [datatype:"WGS",roi:"$WES_ROI"],reads]
        AV1: (meta.panel=~/AV1/)
            return [meta + [datatype:"targeted",roi:"$AV1_ROI"],reads]
        MV1: (meta.panel=~/MV1/)
            return [meta + [datatype:"targeted",roi:"$MV1_ROI"],reads]
        WES: (meta.panel=~/EV8/ ||meta.panel=~/EV7/)
            return [meta + [datatype:"targeted",roi:"$WES_ROI"],reads]
        undetermined: true
            return [meta + [datatype:"unset",analyzed:"NO"],reads]
        [meta, reads]
    }
    | set {readsInputBranched}

    readsInputBranched.MV1.concat(readsInputBranched.AV1).concat(readsInputBranched.WES).concat(readsInputBranched.WGS)
    | set {readsInputReMerged}

readsInputReMerged.view()
    if (params.samplesheet) {
        readsInputReMerged
        | map { meta,reads -> tuple(meta.npn,meta,reads)}
        | set {readsInputForJoin}
  //      readsInputForJoin.view()
        // NBNBNBNBNB: Requires named headers for now!!! (e.g. column with NPN must be named "npn" in samplesheet)
    channel.fromPath(params.samplesheet)
        | splitCsv(sep:'\t',header:true)
        | map { row -> tuple(row.npn, row)}
       // | view
        | set { full_samplesheet }
    full_samplesheet.join(readsInputForJoin)    
        | map {tuple(it[1],it[2],it[3])}
        | map {meta1,meta2,data -> 
          [meta1+meta2,data]}
        | set {readsInputFinal}
    }
    
    if (!params.samplesheet) {
        readsInputReMerged
        | set {readsInputFinal} 
    }
   // readsInputFinal.view()
}
////////////////////////////////////////////////////
////// INPUT DATA: CRAM AS INPUT //////////////////
////////////////////////////////////////////////////

if (!params.fastq && !params.fastqInput && !params.spring){

    if (params.cram) {
        cramfiles="${params.cram}/${reads_pattern_cram}"
        craifiles="${params.cram}/${reads_pattern_crai}"
    }

    if (!params.cram) {
        cramfiles="${dataArchive}/{lnx01,lnx02,tank_kga_external_archive}/**/${reads_pattern_cram}"
        craifiles="${dataArchive}/{lnx01,lnx02,tank_kga_external_archive}/**/${reads_pattern_crai}"
    }
    
    if (params.cram && params.subdirs) {
        cramfiles="${params.cram}/**/${reads_pattern_cram}"
        craifiles="${params.cram}/**/${reads_pattern_crai}"
    }

    Channel.fromPath(cramfiles,checkIfExists:true)
    |map {tuple (it.simpleName,it)}
    |set {cramfiles}
    
    Channel.fromPath(craifiles,checkIfExists:true)
    |map {tuple (it.simpleName,it)}
    |set {craifiles}
    
    cramfiles.join(craifiles)
    |map { id, aln, index ->
        (sample, panel,subpanel)   = id.tokenize("_")
      //  (panel,subpanel)    = ngstype.tokenize("_")
        meta = [id:id, npn:sample, fullpanel:panel+"_"+subpanel, panel:panel, subpanel:subpanel]
        tuple(meta,[aln,index])
    }

    | set {cram_all}
    cram_all
    |branch {meta, aln ->
            WGS: (meta.panel=~/WG/ || meta.panel=~/NGC/)
                return [meta + [datatype:"WGS",roi:"$WES_ROI"],aln]
            AV1: (meta.panel=~/AV1/)
                return [meta + [datatype:"targeted",roi:"$AV1_ROI"],aln]
            MV1: (meta.panel=~/MV1/)
                return [meta + [datatype:"targeted",roi:"$MV1_ROI"],aln]
            WES: (meta.panel=~/EV8/ ||meta.panel=~/EV7/)
                return [meta + [datatype:"targeted",roi:"$WES_ROI"],aln]
            undetermined: true
                return [meta + [datatype:"unset",analyzed:"NO"],aln]
            [meta, aln]
    }
    | set {cramInputBranched}

    cramInputBranched.MV1.concat(cramInputBranched.AV1).concat(cramInputBranched.WES).concat(cramInputBranched.WGS)
    |set {cramInputReMerged}

    if (params.samplesheet) {
        cramInputReMerged
        | map { meta,aln -> tuple(meta.npn,meta,aln)}
        | set {alnInputForJoin}
        // NBNBNBNBNB: Requires named headers for now!!! (e.g. column with NPN must be named "npn" in samplesheet)
         channel.fromPath(params.samplesheet)
        | splitCsv(sep:'\t',header:true)
        | map { row -> tuple(row.npn, row)}
        | set { full_samplesheet }

    full_samplesheet.join(alnInputForJoin)    
        | map {tuple(it[1],it[2],it[3])}
        | map {meta1,meta2,data -> 
          [meta1+meta2,data]}
        | set {alnInputFinal}
    }
    
    if (!params.samplesheet) {
        cramInputReMerged
        | set {alnInputFinal} 
    }
    //alnInputFinal.view()
}


////////////////////////////////////////////////////////////////////
//// NEW June 2024: Add spring as input. ///////////////////////////
////////////////////////////////////////////////////////////////////


if (params.spring && !params.samplesheet) {

    params.spring_reads="${params.spring}/${reads_pattern_spring}"


    Channel
    .fromPath(params.spring_reads, checkIfExists: true)
    .map { tuple(it.baseName.tokenize('-').get(0)+"_"+it.baseName.tokenize('-').get(1),it) }
    .set {spring_input_ch}
}







///// Haplotypecaller splitintervals channel: /////
/*
channel
    .fromPath(params.intervals_list)
    .map { it -> tuple(it.baseName,it)}
    .set { haplotypecallerIntervalList }
*/
////////////////////////////////////////////////////

include { 
         // Symlinks:
         inputFiles_symlinks_cram;
         inputFiles_cramCopy;
         // Preprocess tools:
         //QC tools
         samtools;
         bamtools;
         qualimap;
         fastqc_bam;
         collectWGSmetrics;
         multiQC;
         vntyper_newRef;
         haplotypecallerSplitIntervals;
         //subworkflows:
         SUB_SPRING_DECOMPRESS;
         SUB_PREPROCESS;
         SUB_QC;
         SUB_VARIANTCALL;
         SUB_VARIANTCALL_WGS;
         SUB_CNV_SV;
         SUB_STR;
         SUB_SMN } from "./modules/modules.dna.v1.nf" 

/*
workflow QC {
    take: 
    meta_aln_index
    main:
    samtools(meta_aln_index)

    multiQC(samtools.out.ifEmpty([]).mix(qualimap.out.ifEmpty([])).mix(fastqc_bam.out.ifEmpty([])).collect())

}


process fastq_to_ubam {
    errorStrategy 'ignore'
    tag "$meta.id"
    //publishDir "${outputDir}/unmappedBAM/", mode: 'copy',pattern: '*.{bam,bai}'
    //publishDir "${outputDir}/fastq_symlinks/", mode: 'link', pattern:'*.{fastq,fq}.gz'
    cpus 20
    maxForks 10

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.unmapped.from.fq.bam"), path("${meta.id}.unmapped.from.fq.bam.idx"),emit: testOut
    
    script:
    """
    ${gatk_exec} FastqToSam \
    -F1 ${reads[0]} \
    -F2 ${reads[1]} \
    -SM ${meta.id} \
    -PL illumina \
    -PU KGA_PU \
    -RG KGA_RG \
    --TMP_DIR ${tmpDIR} \
    -O ${meta.id}.unmapped.from.fq.bam
    touch ${meta.id}.unmapped.from.fq.bam.idx
    """
}

*/


workflow {

    if (params.spring) {
        SUB_SPRING_DECOMPRESS(spring_input_ch)
        //SUB_SPRING_DECOMPRESS.out.view()
        readsInputFinal=SUB_SPRING_DECOMPRESS.out.fq_read_input_spring
    }

    if (params.fastqInput||params.fastq||params.spring || params.preprocessOnly) {
        SUB_PREPROCESS(readsInputFinal)
   //     preprocessOutAln=SUB_PREPROCESS.out.finalAln
        
     //   preprocessOutAln
        SUB_PREPROCESS.out.finalAln
        |map {meta, cram,crai ->
            tuple(meta,[cram,crai])}
        |set {alnInputFinal}
    }
    alnInputFinal
        |branch {meta, aln ->
            WGS: (meta.datatype=~/WGS/)
                return [meta ,aln]
            TARGETED: (meta.datatype=~/targeted/)
                return [meta ,aln]
        }
        | set {alnInputFinalBranched}

    if (!params.preprocessOnly) {

        if(!params.skipQC){
            SUB_QC(alnInputFinalBranched.WGS)
        }

        if (!params.skipVariants) {
            SUB_VARIANTCALL(alnInputFinalBranched.TARGETED)
            SUB_VARIANTCALL_WGS(alnInputFinalBranched.WGS)
        }

        if (!params.skipSV) {
            SUB_CNV_SV(alnInputFinalBranched.WGS)
        }
        
        if (!params.skipSTR) {
            SUB_STR(alnInputFinalBranched.WGS)
        }

        if (!params.skipSMN) {
            SUB_SMN(alnInputFinalBranched.WGS)
        }
    }

}
    /*


    haplotypecallerSplitIntervals(alnInputFinalBranched.WGS.combine(haplotypecallerIntervalList))
    haplotypecallerSplitIntervals.out
    //|view
    |set {splitintervalout_test}
    splitintervalout_test
    .groupTuple(size:17)
    | view    



        |branch {meta, aln ->
            WGS: (meta.panel=~/WG/ || meta.panel=~/NGC/)
                return [meta + [datatype:"WGS",roi:"$WES_ROI"],aln]
            AV1: (meta.panel=~/AV1/)
                return [meta + [datatype:"targeted",roi:"$AV1_ROI"],aln]
            MV1: (meta.panel=~/MV1/)
                return [meta + [datatype:"targeted",roi:"$MV1_ROI"],aln]
            WES: (meta.panel=~/EV8/ ||meta.panel=~/EV7/)
                return [meta + [datatype:"targeted",roi:"$WES_ROI"],aln]
            undetermined: true
                return [meta + [datatype:"unset",analyzed:"NO"],aln]
            [meta, aln]
        }
              | set {cramInputBranched}
    
        cramInputBranched.MV1.concat(cramInputBranched.AV1).concat(cramInputBranched.WES).concat(cramInputBranched.WGS)
        |set {cramInputReMerged}
    
        cramInputReMerged.view()
      */
            
       // }
       // |set {alnInputFinal}

       // alnInputFinal.view()
        //SUB_VARIANTCALL(alnInputFinal)
/*
    Working, 240925:
        preprocessOutAln
        |map {meta, cram,crai ->
            tuple(meta,[cram,crai])}
        |set {alnInputFinal}

        alnInputFinal.view()
    End

    if (!params.fastqInput && !params.fastq && !params.spring) {
        inputFiles_symlinks_cram(alnInputFinal) // new channel structure: val(meta), path(data)

    }


    if (!params.panel || params.panel =="WGS_CNV"|| params.panel =="NGC") { //i.e. if WGS data

        if (!params.skipVariants) {
            SUB_VARIANTCALL_WGS(alnInputFinal)
        }
        if (!params.skipSV) {
            SUB_CNV_SV(alnInputFinal)
        }
        if (!params.skipSTR) {
            SUB_STR(alnInputFinal)
        }
        
        if (!params.skipSMN) {
        SUB_SMN(alnInputFinal)
        }

    }

    if (params.panel && params.panel!="WGS_CNV"&& params.panel!="NGC") {

        SUB_VARIANTCALL(alnInputFinal)

        if (params.panel=="MV1") {
            vntyper_newRef(readsInputFinal)
        }
    }
  */  

















/*


workflow.onComplete {
    // only send email if --nomail is not specified, the user is mmaj or raspau and duration is longer than 5 minutes / 300000 milliseconds
    if (!params.nomail && workflow.duration > 300000 && workflow.success) {
        if (System.getenv("USER") in ["raspau", "mmaj"]) {
            def sequencingRun = params.cram ? new File(params.cram).getName().take(6) :
                               params.fastq ? new File(params.fastq).getName().take(6) : 'Not provided'

            // Checks if there are OBS samples in the cram folder
            def obsSampleMessage = ""
            if (params.panel == "AV1" && params.cram) {
                def cramDir = new File(params.cram)
                def obsSamples = cramDir.listFiles().findAll { it.name.contains("OBS") }
                if (obsSamples.size() > 0) {
                    obsSampleMessage = "\nTHERE IS AN OBS SAMPLE IN THIS RUN"
                }
            }

            def workDirMessage = params.keepwork ? "WorkDir             : ${workflow.workDir}" : "WorkDir             : Deleted"

            // Correctly set the outputDir
            def outputDir = "${launchDir}/${launchDir.baseName}.Results"

            def body = """\
            Pipeline execution summary
            ---------------------------
            Pipeline completed  : ${params.panel}
            Sequencing run      : ${sequencingRun}${obsSampleMessage}
            Duration            : ${workflow.duration}
            Success             : ${workflow.success}
            ${workDirMessage}
            OutputDir           : ${outputDir}
            Exit status         : ${workflow.exitStatus}
            ${obsSampleMessage}
            """.stripIndent()

            // Send the email using the built-in sendMail function
            sendMail(to: 'Andreas.Braae.Holmgaard@rsyd.dk,Annabeth.Hogh.Petersen@rsyd.dk,Isabella.Almskou@rsyd.dk,Jesper.Graakjaer@rsyd.dk,Lene.Bjornkjaer@rsyd.dk,Martin.Sokol@rsyd.dk,Mads.Jorgensen@rsyd.dk,Rasmus.Hojrup.Pausgaard@rsyd.dk,Signe.Skou.Tofteng@rsyd.dk', subject: 'GermlineNGS pipeline Update', body: body)

            // Check if --keepwork was specified
            if (!params.keepwork) {
                // If --keepwork was not specified, delete the work directory
                println("Deleting work directory: ${workflow.workDir}")
                "rm -rf ${workflow.workDir}".execute()
            }
        }
    }    
}

workflow.onError {
    // Custom message to be sent when the workflow completes
    def sequencingRun = params.cram ? new File(params.cram).getName().take(6) :
                   params.fastq ? new File(params.fastq).getName().take(6) : 'Not provided'

    def body = """\
    Pipeline execution summary
    ---------------------------
    Pipeline completed  : ${params.panel}
    Sequencing run      : ${sequencingRun}
    Duration            : ${workflow.duration}
    Failed              : ${workflow.failed}
    WorkDir             : ${workflow.workDir}
    Exit status         : ${workflow.exitStatus}
    """.stripIndent()

    // Send the email using the built-in sendMail function
    sendMail(to: 'Mads.Jorgensen@rsyd.dk,Rasmus.Hojrup.Pausgaard@rsyd.dk', subject: 'Pipeline Update', body: body)

}

*/