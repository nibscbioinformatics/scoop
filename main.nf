#!/usr/bin/env nextflow
/*
========================================================================================
                         nibscbioinformatics/scoop
========================================================================================
 nibscbioinformatics/scoop Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nibscbioinformatics/scoop
----------------------------------------------------------------------------------------
*/


/*
=============================================================================
### SET PARAMETERS OTHERWISE UNDEFINED with default false value
============================================================================
*/

/* ##############################################
*  INITIALISE DATABASE FILES
*  ##############################################
*/

// Check if database exists in the config file
if (params.tool && !params.databases.containsKey(params.tool)) {
    exit 1, "The requested tool '${params.tool}' is not available and no database is associated with it. Currently the available databases are ${params.databases.keySet().join(", ")}"
}

// if params.nucleotide_db is empty, it will default to chocophlan
// if not empty, but is not chocophlan we assume it's a URL or local file
params.nucleotide_db = 'chocophlan'
//params.nucleotide_db = params.nucleotide_db ? params.nucleotide_db : 'chocophlan'

// if params.protein_db is empty, it will default to uniref90
params.protein_db = 'uniref90_diamond'
//params.protein_db = params.protein_db ? params.protein_db : 'uniref90_diamond'
params.metaphlan_db = 'bowtie'
//params.metaphlan_db = params.metaphlan_db ? params.metaphlan_db : 'bowtie'

params.mpa_index = 'v20_m200'
//params.mpa_index = params.mpa_index ? params.mpa_index : 'v20_m200'


// no ifs based on selected tool, if database for that tool doesn't exist will set param as null
params.chocophlan = params.tool ? params.databases[params.tool].chocophlan ?: null : null
params.uniref50_diamond = params.tool ? params.databases[params.tool].uniref50_diamond ?: null : null
params.uniref90_diamond = params.tool ? params.databases[params.tool].uniref90_diamond ?: null : null
params.uniref50_ec_filtered_diamond = params.tool ? params.databases[params.tool].uniref50_ec_filtered_diamond ?: null : null
params.uniref90_ec_filtered_diamond = params.tool ? params.databases[params.tool].uniref90_ec_filtered_diamond ?: null : null

// metaphlan default db needs to be downloaded and prepared anyway
params.bowtie = params.databases['metaphlan2'].bowtie
params.mpamdd5 = params.databases['metaphlan2'].md5


// search mode will vary according to the selected database
params.search_mode = 'uniref90'
//params.search_mode = params.search_mode ? params.search_mode : 'uniref50'

// metaphlan db is packed in a quite peculiar way so we need to distinguish when it's custom
// or when it's default
params.mpaType = 'default'
//params.mpaType = params.mpaType ? params.mpaType : 'default'


// =======================================================================

def helpMessage() {

    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nibscbioinformatics/scoop --input 'sample_file.tsv' -profile docker

    Mandatory arguments:
      --input [file]                  Input data TSV file

      -profile [str]                  Configuration profile to use. Can use multiple (comma separated)
                                      Available: conda, docker, singularity, test, awsbatch, <institute> and more

      --tool [str]                    Pipeline to use. Can be:
                                      humann2 - for the entire functional characterisation pipeline
                                      metaphlan2 - to run just the phylogenetic analysis

      Optional arguments:

      --protein_db [str]              Protein database to be used for Diamond if humann2 workflow has been selected.
                                      Available:
                                      uniref50_diamond, uniref90_diamond (default),
                                      uniref50_ec_filtered_diamond, uniref90_ec_filtered_diamond.
                                      If custom, a folder location or the URL of a TAR file should be provided.

      --nucleotide_db [str]           Nucleotide database to be used if humann2 workflow has been selected.
                                      Default: chocophlan.
                                      If custom, a folder location or the URL of a TAR file should be provided.

      --metaphlan_db [str]            Bowtie DB to be used for metaphlan2.
                                      If empty, uses default mpa_v20_m200 database.
                                      If custom, a folder location with properly bowtie indexed db should be provided or URL with a TAR file should be provided.
                                      The file name prefix should always be in the format mpa_index where index is specified with following option.

      --mpa_index [str]               If custom metaphlan2 database is provided, the second part of the prefix in the file name needs to be specificed, otherwise
                                      the software will assume it is v20_m200.
                                      Do not include any preceding underscore which will be added by default between the mandatory mpa previx and the index.

    Other options:
      --outdir [file]                 The output directory where the results will be saved
      --email [email]                 Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail [email]         Same as --email, except only send mail if the workflow is not successful
      --max_multiqc_email_size [str]  Theshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch options:
      --awsqueue [str]                The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]               The AWS Region for your AWS Batch job to run on
      --awscli [str]                  Path to the AWS CLI tool
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)



// channels are set based on parameters
// or as local file or URL if none of the name matches (i.e. using switch default)

ch_nucleotide_db = Channel.empty()
ch_protein_db = Channel.empty()
ch_metaphlan_db = Channel.empty()


if (params.tool == 'humann2'){
  switch (params.nucleotide_db) {
    case 'chocophlan':
          ch_nucleotide_db = params.chocophlan ? Channel.value(file(params.chocophlan)) : "null";
          break;
    default:
          ch_nucleotide_db = params.nucletide_db ? Channel.value(file(params.nucletide_db)) : "null";
          break;
  }
  switch (params.protein_db) {
    case 'uniref50_diamond':
          ch_protein_db = params.uniref50_diamond ? Channel.value(file(params.uniref50_diamond)) : "null";
          params.search_mode = 'uniref50';
          break;

    case 'uniref90_diamond':
          ch_protein_db = params.uniref90_diamond ? Channel.value(file(params.uniref90_diamond)) : "null";
          break;

    case 'uniref50_ec_filtered_diamond':
          ch_protein_db = params.uniref50_ec_filtered_diamond ? Channel.value(file(params.uniref50_ec_filtered_diamond)) : "null";
          break;

    case 'uniref90_ec_filtered_diamond':
          ch_protein_db = params.uniref90_ec_filtered_diamond ? Channel.value(file(params.uniref90_ec_filtered_diamond)) : "null";
          break;

    default:
          ch_protein_db = params.protein_db ? Channel.value(file(params.protein_db)) : "null";
          break;
  }
}

// metaphlan db is always used or custome, but always present
switch(params.metaphlan_db){
  case 'bowtie':
        ch_metaphlan_db = params.bowtie ? Channel.value(file(params.bowtie)) : "null";
        break;
  default:
        ch_metaphlan_db = params.metaphlan_db ? Channel.value(file(params.metaphlan_db)) : "null";
        params.mpaType = 'custom';
        break;
}



/* ############################################
 * Create a channel for input read files
 * ############################################
 */

inputSample = Channel.empty()
if (params.input) {
  tsvFile = file(params.input)
  inputSample = readInputFile(tsvFile)
}
else {
  log.info "No TSV file"
  exit 1, 'No sample were defined, see --help'
}

// split the channel into reading processes

(inputSampleFastqc, inputSampleHumann2, inputSampleMetaphlan2) = inputSample.into(3)


// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName

summary['Input']            = params.input
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
summary['Analysis Type']    = params.tool
summary['Metaphlan DB type'] = params.mpaType
summary['Nucleotide DB']    = params.nucleotide_db
summary['Protein DB']       = params.protein_db
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-scoop-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nibscbioinformatics/scoop Workflow Summary'
    section_href: 'https://github.com/nibscbioinformatics/scoop'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
        saveAs: { filename ->
                      if (filename.indexOf(".csv") > 0) filename
                      else null
                }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file "software_versions.csv"

    script:
    // Get all tools to print their version number here
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$idSample"
    label 'process_medium'
    publishDir "${params.outdir}/${idSample}/fastqc", mode: 'copy',
        saveAs: { filename ->
                      filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
                }

    input:
    set idSample, gender, status, file(read1), file(read2) from inputSampleFastqc

    output:
    file "*_fastqc.{zip,html}" into ch_fastqc_results

    script:
    """
    fastqc --quiet --threads $task.cpus $read1
    fastqc --quiet --threads $task.cpus $read2
    """
}

/*
 * STEP 2 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file (multiqc_config) from ch_multiqc_config
    file (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
    // Add in log files from your new processes for MultiQC to find!
    file ('fastqc/*') from ch_fastqc_results.collect().ifEmpty([])
    file ('software_versions/*') from ch_software_versions_yaml.collect()
    file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")

    output:
    file "*multiqc_report.html" into ch_multiqc_report
    file "*_data"
    file "multiqc_plots"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    // OK for now: Specify which MultiQC modules to use with -m for a faster run time
    """
    multiqc -f $rtitle $rfilename $custom_config_file .
    """
}

ch_nucleotide_db = ch_nucleotide_db.dump(tag:'nucleotide db')

process prepNucleotideDB {

  tag "prepare nucletide db"
  label "process_small"

  input:
  file(nucleotide_database) from ch_nucleotide_db

  output:
  path("nucleotidedb", type: 'dir' ) into ch_nucleotidedb_ready

  when: params.tool == 'humann2'

  script:
  """
  mkdir nucleotidedb
  mv ${nucleotide_database} nucleotidedb/.
  cd nucleotidedb

  tar -xvzf ${nucleotide_database}

  """
}

ch_protein_db = ch_protein_db.dump(tag:'protein db')

process prepProteinDB {
  tag "prepare protein db"
  label "process_small"

  input:
  file(protein_database) from ch_protein_db

  output:
  path("proteindb", type: 'dir' ) into ch_proteindb_ready

  when: params.tool == 'humann2'

  script:
  """
  mkdir proteindb
  mv ${protein_database} proteindb/.
  cd proteindb

  tar -xvzf ${protein_database}
  """


}


process prepMetaphlanDB {
  tag "prepare metaphlan db"
  label "process_medium"

  input:
  file(metaphlan_database) from ch_metaphlan_db
  file(md5) from Channel.value(file(params.mpamdd5))

  output:
  path("mpadb", type: 'dir' ) into (ch_metaphlandb_ready_mpaonly, ch_metaphlandb_ready_humann2)

  when: params.mpaType == 'default'

  script:
  """
  mkdir mpadb
  mv ${md5} mpadb/.
  mv ${metaphlan_database} mpadb/.
  cd mpadb

  tar -xvf ${metaphlan_database}
  bzip2 -d mpa_v20_m200.fna.bz2
  bowtie2-build --threads ${task.cpus} mpa_v20_m200.fna ./mpa_v20_m200

  """
}

//process prepMetaphlanDBCustom {
  // ################################
  // this still needs to be developed
  // TO-DO
  // ################################
//}



process metaphlanOnly {

  tag "metaphlan2 $idSample"

  label 'process_high'
  publishDir "${params.outdir}/${idSample}/metaphlan2", mode: 'copy'

  input:
  set idSample, gender, status, file(read1), file(read2) from inputSampleMetaphlan2
  path(mpadb) from ch_metaphlandb_ready_mpaonly

  output:
  file("${idSample}_metaphlan_bugs_list.tsv") into ch_metaphlan_results

  when: params.tool == 'metaphlan2'

  script:

  """
  cat ${read1} ${read2} >${idSample}_concat.fastq.gz

  metaphlan2.py \
  --input_type fastq \
  --tmp_dir=. \
  --bowtie2out=${idSample}_bt2out.txt \
  --bowtie2db ${mpadb} \
  --index ${params.mpa_index} \
  --nproc ${task.cpus} \
  ${idSample}_concat.fastq.gz \
  ${idSample}_metaphlan_bugs_list.tsv\
  """

}



process mergeMetaphlanOnly {

  label 'process_small'
  tag 'metaphlan2 merging'

  publishDir "${params.outdir}/merged_abundance/", mode: 'copy'

  input:
  file(bugList) from ch_metaphlan_results.collect()

  output:
  file("merged_abundance_table.tsv")

  when: params.tool == 'metaphlan2'


  script:
  """
  merge_metaphlan_tables.py \
  ${bugList} \
  > merged_abundance_table.tsv
  """

}



/*
* ################################################
* ## HUMANN2 PROCESSES ###########################
* ################################################
*/


process characteriseReads {

  tag "humann2 $idSample"
  cpus 8
  queue 'WORK'
  time '360h'
  memory '32 GB'

  publishDir "${params.outdir}", mode: 'copy'

  when: params.tool == 'humann2'

  input:
  set idSample, gender, status, file(read1), file(read2) from inputSampleHumann2
  path(nucleotide_db) from ch_nucleotidedb_ready
  path(protein_db) from ch_proteindb_ready
  path(mpadb) from ch_metaphlandb_ready_humann2

  output:
  file("${idSample}/*genefamilies.tsv") into gene_families_ch
  file("${idSample}/*pathabundance.tsv") into path_abundance_ch
  file("${idSample}/*pathcoverage.tsv")
  file("${idSample}/${idSample}_concat_humann2_temp/${idSample}_concat_metaphlan_bowtie2.txt")
  file("${idSample}/${idSample}_concat_humann2_temp/${idSample}_concat_metaphlan_bugs_list.tsv")

  script:

  """
  cat ${read1} ${read2} >${idSample}_concat.fastq.gz

  humann2 \
  --input ${idSample}_concat.fastq.gz \
  --nucleotide-database ${nucleotide_db} \
  --protein-database ${protein_db} \
  --metaphlan-options \"--bowtie2db ${mpadb} --index ${params.mpa_index}\" \
  --output ${idSample} \
  --threads ${task.cpus}
  """
}


process joinGenes {

  tag "humann2 join genes"
  label 'process_low'

  publishDir "${params.outdir}", mode: 'copy'

  when: params.tool == 'humann2'

  input:
  file genetables from gene_families_ch.collect()

  output:
  file("joined_genefamilies.tsv")
  file("joined_genefamilies_renorm_cpm.tsv")

  script:
  """
  humann2_join_tables \
  -i ./ \
  -o joined_genefamilies.tsv \
  --file_name genefamilies

  humann2_renorm_table \
  -i joined_genefamilies.tsv \
  -o joined_genefamilies_renorm_cpm.tsv \
  --units cpm

  """

}


process joinPathways {

  tag "humann2 join pathways"
  label 'process_low'

  publishDir "${params.outdir}", mode: 'copy'

  when: params.tool == 'humann2'

  input:
  file pathtables from path_abundance_ch.collect()

  output:
  file("joined_pathabundance.tsv")
  file("joined_pathabundance_renorm_cpm.tsv")

  script:
  """
  humann2_join_tables \
  -i ./ \
  -o joined_pathabundance.tsv \
  --file_name pathabundance

  humann2_renorm_table \
  -i joined_pathabundance.tsv \
  -o joined_pathabundance_renorm_cpm.tsv \
  --units cpm

  """

}








/*
 * STEP 3 - Output Description HTML

process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}
*/

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nibscbioinformatics/scoop] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nibscbioinformatics/scoop] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nibscbioinformatics/scoop] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nibscbioinformatics/scoop] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nibscbioinformatics/scoop] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            [ 'mail', '-s', subject, email_address ].execute() << email_txt
            log.info "[nibscbioinformatics/scoop] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nibscbioinformatics/scoop]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nibscbioinformatics/scoop]${c_red} Pipeline completed with errors${c_reset}-"
    }

}


def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nibscbioinformatics/scoop v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}

def readInputFile(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            def idSample  = row[0]
            def gender     = row[1]
            def status     = row[2].toInteger()
            def file1      = returnFile(row[3])
            def file2      = "null"
            if (hasExtension(file1, "fastq.gz") || hasExtension(file1, "fq.gz")) {
                checkNumberOfItem(row, 5)
                file2 = returnFile(row[4])
                if (!hasExtension(file2, "fastq.gz") && !hasExtension(file2, "fq.gz")) exit 1, "File: ${file2} has the wrong extension. See --help for more information"
            }
            // else if (hasExtension(file1, "bam")) checkNumberOfItem(row, 5)
            // here we only use this function for fastq inputs and therefore we suppress bam files
            else "No recognisable extension for input file: ${file1}"
            [idSample, gender, status, file1, file2]
        }
}

// #### SAREK FUNCTIONS #########################
def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "Malformed row in TSV file: ${row}, see --help for more information"
    return true
}

def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Return file if it exists
def returnFile(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    return file(it)
}

// Return status [0,1]
// 0 == Control, 1 == Case
def returnStatus(it) {
    if (!(it in [0, 1])) exit 1, "Status is not recognized in TSV file: ${it}, see --help for more information"
    return it
}

// ############### OTHER UTILS ##########################

// Example usage: defaultIfInexistent({myVar}, "default")
def defaultIfInexistent(varNameExpr, defaultValue) {
    try {
        varNameExpr()
    } catch (exc) {
        defaultValue
    }
}
