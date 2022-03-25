fast5 = Channel.fromPath(params.samples).splitCsv(sep: '\t', header: true)
genome = Channel.fromPath(params.genome)
annotation = Channel.fromPath(params.annotation)

params.help = false

if (params.help) {
    log.info """
Usage:

nextflow run main.nf [input parameters] [system parameters]

Input:
--samples PATH                  Sample definitions file [${params.samples}]
--fast5_root PATH               Main folder that contains subfolders of FAST5 files for each run listed in sample definitions file
--genome STR                    Reference genome FASTA file
--annotation STR                Genome annotation GTF file

System parameters:
--threads INT                   Number of threads available to pipeline [${params.threads}]

Guppy parameters:
--gpu_device STR                GPU device to be used by Guppy [${params.gpu_device}]
--num_barcoding_buffers INT     Number of barcoding buffers to be used by Guppy [${params.num_barcoding_buffers}]

Minimap2 parameters:
--k INT                         K-mer size [${params.k}]
--w INT                         Minimizer window size [${params.w}]
--I STR                         Split index for every ~STR input bases [${params.I}]
--K STR                         Minibatch size for mapping [${params.K}]

FeatureCounts parameters:
--extra_attributes STR          Extra attributes from annotation GTF file to include in counting output, comma separated ['${params.extra_attributes}']
--Q INT                         The minimum mapping quality score a read must satisfy in order to be counted [${params.Q}]

NanoPlot parameters:
--maxlength INT                 Ignore reads longer than INT bases [${params.maxlength}]
"""
    exit 0
}

process defineSamples {
    input:
        val runlist from fast5
    
    output:
        tuple run_name, trim_strategy, fast5_path into samples_basecall
        tuple run_name, q_cutoff into samples_filterReads
        tuple run_name, barcode_kit, barcode_list into samples_demultiplex

    exec:
        run_name = runlist.Run
        q_cutoff = runlist.Q_cutoff
        barcode_kit = runlist.Kit_barcode
        trim_strategy = barcode_kit.size() ? "none" : "dna"
        fast5_path = "${params.fast5_root}/${run_name}"

        barcode_list = [
            runlist.Barcode_00,
            runlist.Barcode_01,
            runlist.Barcode_02,
            runlist.Barcode_03,
            runlist.Barcode_04,
            runlist.Barcode_05,
            runlist.Barcode_06,
            runlist.Barcode_07,
            runlist.Barcode_08,
            runlist.Barcode_09,
            runlist.Barcode_10,
            runlist.Barcode_11,
            runlist.Barcode_12
        ]
}

process basecall {
    tag "$run_name"

    module params.mod_guppy
    cpus 3
    maxForks 1
    memory '3GB'

    stageInMode 'copy'
    afterScript "rm -rf ${run_name}"

    input:
        tuple run_name, trim_strategy, path(fast5_path) from samples_basecall
    
    output:
        tuple run_name, "*.fastq.gz" into fastq_files
        tuple run_name, "sequencing_summary.txt" into sequencing_summary_table
    
    shell:
        '''
        guppy_basecaller \
            --input_path $PWD/!{run_name} \
            --recursive \
            --save_path $PWD \
            --config dna_r9.4.1_450bps_sup.cfg \
            --device !{params.gpu_device} \
            --disable_qscore_filtering \
            --trim_strategy !{trim_strategy} \
            --compress_fastq \
            --calib_detect \
            --disable_pings
        '''
}

sequencing_summary_table
    .into{ sequencing_summary_table_export; sequencing_summary_table_nanoplot}

samples_basecall_filterReads = samples_filterReads.join(fastq_files)

process filterReads {
    tag "$run_name"

    module params.mod_nanopack
    publishDir "data/fastq", mode: 'copy'
    cpus 1
    memory '100MB'

    input:
        tuple run_name, q_cutoff, path(fastq_files) from samples_basecall_filterReads

    output:
        tuple run_name, "${run_name}_unfiltered.fastq.gz" into unfiltered_fastq
        tuple run_name, "${run_name}_filtered.fastq.gz" into filtered_fastq
        path "${run_name}_fastq_files.txt"

    shell:
        '''
        ls *.fastq.gz > !{run_name}_fastq_files.txt
        cat *.fastq.gz > !{run_name}_unfiltered.fastq.gz

        gunzip -c !{run_name}_unfiltered.fastq.gz | \
        NanoFilt --quality !{q_cutoff} | \
        gzip -9 > !{run_name}_filtered.fastq.gz
        '''
}

samples_basecall_demultiplex = samples_demultiplex.join(filtered_fastq)

process demultiplex {
    tag "$run_name"
    
    module params.mod_guppy
    cpus 3
    maxForks 1
    memory '64GB'

    input:
        tuple run_name, barcode_kit, barcode_list, path(fastq_files) from samples_basecall_demultiplex
    
    output:
        tuple run_name, "${run_name}/barcode*", barcode_list into fastq_demultiplexed_folders
        tuple run_name, "${run_name}/barcoding_summary.txt" optional true into barcoding_summary_table
    
    shell:
    if( barcode_kit )
        '''
        guppy_barcoder \
            --input_path ${PWD} \
            --save_path ${PWD}/!{run_name} \
            --barcode_kits !{barcode_kit} \
            --device !{params.gpu_device} \
            --num_barcoding_buffers !{params.num_barcoding_buffers} \
            --worker_threads !{task.cpus} \
            --detect_mid_strand_adapter \
            --detect_mid_strand_barcodes \
            --trim_barcodes \
            --compress_fastq
        
        mkdir -p !{run_name}/barcode{00..12}
        '''
    else
        '''
        mkdir -p !{run_name}/barcode00
        cp -l *.fastq.gz !{run_name}/barcode00/
        '''
}

fastq_demultiplexed_sample_map = fastq_demultiplexed_folders
    .transpose()
    .filter{ it[2].length() }

process mergeDemultiplexedFastq {
    tag "$run_name $sample_name"
    publishDir "data/fastq", mode: 'copy'
    cpus 1
    memory '10MB'

    input:
        tuple run_name, path(barcode_folder), sample_name from fastq_demultiplexed_sample_map

    output:
        tuple run_name, sample_name, "${run_name}_${sample_name}.fastq.gz" into merged_filtered_demultiplexed
    
    shell:
        '''
        cat barcode*/*.fastq.gz > !{run_name}_!{sample_name}.fastq.gz
        '''
}

merged_filtered_demultiplexed.into {
    merged_filtered_demultiplexed_align;
    merged_filtered_demultiplexed_nanoplot
}

process buildIndex {
    tag "$genome.baseName"
    module params.mod_minimap2
    cpus 3
    memory '10GB'

    input:
        path genome

    output:
        path "${genome.baseName}.mmi" into genome_index

    shell:
        '''
        minimap2 \
            -x splice \
            -k !{params.k} \
            -w !{params.w} \
            -I !{params.I} \
            -d !{genome.baseName}.mmi \
            !{genome}
        '''
}

process align {
    tag "$run_name $sample_name"
    
    module params.mod_minimap2
    module params.mod_samtools
    publishDir "results/bam", mode: 'copy'
    cpus 6
    memory '15GB'

    input:
        tuple run_name, sample_name, path(fastq_file) from merged_filtered_demultiplexed_align
        each path(genome_index)

    output:
        tuple run_name, sample_name, "${run_name}_${sample_name}.bam" into bam_files
        tuple run_name, sample_name, "${run_name}_${sample_name}.bam.bai" into bai_files

    shell:
        '''
        minimap2 \
            -x splice \
            -K !{params.K} \
            -a \
            -2 \
            --MD \
            --cs=long \
            --secondary=no \
            -t !{task.cpus} \
            !{genome_index} \
            !{fastq_file} | \
        samtools view \
            -@ !{task.cpus} \
            -b \
            -u \
            - | \
        samtools sort \
            -@ !{task.cpus} \
            -l 9 \
            -o !{run_name}_!{sample_name}.bam \
            -

        samtools index \
            -@ !{task.cpus} \
            !{run_name}_!{sample_name}.bam
        '''
}

bam_files.into {
    bam_files_featureCounts;
    bam_files_alignmentStats;
    bam_files_nanoplot
}

process countFeatures {
    tag "$run_name $sample_name $mode"
    
    module params.mod_subread
    publishDir "data/featureCounts", mode: 'copy'
    cpus 1
    memory '300MB'

    input:
        tuple run_name, sample_name, path(bam_file) from bam_files_featureCounts
        each path(annotation)
        each mode from Channel.fromList( [0, 1, 2] )
    
    output:
        tuple run_name, sample_name, "${run_name}_${sample_name}_s${mode}_featureCounts.tsv", "${run_name}_${sample_name}_s${mode}_featureCounts.tsv.summary" into featureCounts

    shell:
        '''
        featureCounts \
            -o !{run_name}_!{sample_name}_s!{mode}_featureCounts.tsv \
            -s !{mode} \
            -a !{annotation} \
            -T !{task.cpus} \
            -Q !{params.Q} \
            --extraAttributes !{params.extra_attributes} \
            -L \
            --largestOverlap \
            --primary \
            !{bam_file}
        '''
}

grouped_featureCounts = featureCounts.groupTuple(by: [0, 1])

process mergeFeatureCounts {
    tag "$run_name $sample_name"
    publishDir "results/counts", mode: 'copy'
    cpus 1

    input:
        tuple run_name, sample_name, path(featureCounts_table), tmp from grouped_featureCounts
    
    output:
        tuple run_name, sample_name, "${run_name}_${sample_name}_featureCounts.tsv" into merged_featureCounts
    
    shell:
        '''
        !{projectDir}/scripts/merge_featureCounts.awk *_featureCounts.tsv > !{run_name}_!{sample_name}_featureCounts.tsv
        '''
}

// Quality controls
process basecallQC {
    tag "${run_name}"

    module params.mod_nanopack
    publishDir "results/QC", mode: 'copy'
    cpus 1
    memory '5GB'

    errorStrategy 'ignore'

    input:
        tuple run_name, path(sequencing_summary_file) from sequencing_summary_table_nanoplot
    
    output:
        path "${run_name}_basecall"

    shell:
        '''
        NanoPlot \
            --threads !{task.cpus} \
            --summary sequencing_summary.txt \
            --outdir !{run_name}_basecall \
            --title !{run_name}_basecall \
            --plots kde dot \
            --dpi 600 \
            --maxlength !{params.maxlength} \
            --N50 \
            --raw \
            --huge
        '''
}

process fastqQC {
    tag "$run_name $sample_name"

    module params.mod_nanopack
    publishDir "results/QC", mode: 'copy'
    cpus 1
    memory '5GB'

    errorStrategy 'ignore'

    input:
        tuple run_name, sample_name, path(merged_sample_fastq) from merged_filtered_demultiplexed_nanoplot
    
    output:
        path "${run_name}_reads_${sample_name}"

    shell:
        '''
        NanoPlot \
            --threads !{task.cpus} \
            --fastq_rich !{merged_sample_fastq} \
            --outdir !{run_name}_reads_!{sample_name} \
            --title !{run_name}_!{sample_name}_reads \
            --plots kde dot \
            --dpi 600 \
            --maxlength !{params.maxlength} \
            --N50 \
            --raw \
            --huge
        '''
}

process bamQC {
    tag "$run_name $sample_name"

    module params.mod_nanopack
    publishDir "results/QC", mode: 'copy'
    cpus 1
    memory '5GB'

    errorStrategy 'ignore'

   input:
        tuple run_name, sample_name, path(bam_file) from bam_files_nanoplot
    
    output:
        path "${run_name}_alignment_${sample_name}"

    shell:
        '''
        NanoPlot \
            --threads !{task.cpus} \
            --bam !{bam_file} \
            --outdir !{run_name}_alignment_!{sample_name} \
            --title !{run_name}_!{sample_name}_alignment \
            --plots kde dot \
            --dpi 600 \
            --maxlength !{params.maxlength} \
            --N50 \
            --raw \
            --huge
        '''
}

// Stats
process alignmentStats {
    tag "$run_name $sample_name"

    module params.mod_samtools
    publishDir "data/stats", mode: 'copy'
    cpus 2

    input:
        tuple run_name, sample_name, path(bam_file) from bam_files_alignmentStats
    
    output:
        path "${run_name}_${sample_name}_alignmentStats.txt"

    shell:
        '''
        samtools stats -@ !{task.cpus} !{bam_file} > !{run_name}_!{sample_name}_alignmentStats.txt
        '''
}

process extractBasecallingSummary {
    tag "$run_name"
    
    publishDir "data/guppy", mode: 'copy'
    cpus 1
    memory '300MB'

    input:
        tuple run_name, path(sequencing_summary) from sequencing_summary_table_export
    
    output:
        path "${run_name}_sequencing_summary.txt.gz"
    
    shell:
        '''
        gzip -9 -c sequencing_summary.txt > !{run_name}_sequencing_summary.txt.gz
        '''
}

process extractBarcodingSummary {
    tag "$run_name"
    
    publishDir "data/guppy", mode: 'copy'
    cpus 1
    memory '300MB'

    input:
        tuple run_name, path(barcoding_summary) from barcoding_summary_table
    
    output:
        path "${run_name}_barcoding_summary.txt.gz"
    
    shell:
        '''
        gzip -9 -c barcoding_summary.txt > !{run_name}_barcoding_summary.txt.gz
        '''
}