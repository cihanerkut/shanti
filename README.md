# ShaNTi
A Nextflow pipeline to process ONT raw data from cDNA sequencing experiments. This workflow is designed and will be developed in the context of cancer research.

# Before you begin
This workflow was developed on a workstation equipped with 32 virtual cores and an NVIDIA Tesla V100 16Gb GPU. GPU is required for efficient basecalling. All other processes are optimized for 32 cores. If you need to adapt it to your system, you should change the `main.nf` file accordingly.

In the future we will test the pipeline on different platforms and provide singularity / docker containers.

# Usage
The easiest way to use the pipeline is to modify the nextflow.config file and run the main workflow directly.

```{bash}
nextflow run main.nf
```
## Main parameters
You should probably change the following parameters to adapt it to your system:

* Module definitions (e.g. `mod_guppy`, `mod_minimap2`, etc.)
* `threads` depending on number of cores on your system
* Default paths
  * `samples`: Sample definitions as a tab-separated value (TSV) file
  * `fast5_root`: Highest level path to contain FAST5 files
  * `genome`: Path to the reference genome FASTA file
  * `annotation`: Path to the gene annotations GTF file

## Module definitions
We use Lmod to manage our software stack created using Easybuild, therefore it is integrated into the pipeline. We will provide the Easyconfig scripts here.

## Sample definitions file
A sample definitions file should normally look like this:

| Run | Kit_barcode | Q_cutoff | Barcode_00 | Barcode_01 | ... | Barcode_12 | 
| - | - | - | - | - | - | - | 
| 20210526_1353_MN35504_FAP07621_b3ba36f1 |  | 7 | Patient_001 |  | ... |  | 
| 20210611_1045_MN22813_FAP07780_2c9dcd4c | EXP-NBD104 | 7 |   | Patient_002 | ... | Patient_013 | 
| 20210622_1324_MN22813_FAP07698_c9f8bbf3 | SQK-PCB109 | 7 |   | Patient_014 | ... | Patient_025 | 

All columns must exist in the file. Among those, **Run** and **Q_cutoff** are required fields, i.e. you must supply a run name and q-value cutoff.

If the experiment doesn't have a multiplexed design, the name of the single sample must be provided under **Barcode_00** field. Otherwise, the barcoding / sequencing kit (e.g. EXP-NBD104, SQK-PCB109, etc.) must be provided under **Kit_barcode** and sample names under corresponding **Barcode_NN** fields.

For barcodes that are not used in a run, the corresponding cell must be left empty.

## FAST5 files
You should place all FAST5 files under folders with the same name defined under **Run** in the sample definitions file. Each such folder will be recursively searched for FAST5 files.

For each basecalling step, the folder will be copied under the `work` directory until basecalling finishes. Therefore we recommend to place only FAST5 files under these folder for faster processing.

## Genome and annotations
The genome FASTA and annotation GTF files must match. You can get matching files from GENECODE. These files may be GZIPped to save space.

## Questions and comments

Please create issues for your questions, comments, bug reports, feature requests, etc.