# CellRanger step-by-step workflow

> [!NOTE]
> For DeMayo lab members, you will need to use `scigate` for this.

* _Reference: https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-ct_
* Note that this procedure is specific to the DeMayo lab and might need to be altered for others.
* For the sake of example, let's use the following names:
    * Sequencing run name: `SequencingRun`
    * Sample 1 name: `Sample1`
    * Sample 2 name: `Sample2`

## Installing Cell Ranger and reference files

Run the `installCellRanger.sh` file from the `tools` folder of this repository to have Cell Ranger installed locally on your home directory (`~`). This will create two folders:

1. `cellranger-8.0.1`, which will contain the main Cell Ranger application, version 8.0.1.
2. `refdata-gex-GRCh38_and_GRCm39-2024-A`, which contains the reference genomes needed for Cell Ranger for mice (GRCm39) and humans (GRCh38).

It will also add Cell Ranger to your PATH so you can call it from your command line _after
restarting your shell._

## Preparing the sequencing results for Cell Ranger

NOVAseq results typically come in the following format:

```
SequencingRun/
    fastq_path/
        Sample1_S1_L001_I1_001.fastq.gz
        Sample1_S1_L001_I2_001.fastq.gz
        Sample1_S1_L001_R1_001.fastq.gz
        Sample1_S1_L001_R2_001.fastq.gz
        Sample1_S1_L002_I1_001.fastq.gz
        Sample1_S1_L002_I2_001.fastq.gz
        Sample1_S1_L002_R1_001.fastq.gz
        Sample1_S1_L002_R2_001.fastq.gz
        Sample2_S1_L001_I1_001.fastq.gz
        Sample2_S1_L001_I2_001.fastq.gz
        Sample2_S1_L001_R1_001.fastq.gz
        Sample2_S1_L001_R2_001.fastq.gz
        Sample2_S1_L002_I1_001.fastq.gz
        Sample2_S1_L002_I2_001.fastq.gz
        Sample2_S1_L002_R1_001.fastq.gz
        Sample2_S1_L002_R2_001.fastq.gz
        ...
    ...
```

Note a couple of things:

1. The FASTQ files contain the terms 'R1', 'R2', 'I1', and 'I2', which
   distinguish them from combined RNAseq + ATACseq data (which typically
   replaces 'I2' with 'R3').
2. We _only_ want to use the FASTQ files in `fastq_path`, not anywhere else.
3. The sample names can only be found by the starting prefix of each FASTQ file
   in the `fastq_path` folder. For example, the files `Sample1_S1_L001_XX_001.fastq.gz` 
   define the sample `Sample1`.

With that in mind, we can now proceed to setting up the reference data.

## Configure and run Cell Ranger

Execute the following command _for each sample_ you want to analyze:

```sh
cellranger count --id="[ARBITRARY ID NAME]"
    --fastqs="SequencingRun/fastq_path" # or whatever the path is to the fastq_path folder
    --sample="Sample1" # or Sample2 or another name
    --create-bam="[true OR false]"
    --transcriptome="[REFERENCE TRANSCRIPTOME FOLDER]"
```

Some notes:

* `id` can be an arbitrary name that will become the name of the output folder.
* `transcriptome` is the path to the reference data downloaded earlier. If you
  used the install script at the top of this manual, it would likely be located
  at `~/refdata-gex-GRCh38_and_GRCm39-2024-A`.

That's it! Cell Ranger should run okay from this point on.

## Additional information

See the official documentation from Cell Ranger for more information on how to
use and interpret the results of Cell Ranger. Some links are below:

* Tutorials: https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-in
* Analysis: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis
* Manual pages: https://www.10xgenomics.com/support/software/cell-ranger/latest/resources/cr-command-line-arguments