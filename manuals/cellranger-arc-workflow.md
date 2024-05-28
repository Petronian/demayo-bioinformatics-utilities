# Cell Ranger ARC 2.0.2 step-by-step workflow

Last updated by Peter Lais on 05/28/2024.

> [!NOTE]
> For DeMayo lab members, you will need to use `scigate` for this.

* _Reference: https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-ct_
* Note that this procedure is specific to the DeMayo lab and might need to be altered for others.
* For the sake of example, let's use the following names:
    * Sequencing run name: `SequencingRun`
    * Sample 1 name: `Sample1`
    * Sample 2 name: `Sample2`

## Installing Cell Ranger ARC and reference files

Run the `installCellRangerArc.sh` file (once downloaded, execute the command: `bash installCellRangerArc.sh`) from the `tools` folder of this repository to have Cell Ranger installed locally on your home directory (`~`). This will create two folders:

1. `cellranger-arc-2.0.2`, which will contain the main Cell Ranger ARC application, version 2.0.2.
2. `refdata-cellranger-arc-GRCh38-2020-A-2.0.0`, which contains the reference genomes needed by Cell Ranger ARC for humans (GRCh38).
2. `refdata-cellranger-arc-mm10-2020-A-2.0.0`, which contains the reference genomes needed by Cell Ranger ARC for mice (GRCm38/mm10).

It will also add Cell Ranger ARC to your PATH so you can call it from your command line _after
restarting your shell._

## Preparing the sequencing results for Cell Ranger ARC

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

1. The FASTQ files contain the terms 'R1', 'R2', 'R3', and 'I1', which
   distinguish them from solely RNAseq data (which typically
   replaces 'R3' with 'I2').
2. We _only_ want to use the FASTQ files in `fastq_path`, not anywhere else.
3. The sample names can only be found by the starting prefix of each FASTQ file
   in the `fastq_path` folder. For example, the files `Sample1_S1_L001_XX_001.fastq.gz` 
   define the sample `Sample1`.

With that in mind, we can now proceed to setting up the reference data.

## Configure and run Cell Ranger ARC

Execute the following command _for each sample_ you want to analyze:

```sh
cellranger-arc count --id="[ARBITRARY ID NAME]"
    --libraries="./libraries.csv" # Path to libraries file, SEE BELOW!
    --reference="[REFERENCE TRANSCRIPTOME FOLDER]"
```

Some notes:

* `id` can be an arbitrary name that will become the name of the output folder.
* `reference` is the path to the reference data downloaded earlier. If you
  used the install script at the top of this manual, it would likely be located
  at `~/refdata-gex-GRCh38_and_GRCm39-2024-A`.
* `libraries.csv` is a CSV file telling Cell Ranger ARC where the samples to
  process are. See below for more details.

The `libraries.csv` file should be a CSV file _in the same directory as where
you call Cell Ranger ARC_ (assuming its path is `./libraries.csv`) that looks
like the following:

```csv
/absolute/path/to/SequencingRun/fastq_path,Sample1,Chromatin Accessibility
/absolute/path/to/SequencingRun/fastq_path,Sample1,Gene Expression
... (make another file for sample 2 or include Sample2 lines here if desired)
```

This defines three columns, in order: the FASTQ file location (in our case,
`fastq_path`), the name of the sample being supplied, and what type of data are
available in the sample. Since we have a combined RNAseq and ATACseq run, we have
a sample named `Sample1` located at `/absolute/path/to/SequencingRun` with both
`Chromatin Accessibility` and `Gene Expression` data. We indicate such
information using two lines on our CSV file.

If desired, more samples can be included on the same CSV file, but the results
will all be aggregated together. As an alternative, the runs can be done
separately and `cellranger-arc aggr` can be used to aggregate the results. (See
links below for more information about that.)

That's it! Cell Ranger ARC should run okay from this point on.

## Additional information

See the official documentation from Cell Ranger for more information on how to
use and interpret the results of Cell Ranger. Some links are below:

* Tutorials: https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/tutorials
* Analysis: https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis