# Using the CIBERSORTx tool to impute bulk cell fractions

Last updated by Peter Lais on 05/28/2024.

> [!NOTE] For DeMayo lab members, you will need to use `scigate` for this.

* Reference: CIBERSORTx documentation. Please sign up at
  https://cibersortx.stanford.edu/, navigate to the download page, and then
  download the access instructions folder. Relevant instructions will be
  contained within the 'CIBERSORTx-Fractions.txt' file of the zipped archive.
  This file contains documentation about how to use CIBERSORTx Fractions
  offline as well as informations about parameters you can use (for example,
  how to impute absolute vs. relative cell fractions; see the end of this
  tutorial document).

## General CIBERSORTx workflow for imputing fractions

CIBERSORTx has a number of modes, with `fractions` being one of them. This
allows us to use a **single-cell reference** file to call the proportions of
cells from a series of **bulk samples.** Before starting, we will need the
following two files. For reference, we will declare the following variables:

* $M$: The number of genes in both our bulk and single-cell experiment files.
* $N$: The number of cells in our single-cell experiment files.
* $P$: The number of bulk samples that we have.

With this in mind, we need the following files before beginning anything:

1. **A single-cell reference matrix.** This matrix should be of shape $(M,N)$,
   having one column for each cell with _raw_ (or, if scaled, at least 
   linear-space) counts for each of our $M$ queried genes. **These data
   should be pre-clustered; see the 'Getting Started with Seurat' vignette
   for more information.** The column names of this matrix should be the
   cluster from which each cell originated, and the row names should be the
   names of the genes being studied. Here is an example single-cell reference
   matrix:

   | | Cluster 1 | Cluster 1 | Cluster 2 | Cluster 3 | ... |
   | --- | --- | --- | --- | --- | --- |
   | Pgr | 1 | 23 | 0 | 1 | ... | 
   | Srf | 0 | 2 | 0 | 1 | ... | 
   | Trim28 | 1 | 3 | 2 | 1 | ... | 
   | ... | ... | ... | ... | ... | ... |

   This matrix should be saved as a tab-delimited text file (.txt in Excel,
   also can have the suffix '.tsv'). Note how more than one column can come
   from the same cluster; this represents the idea that there are multiple
   cells for each cluster present.

2. **A matrix of bulk samples.** This matrix should be of shape $(M,P)$, having
   one column for each _sample_ with _raw_ (or, if scaled, at least linear-space)
   counts for each of our $M$ queried genes. Here is an example bulk matrix:

   | | Sample 1 | Sample 2 | Sample 3 | Sample 4 | ... |
   | --- | --- | --- | --- | --- | --- |
   | Pgr | 5 | 2 | 1 | 1 | ... | 
   | Srf | 3 | 2 | 0 | 0 | ... | 
   | Trim28 | 0 | 0 | 2 | 1 | ... | 
   | ... | ... | ... | ... | ... | ... |

From here, there are two ways you can process these data.

### Process data online at https://cibersortx.stanford.edu/

See ['Tutorial 2'](https://cibersortx.stanford.edu/tutorial.php) on the main
CIBERSORTx website for information about how to upload your
two tab-delimited files onto the CIBERSORTx server and perform the
analysis. There will be two processing steps in total:

1. Convert the single-cell reference file into a signature matrix file,
   which has shape $(M,N_\textrm{unique})$. Here, $N_\textrm{unique}$
   denotes the number of unique cell types (i.e., number of clusters)
   present in the single-cell data.
2. Using the signature matrix, impute cell type fractions on the bulk
   sample matrix.

### Process data offline using Docker/Singularity

The main CIBERSORTx website sometimes doesn't cooperate properly with users.
If problems arise here—for example, your files are too big, the server isn't
responding, maintenance is being done—you can also do your analysis offline.
Follow the steps below.

1. For this analysis, we will use **Apptainer** and **conda**. The first step
will be to install each of these programs.
    
    a. Install conda by either following the tutorial [here](https://github.com/conda-forge/miniforge)
       or executing the
       `installMiniforge.sh` script in the `tools` folder of this repository.

    b. Install Apptainer using conda by executing the command: `mamba install
       conda-forge::apptainer`. If you don't have `mamba` installed, replace
       `mamba` with `conda` in the original command.

    c. Obtain an API token by following the instructions for running CIBERSORTx
       locally, present on the [Downloads](https://cibersortx.stanford.edu/download.php)
       page of CIBERSORTx's website.

    d. Install the CIBERSORTx fractions image on your machine using the command
       `apptainer pull docker://cibersortx/fractions`.

    e. Execute the following two commands in your target folder:

       ```sh
       mkdir CIBERSORTxOutput
       apptainer run
           --no-mount home,cwd
           --home :src/
           -B $(realpath .):/src/data
           -B $(realpath .)/CIBERSORTxOutput:/src/outdir
           docker://cibersortx/fractions
           --username [CIBERSORTX ACCOUNT USERNAME]
           --token [CIBERSORTX API ACCESS TOKEN]
           --single_cell TRUE
           --refsample /path/to/single_cell_reference_mtx.txt
           --mixture /path/to/bulk_mixture_matrix.txt
           --perm [NUMBER OF PERMUTATIONS FOR P-VALUES]
       ```

Results will be output in the `CIBERSORTxOutput` folder that you created.
Make sure that you clear this folder if you want to re-run your analysis.

## Expected results

You will have a number of output files once completing deconvolution/cell-type
imputation using CIBERSORTx, but the most important result is a single matrix
that takes the shape $(P,N_\textrm{unique})$:

| | Cluster 1 | Cluster 2 | Cluster 3 | Cluster 4 | ... |
| --- | --- | --- | --- | --- | --- |
| Sample 1 | 0.5 | 0.2 | 0.1 | 0.1 | ... | 
| Sample 2 | 0.3 | 2 | 0 | 0 | ... | 
| Sample 3 | 0 | 0 | 0.2 | 0.1 | ... | 
| ... | ... | ... | ... | ... | ... |

This tells you the **relative** proportions of each cell type across each of your
samples. The relative fractioning is done by default; see the reference sources
at the top of the slide for information about how to impute **absolute** cell fractions.