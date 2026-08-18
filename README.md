# MycoAssembler 🍄 → 🧬 → 💻 
This is a wrapper fungal assembly pipeline developed for the Borneman lab utilizing the following programs: 

* **[AAFTF](https://github.com/stajichlab/AAFTF)** 
* **[funannotate](https://github.com/nextgenusfs/funannotate)**

It takes an input metadata file describing the data and options for each genome, handles assembly and annotation of long and/or short reads, and optionally utilizes RNA-seq data for training in funannotate. The pipeline is configured to run on HPCC and utilize the modules installed there, with some required modifications described below.

## Basic Usage 🏁
There is a folder called `put_in_path` that contains the files you need to add to your `$PATH` to run the pipeline. You also need to make the modifications described below under **[Required Modifications and Notes](#important-required-modifications-and-notes-️)**.

Once everything is in your path, there are three options that need input:

1. Create a **tab-separated input file**. See [`example_infiles.md`](https://github.com/bpeacock44/FungiAssemblyPipeline/blob/main/example_infiles.md) for instructions and examples.

2. See the **[CodingQuarry](#codingquarry)** section below for information on what to put here. (It's **required**.)

3. Optionally create a `slurm.opts` file to specify slurm resources. (There are defaults.) See [`example_slurm.opts`](https://github.com/bpeacock44/FungiAssemblyPipeline/blob/main/example_slurm.opts) for an example.

```sh
MycoForge.sh --input infile.tsv \
    --codingquarry /path/to/CodingQuarry_v1.3 \
    --slurm-opts slurm.opts   # slurm.opts is optional
```

The pipeline will first validate all lines in your input file and notify you of any inconsistencies (e.g. missing data, incompatible assembler choices). Once validated, each line is submitted as a separate batch job to the cluster.

If you do not provide a `slurm.opts` file, the default resources are:

* 32 CPUs
* 128 GB RAM
* 14 days walltime

These defaults have been sufficient for all genomes tested so far. If you want email notifications or custom job names, use `--slurm-opts`. The default job name is set to "fass_#" where # is the line in your infile. If you don't set the job name in your slurm.opts file it'll use this default.

## REQUIRED Modifications and Helpful Notes 🪚
Some issues arise due to module configuration on HPCC. These are roughly documented in
[`issues_log.sh`](https://github.com/bpeacock44/FungiAssemblyPipeline/blob/main/archive/issues_log.sh).
Below is a summary of what you need to do to run the pipeline successfully.

### RNA-seq data format (SRA downloads)
If you download RNA-seq data from SRA, **Trimmomatic expects specific FASTQ headers**:

* Paired reads must include `1` and `2`
* The `+` line must match the header

Use the following `fastq-dump` format to ensure compatibility:

```sh
fastq-dump \
  --defline-seq '@$sn[_$rn]/$ri' \
  --defline-qual '+$sn[_$rn]/$ri' \
  --split-files SRR5488375
```

### PASA with SQLite
Running funannotate with SQLite on HPCC can fail. To resolve:

1. Create a personal PASA folder and symlink required files:

```sh
mkdir -p ~/.pasa/pasa_conf
cd ~/.pasa
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/bin bin
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/Launch_PASA_pipeline.pl Launch_PASA_pipeline.pl
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/misc_utilities misc_utilities
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/pasa-plugins pasa-plugins
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/PerlLib PerlLib
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/SAMPLE_HOOKS SAMPLE_HOOKS
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/schema schema
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/scripts scripts
cd pasa_conf
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/pasa_conf/pasa.alignAssembly.Template.txt pasa.alignAssembly.Template.txt
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/pasa_conf/pasa.annotationCompare.Template.txt pasa.annotationCompare.Template.txt
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/pasa_conf/pasa.CONFIG.template pasa.CONFIG.template
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/pasa_conf/README.conf README.conf
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/pasa_conf/sample_test.conf sample_test.conf
```

2. Create a custom config file:

```sh
cp $PASAHOME/pasa_conf/pasa.alignAssembly.Template.txt \
   $PASAHOME/pasa_conf/alignAssembly.sqlite.config
```

3. Edit to include:

```
DATABASE=pasa.sqlite

# Example parameters
validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=<__MIN_PERCENT_ALIGNED__>
validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=<__MIN_AVG_PER_ID__>
subcluster_builder.dbi:-m=50
```

Note that in the pipeline, the following variables will be automatically set. If you haven't done steps 1-3, the pipeline will fail!
```sh
export PASAHOME=~/.pasa
export PASACONF=$PASAHOME/pasa_conf/alignAssembly.sqlite.config
```

### CodingQuarry
Funannotate apparently no longer works with CodingQuarry directly. Use v1.3:

1. Download [CodingQuarry v1.3](https://sourceforge.net/projects/codingquarry/files/)
2. Extract: `tar -xvzf CodingQuarry_v1.3.tar.gz`
3. Build: `cd CodingQuarry_v1.3 && make`

The pipeline requires the path to CodingQuarry as an option: `--codingquarry /path/to/CodingQuarry_v1.3`

### Eggnog-mapper
Eggnog-mapper has been failing for me within funannotate. My current recommended approach is to run it manually using `eggnog_only.py`, which automates placement of results. Ensure it is in your `$PATH`. This will be done automatically in the pipeline and you don't need to do anything particular to set it up.

