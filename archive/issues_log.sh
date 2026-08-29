# trimmomatic expected headers of RNAseq files to include "1 and 2" for each pair and some were just numbered (e.g. each pair was given a number between the two files) when I fixed the files, it kept trying to use the old version.
# So now the code in the script downloads them correctly but make sure to delete old output files!

#RNAseq has to be formatted in a very specific way and downloading from SRA causes issues. This is why we use the following:
#fastq-dump --defline-seq '@$sn[_$rn]/$ri' --defline-qual '+$sn[_$rn]/$ri' --split-files SRR5488375

# - - 

# # # # # NEXT ISSUE # # # # #

# - -


# Another error: "$PASAHOME environmental variable not found, PASA is not properly configured.  You can use the --PASAHOME argument to specifiy a path at runtime" - this one is related to how the module is configured. Funannotate can use the PASA installation but you have to point it directly there via "export PASAHOME=/opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2" after you load the funannotate module
# now I am getting: 
#-------------------------------------------------------
#Traceback (most recent call last):
  #File "/opt/linux/rocky/8.x/x86_64/pkgs/funannotate/1.8.x/bin/funannotate", line 8, in <module>
    #sys.exit(main())
  #File "/bigdata/operations/pkgadmin/opt/linux/centos/8.x/x86_64/pkgs/funannotate/1.8.x/lib/python3.8/site-packages/funannotate/funannotate.py", line 717, in main
    #mod.main(arguments)
  #File "/bigdata/operations/pkgadmin/opt/linux/centos/8.x/x86_64/pkgs/funannotate/1.8.x/lib/python3.8/site-packages/funannotate/train.py", line 1178, in main
    #runPASAtrain(genome,
  #File "/bigdata/operations/pkgadmin/opt/linux/centos/8.x/x86_64/pkgs/funannotate/1.8.x/lib/python3.8/site-packages/funannotate/train.py", line 425, in runPASAtrain
    #if os.environ['PASACONF']:
  #File "/bigdata/operations/pkgadmin/opt/linux/centos/8.x/x86_64/pkgs/funannotate/1.8.x/lib/python3.8/os.py", line 675, in __getitem__
    #raise KeyError(key) from None
#KeyError: 'PASACONF'

# you essentially need to recreate the PASA folder in your own account and symlink to the official one and then create a config file in there.

# set up symlinks
mkdir ~/.pasa
cd ~/.pasa
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/bin bin
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/Launch_PASA_pipeline.pl Launch_PASA_pipeline.pl
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/misc_utilities misc_utilities
mkdir -p pasa_conf
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/pasa-plugins pasa-plugins
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/PerlLib PerlLib
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/SAMPLE_HOOKS SAMPLE_HOOKS
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/schema schema
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/scripts scripts
cd ~/.pasa/pasa_conf
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/pasa_conf/pasa.alignAssembly.Template.txt pasa.alignAssembly.Template.txt
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/pasa_conf/pasa.annotationCompare.Template.txt pasa.annotationCompare.Template.txt
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/pasa_conf/pasa.CONFIG.template pasa.CONFIG.template
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/pasa_conf/README.conf README.conf
ln -s /opt/linux/rocky/8.x/x86_64/pkgs/PASA/2.5.2/pasa_conf/sample_test.conf sample_test.conf

# create the config file
cp $PASAHOME/pasa_conf/pasa.alignAssembly.Template.txt \
   $PASAHOME/pasa_conf/alignAssembly.sqlite.config 

# And modify to include the following:

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#

## templated variables to be replaced exist as <__var_name__>
# database settings
DATABASE=pasa.sqlite

#######################################################
# Parameters to specify to specific scripts in pipeline
# create a key = "script_name" + ":" + "parameter"
# assign a value as done above.

#script validate_alignments_in_db.dbi
validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=<__MIN_PERCENT_ALIGNED__>
validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=<__MIN_AVG_PER_ID__>

#script subcluster_builder.dbi
subcluster_builder.dbi:-m=50

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#

# must change these variables before running part 4:
module load funannotate
export PASACONF=$PASAHOME/pasa_conf/alignAssembly.sqlite.config

# - - 

# # # # # NEXT ISSUE # # # # #

# - -

#[Jan 03 08:54 PM]: CMD ERROR: CodingQuarry -p 32 -f /bigdata/bornemanlab/bpeacock/PN113_long_read_assembly/outgroups/funannotate_genomes/acst.flye/acst.funannotate04/predict_misc /genome.softmasked.fa -t /bigdata/bornemanlab/bpeacock/PN113_long_read_assembly/outgroups/funannotate_genomes/acst.flye/acst.funannotate04/predict_misc/stringtie.gff3
#There is a long discussion of this issue here: https://github.com/nextgenusfs/funannotate/issues/807
#Go to SourceForge
#Download CodingQuarry_v1.3.tar.gz from Files to a user accessible directory in your computer
#Extract the compressed file using tar -xvzf CodingQuarry_v1.3.tar.gz
#Change directory cd CodingQuarry_v1.3 and use make
#Now you have a executable for CodingQuarry v1.3 called CodingQuarry
# After module load funannotate:
export PATH=/rhome/bpeacock/bigdata/programs/CodingQuarry_v1.3:$PATH

# - - 

# # # # # NEXT ISSUE # # # # #

# - -

#[Jan 08 08:34 AM]: Running Eggnog-mapper
#-------------------------------------------------------
#Traceback (most recent call last):
  #File "/opt/linux/rocky/8.x/x86_64/pkgs/funannotate/1.8.x/bin/funannotate", line 8, in <module>
    #sys.exit(main())
  #File "/bigdata/operations/pkgadmin/opt/linux/centos/8.x/x86_64/pkgs/funannotate/1.8.x/lib/python3.8/site-packages/funannotate/funannotate.py", line 720, in main
    #mod.main(arguments)
  #File "/bigdata/operations/pkgadmin/opt/linux/centos/8.x/x86_64/pkgs/funannotate/1.8.x/lib/python3.8/site-packages/funannotate/annotate.py", line 1041, in main
    #if packaging.version.parse(get_emapper_version()) >= packaging.version.parse("2.1.0"):
#NameError: name 'packaging' is not defined

# The only way I could figure out to get past this was to run eggnog mapper manually before running annotate.py. 
# So I wrote a script called "eggnog_only.py" and changed part 4 so it runs it before even starting annotate and 
# puts the results in the correct location so the annotate command can find them. 

# - - 

# # # # # NEXT ISSUE # # # # #

# - -

#[Jan 14 09:50 PM]: Evidence modeler has failed, exiting
#Traceback (most recent call last):
  #File "/opt/linux/rocky/8.x/x86_64/pkgs/funannotate/1.8.x/bin/funannotate", line 8, in <module>
    #sys.exit(main())
  #File "/bigdata/operations/pkgadmin/opt/linux/centos/8.x/x86_64/pkgs/funannotate/1.8.x/lib/python3.8/site-packages/funannotate/funanno
#tate.py", line 720, in main
    #mod.main(arguments)
  #File "/bigdata/operations/pkgadmin/opt/linux/centos/8.x/x86_64/pkgs/funannotate/1.8.x/lib/python3.8/site-packages/funannotate/predict
#.py", line 2644, in main
    #os.remove(EVM_out)
#FileNotFoundError: [Errno 2] No such file or directory: '/rhome/bpeacock/bigdata/PN113_long_read_assembly/running_all_again/04_funannot
#ate_genomes/DoUCR50.canu/DoUCR50.funannotate04/predict_misc/evm.round1.gff3'
#Prediction exists; skipping.
#Error: Prediction stats JSON (/rhome/bpeacock/bigdata/PN113_long_read_assembly/running_all_again/04_funannotate_genomes/DoUCR50.canu/Do
#UCR50.funannotate04/predict_results/Hyalorbilia_sp_DoUCR50.stats.json) not found or empty. Exiting.

# During predict, I got the above error. I tried running it again without changing anything and the error didn't occur. So give it another go if you run into this one.

# - - 

# # # # # NEXT ISSUE # # # # #

# - -

#Once I added the update step (following predict) I started getting this error only when I had the PASAHOME/PASACONF modifications described earlier. 
#[Jan 14 07:17 PM]: CMD ERROR: /rhome/bpeacock/.pasa/Launch_PASA_pipeline.pl -c /rhome/bpeacock/bigdata/PN113_long_read_assembly/running
#_all_again/04_funannotate_genomes/acst.flye/acst.funannotate04/update_misc/pasa/annotCompare.txt -g /rhome/bpeacock/bigdata/PN113_long_
#read_assembly/running_all_again/04_funannotate_genomes/acst.flye/acst.funannotate04/update_misc/genome.fa -t /rhome/bpeacock/bigdata/PN
#113_long_read_assembly/running_all_again/04_funannotate_genomes/acst.flye/acst.funannotate04/update_misc/trinity.fasta.clean -A -L --CP
#U 2 --annots /rhome/bpeacock/bigdata/PN113_long_read_assembly/running_all_again/04_funannotate_genomes/acst.flye/acst.funannotate04/upd
#ate_misc/genome.gff3
#
# Turns out adding the config file (setting PASACONF) is needed for the train step but causes issues with the update step. The solution is to unset PASACONF before running update. 

# - - 

# # # # # NEXT ISSUE # # # # #

# - -

# Getting the same error as the previous set. But I noticed the following:
# [01/28/26 12:56:50]: 2 transcripts derived from Trinity
# Apparently I was using bad data that couldn't be assembled. Once I changed the data it ran without issue. 