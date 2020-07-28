[General]
job_type = SGE
job_queue = bigmem

[Unzip]
input_fofn= input.fofn
input_bam_fofn=input_bam.fofn

#smrt_bin=/mnt/secondary/builds/full/3.0.0/prod/current-build_smrtanalysis/smrtcmds/bin/
#smrt_bin=/pbi/dept/secondary/builds/mainline/current_smrttools_prebuilt_installdir/smrtcmds/bin
#smrt_bin=/mnt/software/s/smrttools/4.0.0/private/otherbins/all/bin
#smrt_bin==/home/UNIXHOME/cdunn/work/hops/VENV/bin
smrt_bin=/pbi/dept/secondary/builds/5.1.0.SNAPSHOT/current_smrttools-release_installdir/private/otherbins/all/bin
jobqueue = bigmem
sge_phasing= -pe smp 24 -q %(jobqueue)s
sge_quiver= -pe smp 24 -q %(jobqueue)s
sge_track_reads= -pe smp 48 -q %(jobqueue)s
sge_blasr_aln=  -pe smp 24 -q %(jobqueue)s
sge_hasm=  -pe smp 64 -q %(jobqueue)s
unzip_blasr_concurrent_jobs = 60
unzip_phasing_concurrent_jobs = 60
quiver_concurrent_jobs = 60
