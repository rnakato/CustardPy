# from .hic
hic=CustardPyResults/Cooler_hg38/Control/hic/contact.bwa.q30.hic
# singularity exec --bind /work,/work2,/work3 /work/SingularityImages/custardpy.3.5.0.sif  \
#             mustache -f $hic -norm SCALE \
#	     -p 12 -r 5kb -pt 0.05 -o mustache.hic.5kb.tsv

#exit
# from .cool
cool=CustardPyResults/Cooler_hg38/Control/cool/contact.bwa.q30.mcool
singularity exec --bind /work,/work2,/work3 /work/SingularityImages/custardpy.3.5.0.sif  \
            mustache -f $cool \
	    -r 5kb -norm weight \
	    -pt 0.05 -p 12 \
	    -o mustache.cool.5kb.tsv
