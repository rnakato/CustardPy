#sing="singularity exec --bind /work,/work2,/work3 /work/SingularityImages/custardpy.2.2.2.sif"
sing="singularity exec custardpy.sif"

$sing hictk --version

hicdir=CustardPyResults_Hi-C/Cooler_hg38/Control

# Dump
$sing hictk dump $hicdir/cool/Control.25000.cool
$sing hictk dump $hicdir/cool/Control.multires.cool
$sing hictk dump $hicdir/hic/contact_map.q30.hic --table normalizations
$sing hictk dump $hicdir/hic/contact_map.q30.hic --table resolutions

$sing hictk dump $hicdir/cool/Control.25000.cool --join
$sing hictk dump $hicdir/cool/Control.25000.cool --cis-only --join
$sing hictk dump $hicdir/cool/Control.25000.cool --trans-only --join

# Validate
$sing hictk validate $hicdir/cool/Control.multires.cool --validate-index
$sing hictk validate $hicdir/hic/contact_map.q30.hic

# Metadata
$sing hictk metadata $hicdir/cool/Control.multires.cool::/resolutions/5000000
$sing hictk metadata $hicdir/cool/Control.multires.cool --recursive
$sing hictk metadata $hicdir/hic/contact_map.q30.hic

# Convert
$sing hictk convert $hicdir/hic/contact_map.q30.hic test.mcool
$sing hictk convert $hicdir/hic/contact_map.q30.hic test.5000.cool --resolutions 5khp
$sing hictk convert -t 8 $hicdir/cool/Control.multires.cool test.hic

# Load
$sing hictk load --help
$sing hictk load --format 4dn --bin-size 10kbp $hicdir/pairs/dedup.bwa.q30.pairs.gz load.10000.cool --chrom-sizes=genometable.hg38.txt
$sing hictk load --format 4dn --bin-size 10kbp $hicdir/pairs/dedup.bwa.q30.pairs.gz load.10000.hic --chrom-sizes=genometable.hg38.txt -t 8

# Merge
$sing hictk merge data/4DNFIZ1ZVXC8.mcool::/resolutions/10000 data/4DNFIZ1ZVXC8.mcool::/resolutions/10000 -o 4DNFIZ1ZVXC8.merged.10000.cool

# Zoomify (multires)
$sing hictk zoomify load.10000.cool load.mcool
$sing hictk zoomify load.10000.hic load.hic -t 8

# Balancing
for norm in ice scale vc
do
    $sing hictk balance $norm load.mcool
    $sing hictk balance $norm load.hic
done
