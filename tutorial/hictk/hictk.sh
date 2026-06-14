sing="apptainer exec --bind /work,/work2,/work3 /work/SingularityImages/custardpy.3.5.0.sif"

$sing hictk --version

hicdir=../01.Cooler_Hi-C/CustardPyResults/Cooler_hg38/Control

# Dump
$sing hictk dump $hicdir/cool/contact.bwa.q30.mcool
$sing hictk dump $hicdir/hic/contact.bwa.q30.hic --table normalizations
$sing hictk dump $hicdir/hic/contact.bwa.q30.hic --table resolutions

$sing hictk dump $hicdir/cool/contact.bwa.q30.5000.cool --join
$sing hictk dump $hicdir/cool/contact.bwa.q30.5000.cool --cis-only --join
$sing hictk dump $hicdir/cool/contact.bwa.q30.5000.cool --trans-only --join

# Validate
$sing hictk validate $hicdir/cool/contact.bwa.q30.mcool --validate-index
$sing hictk validate $hicdir/hic/contact.bwa.q30.hic

# Metadata
$sing hictk metadata $hicdir/cool/contact.bwa.q30.mcool::/resolutions/5000000
$sing hictk metadata $hicdir/cool/contact.bwa.q30.mcool --recursive
$sing hictk metadata $hicdir/hic/contact.bwa.q30.hic

# Convert
$sing hictk convert $hicdir/hic/contact.bwa.q30.hic contact.hicto.mcool
$sing hictk convert $hicdir/hic/contact.bwa.q30.hic contact.hicto.5000.cool --resolutions 5kbp
$sing hictk convert -t 8 $hicdir/cool/contact.bwa.q30.mcool contact.mcoolto.hic

# Load
$sing hictk load --help
$sing hictk load --format 4dn --bin-size 5kbp $hicdir/pairs/dedup.bwa.q30.pairs.gz dedup.bwa.pairto.5000.cool --chrom-sizes=genometable.hg38.txt
$sing hictk load --format 4dn --bin-size 5kbp $hicdir/pairs/dedup.bwa.q30.pairs.gz load.5000.pairto.hic --chrom-sizes=genometable.hg38.txt -t 8

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
