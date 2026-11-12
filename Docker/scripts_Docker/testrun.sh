#!/bin/bash
#
# Smoke test for the tools installed in the CustardPy Docker/Apptainer image.
#
# It invokes each tool with a harmless flag (--help/-h, or no arguments for
# scripts that print usage on missing arguments) and checks that the process
# actually starts and runs as the installed tool -- not that it succeeds.
# A script printing its own usage and exiting non-zero is a PASS; a missing
# binary, a broken shared-library link, or a Python import error is a FAIL.
#
# Usage:
#   testrun.sh [-v]
#     -v : also print the captured output of each check
#
set -u

verbose=0
if [ "${1:-}" = "-v" ]; then
    verbose=1
fi

npass=0
nfail=0
nskip=0
failed_checks=()

workdir=$(mktemp -d)
trap 'rm -rf "$workdir"' EXIT
cd "$workdir"

section() {
    echo
    echo "== $1 =="
}

# check <label> <command...>
# Runs <command...>, treats "not found / not executable" exit codes and
# well-known failure signatures (missing shared lib, Python import error,
# uncaught traceback) as FAIL. Any other outcome (including a usage message
# and a non-zero exit from argument parsing) is a PASS.
check() {
    local label=$1; shift
    local out rc
    out=$("$@" 2>&1)
    rc=$?

    if [ $rc -eq 127 ] || [ $rc -eq 126 ]; then
        echo "[FAIL] $label (exit $rc)"
        failed_checks+=("$label")
        nfail=$((nfail + 1))
    elif echo "$out" | grep -qE \
        'error while loading shared libraries|cannot execute binary file|command not found|ModuleNotFoundError|ImportError:|Traceback \(most recent call last\)|there is no package called|Execution halted'; then
        echo "[FAIL] $label"
        failed_checks+=("$label")
        nfail=$((nfail + 1))
    else
        echo "[ OK ] $label"
        npass=$((npass + 1))
    fi

    if [ $verbose -eq 1 ]; then
        echo "$out" | sed 's/^/       /'
    fi
}

# check_env <label> <micromamba-env> <command...>
check_env() {
    local label=$1; local envname=$2; shift 2
    check "$label" run_env.sh "$envname" "$@"
}

# check_cmd <label> <command>: PASS if the command is on PATH and executable.
# Used for scripts that do not handle --help and fail on dummy arguments.
check_cmd() {
    local p
    p=$(command -v "$2")
    if [ -n "$p" ] && [ -x "$p" ]; then
        echo "[ OK ] $1"
        npass=$((npass + 1))
    else
        echo "[FAIL] $1 (not found on PATH)"
        failed_checks+=("$1")
        nfail=$((nfail + 1))
    fi
}

# check_path <label> <path>: PASS if the file or directory exists
check_path() {
    if [ -e "$2" ]; then
        echo "[ OK ] $1"
        npass=$((npass + 1))
    else
        echo "[FAIL] $1 ($2 not found)"
        failed_checks+=("$1")
        nfail=$((nfail + 1))
    fi
}

skip() {
    echo "[SKIP] $1 ($2)"
    nskip=$((nskip + 1))
}

# ---------------------------------------------------------------------------
section "Mapping tools"
check "bwa"             bwa
check "bowtie"          bowtie --version
check "bowtie2"         bowtie2 --version
check "chromap"         chromap -h
check "samtools"        samtools --version
check "bedtools"        bedtools --version
check "HiCUP"           hicup --help

# ---------------------------------------------------------------------------
section "Juicer pipeline"
check "java"                       java -version
check "juicertools.sh"             juicertools.sh
check "custardpy_juicer"           custardpy_juicer
check "juicer_map.sh"              juicer_map.sh
check "juicer_pigz.sh"             juicer_pigz.sh nonexistent_dir
check "juicer_unpigz.sh"           juicer_unpigz.sh nonexistent_dir
check_cmd "plot_distance_count.sh" plot_distance_count.sh
check "juicer_callTAD.sh"          juicer_callTAD.sh
check "makeMatrix_intra.sh"        makeMatrix_intra.sh
check "makeEigen.sh"               makeEigen.sh
check "makeInsulationScore.sh"     makeInsulationScore.sh
check "call_MotifFinder.sh"        call_MotifFinder.sh
check "Juicerstats.sh"             Juicerstats.sh
check "custardpy_process_hic"      custardpy_process_hic -h
check "Juicer_statistics"          Juicer_statistics
check "Juicer_remove_duplicate"    Juicer_remove_duplicate
check "Juicer_fragment"            Juicer_fragment
check "Juicer_chimeric_blacklist"  Juicer_chimeric_blacklist
check "distance_vs_count.Juicer"     distance_vs_count.Juicer
check "distance_vs_count.Juicer.log" distance_vs_count.Juicer.log

if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi >/dev/null 2>&1; then
    check "call_HiCCUPS.sh" call_HiCCUPS.sh
else
    skip "call_HiCCUPS.sh" "no GPU detected; supply --gpus all / --nv to test"
fi

# ---------------------------------------------------------------------------
section "Cooler/cooltools pipeline"
check "cooler"                  cooler --version
check "cooltools"               cooltools --version
check "pairtools"                pairtools --version
check "hic2cool"                 hic2cool -h
check "coolpup.py"               coolpup.py --help
check "hictk"                    hictk --version
check "custardpy_cooler"         custardpy_cooler -h
check "custardpy_process_cool"   custardpy_process_cool -h
check "run-cool2multirescool.sh" run-cool2multirescool.sh
check "run_fithic.sh"            run_fithic.sh
check "custardpy_phic"           custardpy_phic -h
check "custardpy_pastis"         custardpy_pastis
check "makeEigen_cool.sh"        makeEigen_cool.sh
check "Coolerstats.sh"           Coolerstats.sh

# ---------------------------------------------------------------------------
section "custardpy (PyPI) scripts"
check "checkHiCfile.py"              checkHiCfile.py --help
check "convert_JuicerDump_to_dense.py" convert_JuicerDump_to_dense.py --help
check "custardpy_clustering_boundary" custardpy_clustering_boundary --help
check "custardpy_differential_DRF"    custardpy_differential_DRF --help
check "DEG_boundary_analysis"         DEG_boundary_analysis --help
check "drawSquareMulti"               drawSquareMulti --help
check "drawSquarePair"                drawSquarePair --help
check "drawSquareRatioMulti"          drawSquareRatioMulti --help
check "drawSquareRatioPair"           drawSquareRatioPair --help
check "drawTriangleMulti"             drawTriangleMulti --help
check "drawTrianglePair"              drawTrianglePair --help
check "getBoundaryfromInsulationScore" getBoundaryfromInsulationScore --help
check "plotCompartmentGenome"         plotCompartmentGenome --help
check "plotInsulationScore"           plotInsulationScore --help
check "plotMultiScaleInsulationScore" plotMultiScaleInsulationScore --help
check_cmd "plotHiCMatrix"             plotHiCMatrix
check "plotHiCfeature"                plotHiCfeature --help

# ---------------------------------------------------------------------------
section "Utility tools (base environment)"
check "calculate_compartment_strength" calculate_compartment_strength
check "run_3DChromatin_ReplicateQC.sh" run_3DChromatin_ReplicateQC.sh -h
check_cmd "visualize_QC.py"            visualize_QC.py
check_cmd "generate_binlist_from_gtfile.py" generate_binlist_from_gtfile.py
check "get_qc.py"                      get_qc.py --help
check "OnTAD"                          OnTAD
check "bedToBigBed"                    bedToBigBed
check "chess"                          chess --help
check "mustache"                       mustache --help
check "chromosight"                    chromosight --help
check "HOMER (homerTools)"             homerTools
check "h1d"                            h1d --help
check "fithic"                         fithic --help
check "HiCLift"                        python -c 'import importlib.metadata as m; print(m.version("HiCLift"))'

# ---------------------------------------------------------------------------
section "Virtual environments (micromamba)"
check_env "pastis (pastis-pm2)"     pastis pastis-pm2 --help
check_env "hicexplorer (hicInfo)"   hicexplorer hicInfo --version
check_env "stripenn"                stripenn stripenn --help
check "multiqc"                     multiqc --version
check "hicpro"                      hicpro --help
check_env "macs2"                   hic-pro macs2 --version
check "run_fithichip.sh"            run_fithichip.sh FitHiChIP_HiCPro.sh -h
check_path "ChIA-PIPE (source)"     /opt/ChIA-PIPE
check_env "ChIA-PIPE env (python2 modules)" chiapipe python -c 'import numpy, Bio, pysam, regex'
check_env "ChIA-PIPE env (java)"    chiapipe java -version

# ---------------------------------------------------------------------------
section "R packages"
check "R"               Rscript --version
check "R: GENOVA"        Rscript -e 'library(GENOVA)'
check "R: FIREcaller"    Rscript -e 'library(FIREcaller)'
check "R: CALDER2"       Rscript -e 'library(CALDER)'
check "R: mango"         Rscript -e 'library(mango)'
check "R: hicrep"        Rscript -e 'library(hicrep)'
check "R: strawr"        Rscript -e 'library(strawr)'
check "R: edgeR"         Rscript -e 'library(edgeR)'
check "R: GenomicRanges" Rscript -e 'library(GenomicRanges)'

# ---------------------------------------------------------------------------
section "Data files"
check_path "Cooler restriction sites"  /opt/restriction_sites/Cooler-restriction_sites
check_path "HiC-Pro restriction sites" /opt/restriction_sites/HiCPro-restriction_sites
check_path "Juicer restriction sites"  /opt/restriction_sites/Juicer-restriction_sites
check_path "motif files"               /opt/motiffiles
check_path "juicer_tools 1.22.01"      /opt/juicer/juicer_tools.1.22.01.jar
check_path "juicer_tools 1.9.9 (jcuda)" /opt/juicer/juicer_tools.1.9.9_jcuda.0.8.jar
check_path "juicer_tools 2.20.00"      /opt/juicer/juicer_tools.2.20.00.jar
check_path "Juicebox.jar"              /opt/Juicebox.jar

# ---------------------------------------------------------------------------
echo
echo "=================================================="
echo "PASS: $npass   FAIL: $nfail   SKIP: $nskip"
if [ $nfail -gt 0 ]; then
    echo "Failed checks:"
    for c in "${failed_checks[@]}"; do
        echo "  - $c"
    done
fi
echo "=================================================="

exit $nfail
