# CRYPTAD — Ligand COM drift analysis for all complex MD replicas
#
# Usage:
#   vmd -dispdev none -e 09_scripts/09_complex_md/05_drift_analysis.tcl
#
# Override paths via environment variables:
#   CRYPTAD_ROOT  — project root (default: two levels above this script)
#   CRYPTAD_OUT   — output file path (default: <root>/04_virtual_screening/complex_md_results/drift_results.txt)

# ── Resolve project root ──────────────────────────────────────────────────────
# Default: two directories above this script's location
# (09_scripts/09_complex_md/ → 09_scripts/ → project root)
set script_dir [file dirname [file normalize [info script]]]
set default_root [file normalize "$script_dir/../.."]
set base_root [expr {[info exists env(CRYPTAD_ROOT)] ? $env(CRYPTAD_ROOT) : $default_root}]

set complexmd_base "$base_root/02_md_simulations/complex_md"

# ── Output file ───────────────────────────────────────────────────────────────
set default_out "$base_root/04_virtual_screening/complex_md_results/drift_results.txt"
set outpath [expr {[info exists env(CRYPTAD_OUT)] ? $env(CRYPTAD_OUT) : $default_out}]
file mkdir [file dirname $outpath]

set outfile [open $outpath w]

puts $outfile "# CRYPTAD ligand COM drift analysis"
puts $outfile "# Drift = displacement of ligand COM from frame-0 position (Angstrom)"
puts $outfile "# Columns: system  rep  mean_drift  max_drift  final_drift  n_frames"
puts $outfile [string repeat "-" 80]

# ── Auto-discover systems and replicas ───────────────────────────────────────
set systems [lsort [glob -nocomplain -tails -directory $complexmd_base -type d *]]

if {[llength $systems] == 0} {
    puts stderr "ERROR: no system directories found under $complexmd_base"
    close $outfile
    exit 1
}

foreach sys $systems {
    set sys_dir "$complexmd_base/$sys"
    set reps [lsort [glob -nocomplain -tails -directory $sys_dir -type d rep*]]

    foreach rep $reps {
        set dir "$sys_dir/$rep"
        set ref "$dir/ref_nw.pdb"
        set xtc "$dir/complexmd_fit_nw.xtc"

        if {![file exists $ref] || ![file exists $xtc]} {
            puts "SKIP $sys / $rep — ref_nw.pdb or complexmd_fit_nw.xtc missing"
            continue
        }

        puts "Processing $sys / $rep ..."

        mol new $ref type pdb waitfor all
        mol addfile $xtc type xtc waitfor all

        set lig [atomselect top "resname MOL"]
        set nframes [molinfo top get numframes]

        # Frame 0 reference COM
        $lig frame 0
        $lig update
        set com0 [measure center $lig weight mass]

        set total_drift 0.0
        set max_drift   0.0
        set final_drift 0.0

        for {set f 0} {$f < $nframes} {incr f} {
            $lig frame $f
            $lig update
            set com [measure center $lig weight mass]
            set dx [expr {[lindex $com 0] - [lindex $com0 0]}]
            set dy [expr {[lindex $com 1] - [lindex $com0 1]}]
            set dz [expr {[lindex $com 2] - [lindex $com0 2]}]
            set d  [expr {sqrt($dx*$dx + $dy*$dy + $dz*$dz)}]
            set total_drift [expr {$total_drift + $d}]
            if {$d > $max_drift} { set max_drift $d }
            if {$f == $nframes - 1} { set final_drift $d }
        }

        set mean_drift [expr {$total_drift / $nframes}]

        set line [format "%-35s  %-4s  %7.2f  %7.2f  %7.2f  %d" \
            $sys $rep $mean_drift $max_drift $final_drift $nframes]
        puts $outfile $line
        puts $line

        $lig delete
        mol delete top
    }
}

close $outfile
puts "Done. Results written to $outpath"
