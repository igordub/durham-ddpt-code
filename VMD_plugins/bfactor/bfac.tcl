
proc bfactor { {scale 1} {molid top} {plottype "NewCartoon"} } { 

    set all1 [atomselect $molid "protein"]
    set atomid [$all1 get index]
    set resid [$all1 get residue]
    set bfac [$all1 get beta]

    foreach b $bfac res $resid {
	set num($res) [expr $b * $scale ] 
    }

    set numlist {}

    foreach ind1 $atomid res $resid {
        lappend numlist $num($res)
    } 

    color scale midpoint 0.00
#    color scale offset 0.00 
    color scale method BWR
    
    $all1 set beta $numlist
    mol modstyle 0 $molid $plottype
    mol modcolor 0 $molid "Beta"
    mol scaleminmax $molid 0 10.0 100.0

    return
}
