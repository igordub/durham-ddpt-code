####--------------------------------------------------------------------
####-----------  CONNECTIVITY  - TCL PLUG-IN FOR VMD --------------
####--------------------------------------------------------------------
#   This script has be produced by Tom Rodgers 2011 (c)2011 
#
#### VMD ---------------------------------------------------------------
#    this program requires the free molecular graphics package "vmd".
#    web site: http://www.ks.uiuc.edu/Research/vmd/
#
#### PARTITION --------------------------------------------------------- 
proc connectivity {cutoff {molid top} {plottype "Tube"} } { 
#   - called by the user.

    set all1 [atomselect $molid "protein"]
    set ca1 [atomselect $molid "name CA"]

    set coord [$all1 get {x y z}]
    set coordca [$ca1 get {x y z}]

    set atomid [$all1 get index]
    set atomidca [$ca1 get index]

    set resid [$all1 get residue]
    set residca [$ca1 get residue]

    foreach cd1 $coordca ind1 $atomidca res $residca {
	set num1 -1
	foreach cd2 $coordca ind2 $atomidca {
	    set dist($ind1,$ind2) [vecdist $cd1 $cd2]
	    if {$dist($ind1,$ind2) <= $cutoff} {
		set num1 [expr $num1 + 1] 
	    }
	}
	set num($res) [expr $num1 / 2]
    } 
    
    set fo [open "./neighbours.dat" "w"]
    set numlist {}

    foreach ind1 $atomid res $resid {
        lappend numlist $num($res)
    } 

    foreach res $residca {
	puts $fo "$res $num($res)"
    } 

   close $fo

    $all1 set beta $numlist
    mol modstyle 0 $molid $plottype
    mol modcolor 0 $molid "Beta"
    mol scaleminmax $molid 0 2.0 8.0

    return
}

proc get_total_mass {{molid top}} {
	eval "vecadd [[atomselect $molid all] get mass]"
}