set out [open "SPpnpg328_.dat" w]

set scale 180.0
set sysref [mol new em_pBC.pdb]
set sys [mol new em_pBC.pdb]
mol addfile pnpg328_1us_pro-ligski1000pbC.xtc first 0 last -1 step 1 waitfor all


set nf [molinfo $sys get numframes]
puts $nf

for {set f 0} {$f <= $nf} {incr f} {
	molinfo $sys set frame $f
#*	set t [expr {($f*1.0)/1000.0}]
	set t [expr {($f*1.0)/1.0}]
	puts "Frame: $t"
		
	set sum 0.0
        set fres_b1 {22 316 243 167 166 356 357 224 225 223 169 296 298 174 181 409 410 402 122 121 326 300 328 412 170 123 22 411 299}
	foreach r $fres_b1 { 
		set a [atomselect $sysref "resid $r and name CA"] 
	        set n [$a get resname]
		 if {$n != "GLY"} {
	        	set ref_phi1 [$a get phi]
        		set ref_psi1 [$a get psi]
		
                	set sel_b1 [atomselect $sys "resid $r and name CA"]
		        set curr_phi1 [$sel_b1 get phi]
        	        set curr_psi1 [$sel_b1 get psi]
                	set delphi1 [expr {$curr_phi1 - $ref_phi1}]
                set delpsi1 [expr {$curr_psi1 - $ref_psi1}]
                set mphi1 [expr {$delphi1/$scale}]
                set mpsi1 [expr {$delpsi1/$scale}]
                set np1 [expr {abs($mphi1)}]
                set ns1 [expr {abs($mpsi1)}]
                set sum [expr {$sum + ((exp(-1*$np1))*(exp(-1*$ns1)))}]
		$sel_b1 delete
		unset ref_phi1 ref_psi1 curr_phi1 curr_psi1 delphi1 delpsi1 mphi1 mpsi1 np1 ns1 
	}
}
        set pval [expr {$sum/(1.0*29)}] 
        puts $out "$t $pval"
	#unset sum_b1


}
close $out
exit
