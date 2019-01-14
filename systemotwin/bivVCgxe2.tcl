proc bivVC_CE { phen1 phen2 pedFile phenoFile outdirbase } {

## Calculates the variance components if the corretaion between two traits giving a genetic and an environmental components of the correlation
##
## In:       phen1          name of the trait1
##           phen2          name of the trait2  
##           phenoFile      The file with the phenotypes
##           outdirbase     The output directory   
##


    verbosity min
    file mkdir $outdirbase 
    set f [open "$outdirbase/gxeC2.out" "a"]
            
#   puts $f "$phen1 $phen2 h21 h22 rhog lrtg pvalg rhoe lrte pvale"
   
    ped load $pedFile
    pheno load $phenoFile


          model new
          define pp1 = inorm_$phen1
          define pp2 = inorm_$phen2

## Univariate models
          set outd "$outdirbase"
          outdir $outd
          trait pp1
          option StandErr 0
          polymod
          maximize 
 	  house
          option StandErr 1
	  maximize
          set h2r1 [parameter h2r =]
          set h2r1se [parameter h2r se]
          set c21 [parameter c2 =]
          set c21se [parameter c2 se]
          if {$h2r1se > 0} {
   	      set chi1 [expr $h2r1 * $h2r1 / $h2r1se / $h2r1se]        
              set test1 [catch {set ph1 [chi -number $chi1 1] } ]
              if {$test1 != 0 } {
                 set ph1 99.0
              } else {
                 set ph1 [expr $ph1 / 2]
              }
          } else {
             set ph1 99.0
          }
          if {$c21se > 0} {
              set chi1 [expr $c21 * $c21 / $c21se / $c21se]
              set test1 [catch {set pc1 [chi -number $chi1 1] } ]
              if {$test1 != 0 } {
                 set pc1 99.0
              } else {
                 set pc1 [expr $pc1 / 2]
              }  
          } else {
             set pc1 99.0
          }



          model new
          outdir $outd
          trait pp2
          option StandErr 0
          polymod
          maximize
          house
          option StandErr 1
          maximize
          set h2r2 [parameter h2r =]
          set h2r2se [parameter h2r se]
          set c22 [parameter c2 =]
          set c22se [parameter c2 se]
          if {$h2r2se > 0} {
              set chi2 [expr $h2r2 * $h2r2 / $h2r2se / $h2r2se]
              set test2 [catch {set ph2 [chi -number $chi2 1] } ]
              if {$test2 != 0 } {
                 set ph2 99.0
              } else {
                 set ph2 [expr $ph2 / 2]
              }
          } else {
             set ph2 99.0
          }
          if {$c22se > 0} {
              set chi2 [expr $c22 * $c22 / $c22se / $c22se]
              set test2 [catch {set pc2 [chi -number $chi2 1] } ]
              if {$test2 != 0 } {
                 set pc2 99.0
              } else {
                 set pc2 [expr $pc2 / 2]
              }
          } else {
             set pc2 99.0
          }
 
## Bivariate model
          model new
          trait pp1 pp2
          outdir $outd
          option StandErr 0
          polygsd
          maximize

## results
	  set m1_0 [parameter mean(pp1) =]
 	  set esd1_0 [parameter esd(pp1) =]
 	  set gsd1_0 [parameter gsd(pp1) =]
    	  set h21_0 [expr $gsd1_0 / ($gsd1_0 + $esd1_0)]
	  set m2_0 [parameter mean(pp2) =]
 	  set esd2_0 [parameter esd(pp2) =]
 	  set gsd2_0 [parameter gsd(pp2) =]
    	  set h22_0 [expr $gsd2_0 / ($gsd2_0 + $esd2_0)]
          set rhog_0 [parameter rhog =]
          set rhoe_0 [parameter rhoe =]

	  set gxe_0 [expr $gsd1_0 * $gsd1_0 + $gsd2_0 * $gsd2_0 - (2 * $rhog_0 * $gsd1_0 * $gsd2_0)]

          set loglike_0 [loglike]

####################################################################################3
## Add C2
	  linkgsd house.gz
	  maximize

          model save "$outdirbase/null0.mod"

          set esd1 [parameter esd(pp1) =]
          set gsd1 [parameter gsd(pp1) =]
          set esd2 [parameter esd(pp2) =]
          set gsd2 [parameter gsd(pp2) =]
          set rhog [parameter rhog =]
          set rhoe [parameter rhoe =]

          set c2sd1 [parameter qsd1(pp1) =]
          set c2sd2 [parameter qsd1(pp2) =]
	  set rhoc [parameter rhoq1 =]

          set h21 [expr $gsd1 / ($gsd1 + $c2sd1 + $esd1)]
          set h22 [expr $gsd2 / ($gsd2 + $c2sd2 + $esd2)]

          set gxe [expr $gsd1 * $gsd1 + $gsd2 * $gsd2 - (2 * $rhog * $gsd1 * $gsd2)]
          set cxe [expr $c2sd1 * $c2sd1 + $c2sd2 * $c2sd2 - (2 * $rhoc * $c2sd1 * $c2sd2)]
          set exe [expr $esd1 * $esd1 + $esd2 * $esd2 - (2 * $rhoe * $esd1 * $esd2)]

          set loglike1 [loglike]


####################################################################################3
## Calculate LRT test for genetic correlation

           parameter rhog = 0
           constraint rhog = 0
           maximize
           set loglike2 [loglike]

           set lrtg [expr 2 * ($loglike1 - $loglike2)]
           ## for prevent error for certain out of bound conditions
           set test [catch {set pvalg [chi -number $lrtg 1] } ]
           if {$test != 0 } {set pvalg 99.0}

#####################################################################################3
## Calculate LRT test for environmental correlation
           model load "$outdirbase/null0.mod"

           parameter rhoe = 0
           constraint rhoe = 0
           maximize
           set loglike3 [loglike]

           set lrte [expr 2 * ($loglike1 - $loglike3)]
           ## for prevent error for certain out of bound conditions
           set test [catch {set pvale [chi -number $lrte 1] } ]
           if {$test != 0 } {set pvale 99.0}

####################################################################################3
## Calculate LRT test for genetic correlation = 1
           model load "$outdirbase/null0.mod"

           parameter rhog = 1
           constraint rhog = 1
           maximize
           set loglike4 [loglike]

           set lrtg1 [expr 2 * ($loglike1 - $loglike4)]
           ## for prevent error for certain out of bound conditions
           set test [catch {set pvalg1 [chi -number $lrtg1 1] } ]
           if {$test != 0 } {set pvalg1 99.0}

####################################################################################3
## Calculate LRT test for equal variances 
           model load "$outdirbase/null0.mod"

           parameter rhog = 1
	   parameter esd(pp2) = $esd1
	   parameter gsd(pp2) = $gsd1           
 
   	   constraint gsd(pp1) = gsd(pp2)
           maximize
           model save "$outdirbase/null0equalVar.mod"

           set gsdequal [parameter gsd(pp1) =]
           set gxeEqualVars [expr 2 * $gsdequal * $gsdequal * (1 - $rhog)]

           set loglike5 [loglike]

           set lrt5 [expr 2 * ($loglike1 - $loglike5)]
           ## for prevent error for certain out of bound conditions
           set test [catch {set pval5 [chi -number $lrt5 1] } ]
           if {$test != 0 } {set pval5 99.0}

####################################################################################3
## Calculate LRT test for equal variances and genetic correlation = 1

           parameter rhog = 1
           constraint rhog = 1
           maximize

           set loglike6 [loglike]

           set lrt6 [expr 2 * ($loglike1 - $loglike6)]
           ## for prevent error for certain out of bound conditions
           set test [catch {set pval6 [chi -number $lrt6 2] } ]
           if {$test != 0 } {set pval6 99.0}

####################################################################################3
## Calculate LRT test for equal variances and environmental correlation = 1
           model load "$outdirbase/null0equalVar.mod"

           parameter rhoe = 0.99
           constraint rhoe = 0.99
        
           maximize

           set loglike7 [loglike]

           set lrt7 [expr 2 * ($loglike1 - $loglike7)]
           ## for prevent error for certain out of bound conditions
           set test [catch {set pval7 [chi -number $lrt7 2] } ]
           if {$test != 0 } {set pval7 99.0}


####################################################################################3
## Calculate LRT test for genetic variances = 0
           model load "$outdirbase/null0.mod"

           parameter rhog = 0
   	   constraint rhog = 0
           parameter esd(pp1) = 1
           parameter gsd(pp1) = 0
           parameter esd(pp2) = 1
           parameter gsd(pp2) = 0
           constraint gsd(pp1) = 0
  	   constraint gsd(pp2) = 0
        
           maximize

           set loglike8 [loglike]

           set lrt8 [expr 2 * ($loglike1 - $loglike8)]
           ## for prevent error for certain out of bound conditions
           set test [catch {set pval8 [expr [chi -number $lrt8 2] / 4 + [chi -number $lrt8 1] / 2] } ]
           if {$test != 0 } {set pval8 99.0}

###################################################################################3
## Calculate LRT test for equal c2 variances and c2 correlation = 1
           model load "$outdirbase/null0.mod"

	   parameter esd(pp2) = $esd1
           parameter gsd(pp2) = $gsd1
           parameter qsd1(pp2) = $c2sd1
           constraint qsd1(pp1) = qsd1(pp2)
         
           parameter rhoq1 = 1
           constraint rhoq1 = 1

           maximize

           set loglike9 [loglike]

           set lrt9 [expr 2 * ($loglike1 - $loglike9)]
           ## for prevent error for certain out of bound conditions
           set test [catch {set pval9 [chi -number $lrt9 2] } ]
           if {$test != 0 } {set pval9 99.0}


####################################################################################3
## Calculate LRT test for c2 variances = 0
           model load "$outdirbase/null0.mod"

           parameter rhoq1 = 0
           constraint rhoq1 = 0
           parameter esd(pp1) = 0.9
           parameter gsd(pp1) = 0.1
 	   parameter qsd1(pp1) = 0
           parameter esd(pp2) = 0.0
           parameter gsd(pp2) = 0.1
 	   parameter qsd1(pp2) = 0
           constraint qsd1(pp1) = 0
           constraint qsd1(pp2) = 0

           maximize

           set loglike10 [loglike]

           set lrt10 [expr 2 * ($loglike1 - $loglike10)]
           ## for prevent error for certain out of bound conditions
           set test [catch {set pval10 [expr [chi -number $lrt10 2] / 4 + [chi -number $lrt10 1] / 2] } ]
           if {$test != 0 } {set pval10 99.0}


#####################################################################################3
## Calculate LRT test for c2 correlation
           model load "$outdirbase/null0.mod"

           parameter rhoq1 = 0
           constraint rhoq1 = 0
           maximize
           set loglike11 [loglike]

           set lrt11 [expr 2 * ($loglike1 - $loglike11)]
           ## for prevent error for certain out of bound conditions
           set test [catch {set pval11 [chi -number $lrt11 1] } ]
           if {$test != 0 } {set pval11 99.0}


## Write results
#           puts $f "$phenoFile $h2r1 h2r1 $h21 $h22 $rhog $lrtg $pvalg $rhoe $lrte $pvale"

           puts $f "$phen1 $phen2 $h2r1 $h2r1se $c21 $c21se $ph1 $pc1 $h2r2 $h2r2se $c22 $c22se $ph2 $pc2 $m1_0 $m2_0 $h21_0 $h22_0 $gsd1_0 $esd1_0 $gsd2_0 $esd2_0 $rhog_0 $rhoe_0 $gxe_0 $loglike_0 $h21 $h22 $gsd1 $c2sd1 $esd1 $gsd2 $c2sd2 $esd2 $rhog $rhoc $rhoe $gxe $cxe $exe $loglike1 $loglike2 $pvalg $loglike3 $pvale $loglike4 $pvalg1 $gsdequal $gxeEqualVars $loglike5 $pval5 $loglike6 $pval6 $loglike7 $pval7 $loglike8 $pval8 $loglike9 $pval9 $loglike10 $pval10 $loglike11 $pval11"

     close $f
}





