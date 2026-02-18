import delimited "/Users/nhacker/Library/CloudStorage/OneDrive-TexasA&MUniversity/TAMU/VERO_Research/Adaptive_Sampling_M_Haemolytica/Results/Final/All_Mh_base_stats.csv", clear

gen pool_str = regexs(1) if regexm(sample, "Pool([0-9]+)")


gen time_str = regexs(1) if regexm(sample,"Pool[0-9]+_(0-[0-9]+)_barcode")
gen time_hr = real(substr(time_str, 3, .))


gen barcode_str = regexs(1) if regexm(sample, "barcode([0-9]+)")
egen sample_id = group(pool_str barcode_str)


generate conc = "High" 
replace conc = "Low" if inlist(pool_str, "4", "5", "6")
replace conc = "Neg" if inlist(pool_str, "7", "8", "9")


histogram mapped_bases

histogram proportion_mapped, by(conc)

twoway ///
    (lowess mapped_bases time_hr if inlist(real(pool_str),1,2,3), ///
        lcolor(blue) lwidth(medium) legend(label(1 "High Mh"))) ///
    (lowess mapped_bases time_hr if inlist(real(pool_str),4,5,6), ///
        lcolor(red) lwidth(medium) legend(label(2 "Low Mh"))) ///
    (lowess mapped_bases time_hr if inlist(real(pool_str),7,8,9), ///
        lcolor(green) lwidth(medium) legend(label(3 "Neg Mh"))), ///
    by(pool_str, note("") title("") legend(pos(6))) ///
    xtitle("Sequencing time") ///
    ytitle("Bases mapping to {it: Mannheimia haemolytica}")

lowess mapped_bases time_hr, by(pool_str, note("")legend(at(4)) title("") ) xtitle("Sequencing time") mcolor(#4B0082) ytitle("Reads mapping to{it: Mannheimia haemolytica}") legend(order(1 2) label(1 "Mapped bases") label(2 "Lowess smoother"))

lowess total_bases time_hr, by(pool_str, note("")legend(at(4)) title("") ) xtitle("Sequencing time")  ytitle("Total bases") legend(order(1 2) label(1 "Total reads") label(2 "Lowess smoother"))


sort sample_id time_hr

twoway line mapped_bases time_hr if pool_str=="1", by(sample_id)

twoway line mapped_bases time_hr if pool_str=="2", by(sample_id)

twoway line mapped_bases time_hr if pool_str=="3", by(sample_id)

twoway line mapped_bases time_hr if pool_str=="4", by(sample_id)

twoway line mapped_bases time_hr if pool_str=="5", by(sample_id)

twoway line mapped_bases time_hr if pool_str=="6", by(sample_id)

twoway line mapped_bases time_hr if pool_str=="7", by(sample_id)

twoway line mapped_bases time_hr if pool_str=="8", by(sample_id)

twoway line mapped_bases time_hr if pool_str=="9", by(sample_id)


xtset sample_id time_hr, delta(6)

encode(conc), generate(conc_cat)



* xtgee mapped i.time_hr i.conc_cat, family(binomial total_reads) corr(unstructured) vce(robust) /* doesn't converge*/

* xtgee mapped i.time_hr i.conc_cat, family(binomial total_reads) corr(ar 1) vce(robust) /* doesn't converge*/


* xtgee mapped i.time_hr i.conc_cat, family(binomial total_reads) corr(independent) vce(robust) /*main effects model*/


xtgee mapped_bases i.time_hr##i.conc_cat, family(binomial total_bases) corr(independent) vce(robust) 

* xtgee mapped i.time_hr##i.conc_cat if conc!="Neg", family(binomial total_reads) corr(independent) vce(robust) /* does not include the Neg group*/


testparm i.time_hr

testparm i.conc_cat

testparm i.time_hr#i.conc_cat /* we cant use with main effects model*/

testparm  i.conc_cat i.time_hr#i.conc_cat /*cant use with main effects model*/

margins time_hr#conc_cat, expression(exp(predict(xb))/(1+exp(predict(xb))))
marginsplot, yscale(range(0 1)) ytitle("Proportion of mapped bases") ytick(0 0.2 0.4 0.6 0.8 1) ylabel(0 "0" 0.2 "0.2" 0.4 "0.4" 0.6 "0.6" 0.8 "0.8" 1 "1") title("")

contrast ar.time_hr@i.conc_cat, mcompare(sidak)



* xtgee mapped_reads i.time_hr i.conc_cat, corr(unstructured) vce(robust) /* doesn't converge*/

* xtgee mapped i.time_hr i.conc_cat,  corr(ar 1) vce(robust) /* doesn't converge*/


* xtgee mapped i.time_hr i.conc_cat, corr(independent) vce(robust) /* Simple main effects model */

xtgee mapped i.time_hr##i.conc_cat, corr(independent) vce(robust) /* overparameterized model*/


* xtgee mapped i.time_h##i.conc_cat if conc!="Neg", corr(independent) vce(robust)


testparm i.time_hr

testparm i.conc_cat

testparm i.time_hr#i.conc_cat

testparm i.conc_cat i.time_hr#i.conc_cat

margins i.time_hr#conc_cat
marginsplot, by(conc_cat) ///
 ytitle("Bases mapping to{it: Mannheimia haemolytica}") ///
 xtitle("Sequencing time") ///
 plotopts(mcolor(#blue) lcolor(#blue)) ///
 recastci(rarea) ciopts(color(blue%10))
	
contrast ar.time_hr@i.conc_cat, mcompare(sidak)






/* Creation of incident reads (mapped and total): new reads from a specific time span*/

sort sample_id time_hr

by sample_id: generate new_mapped=mapped[_n] - mapped[_n-1]


by sample_id: replace new_mapped=mapped if time_hr==6



by sample_id: generate new_total=total_bases[_n] - total_bases[_n-1]


by sample_id: replace new_total=total_bases if time_hr==6


*/Generalized Estimating Equations: logistic regression on proportion of new reads (incident reads) */


* xtgee new_mapped i.time_hr i.conc_c, family(binomial new_total) corr(independent) vce(robust) /* zeros in the new_total variable*/




*/Generalized Estimating Equations: linear regression on proportion of new reads (incident mapped and total reads) */

xtgee new_mapped i.time_hr##i.conc_cat, corr(independent) vce(robust)

testparm i.time_hr

testparm i.conc_cat

testparm i.conc_cat#i.time_hr

margins i.time_hr#i.conc_cat

marginsplot, by(conc_cat) title("") ytitle("New mapped bases") plotopts(mcolor(black) lcolor(black))recastci(rarea) ciopts(color(blue%10))

contrast ar.time_hr@conc_cat, mcompare(sidak)
contrast r.conc_cat@time_hr, mcompare(sidak)




* =========================================
*/ 16S Models
* =========================================

import delimited "/Users/nhacker/Library/CloudStorage/OneDrive-TexasA&MUniversity/TAMU/VERO_Research/Adaptive_Sampling_M_Haemolytica/Results/Final/All_16S_base_stats.csv", clear

gen pool_str = regexs(1) if regexm(sample, "Pool([0-9]+)")


gen time_str = regexs(1) if regexm(sample,"Pool[0-9]+_(0-[0-9]+)_barcode")
gen time_hr = real(substr(time_str, 3, .))


gen barcode_str = regexs(1) if regexm(sample, "barcode([0-9]+)")
egen sample_id = group(pool_str barcode_str)


generate conc = "High" 
replace conc = "Low" if inlist(pool_str, "4", "5", "6")
replace conc = "Neg" if inlist(pool_str, "7", "8", "9")


histogram mapped_bases


histogram proportion_mapped, by(conc)

twoway ///
    (lowess mapped_bases time_hr if inlist(real(pool_str),1,2,3), ///
        lcolor(blue) lwidth(medium) legend(label(1 "High Mh"))) ///
    (lowess mapped_bases time_hr if inlist(real(pool_str),4,5,6), ///
        lcolor(red) lwidth(medium) legend(label(2 "Low Mh"))) ///
    (lowess mapped_bases time_hr if inlist(real(pool_str),7,8,9), ///
        lcolor(green) lwidth(medium) legend(label(3 "Neg Mh"))), ///
    by(pool_str, note("") title("") legend(pos(6))) ///
    xtitle("Sequencing time") ///
    ytitle("Bases mapping to 16S Region")

lowess mapped_bases time_hr, by(pool_str, note("")legend(at(4)) title("") ) xtitle("Sequencing time") mcolor(#4B0082) ytitle("Bases mapping to 16S Region") legend(order(1 2) label(1 "Mapped bases") label(2 "Lowess smoother"))


lowess total_bases time_hr, by(pool_str, note("")legend(at(4)) title("") ) xtitle("Sequencing time")  ytitle("Total bases") legend(order(1 2) label(1 "Total bases") label(2 "Lowess smoother"))


sort sample_id time_hr

twoway line mapped_bases time_hr if pool_str=="1", by(sample_id)

twoway line mapped_bases time_hr if pool_str=="2", by(sample_id)

twoway line mapped_bases time_hr if pool_str=="3", by(sample_id)

twoway line mapped_bases time_hr if pool_str=="4", by(sample_id)

twoway line mapped_bases time_hr if pool_str=="5", by(sample_id)

twoway line mapped_bases time_hr if pool_str=="6", by(sample_id)

twoway line mapped_bases time_hr if pool_str=="7", by(sample_id)

twoway line mapped_bases time_hr if pool_str=="8", by(sample_id)

twoway line mapped_bases time_hr if pool_str=="9", by(sample_id)


xtset sample_id time_hr, delta(6)

encode(conc), generate(conc_cat)



*xtgee mapped i.time_hr i.conc_cat, family(binomial total_bases) corr(unstructured) vce(robust) /* doesn't converge*/

*xtgee mapped i.time_hr i.conc_cat, family(binomial total_bases) corr(ar 1) vce(robust) /* doesn't converge*/


* xtgee mapped i.time_hr i.conc_cat, family(binomial total_bases) corr(independent) vce(robust) /*main effects model*/


xtgee mapped i.time_hr##i.conc_cat, family(binomial total_bases) corr(independent) vce(robust) 

* xtgee mapped i.time_hr##i.conc_cat if conc!="Neg", family(binomial total_bases) corr(independent) vce(robust) /* does not include the Neg group*/


testparm i.time_hr

testparm i.conc_cat

testparm i.time_hr#i.conc_cat /* we cant use with main effects model*/

testparm  i.conc_cat i.time_hr#i.conc_cat /*cant use with main effects model*/

margins time_hr#conc_cat, expression(exp(predict(xb))/(1+exp(predict(xb))))
marginsplot, yscale(range(0 1)) ytitle("Proportion of 16S mapped bases") ytick(0 0.2 0.4 0.6 0.8 1) ylabel(0 "0" 0.2 "0.2" 0.4 "0.4" 0.6 "0.6" 0.8 "0.8" 1 "1")

contrast ar.time_hr@i.conc_cat, mcompare(sidak)



* xtgee mapped_reads i.time_hr i.conc_cat, corr(unstructured) vce(robust) /* doesn't converge*/

* xtgee mapped i.time_hr i.conc_cat,  corr(ar 1) vce(robust) /* doesn't converge*/


* xtgee mapped i.time_hr i.conc_cat, corr(independent) vce(robust) /* Simple main effects model */

xtgee mapped i.time_hr##i.conc_cat, corr(independent) vce(robust) /* overparameterized model*/


* xtgee mapped i.time_h##i.conc_cat if conc!="Neg", corr(independent) vce(robust)


testparm i.time_hr

testparm i.conc_cat

testparm i.time_hr#i.conc_cat

testparm i.conc_cat i.time_hr#i.conc_cat

margins i.time_hr#conc_cat
marginsplot, by(conc_cat) ///
 ytitle("Bases mapping to 16S Region") ///
 xtitle("Sequencing time") ///
 plotopts(mcolor(#4B0082) lcolor(#4B0082)) ///
 recastci(rarea) ciopts(color(green%10))

contrast ar.time_hr@i.conc_cat, mcompare(sidak)






/* Creation of incident reads (mapped and total): new reads from a specific time span*/

sort sample_id time_hr

by sample_id: generate new_mapped=mapped[_n] - mapped[_n-1]


by sample_id: replace new_mapped=mapped if time_hr==6



by sample_id: generate new_total=total_bases[_n] - total_bases[_n-1]


by sample_id: replace new_total=total_bases if time_hr==6


*/Generalized Estimating Equations: logistic regression on proportion of new reads (incident reads) */


* xtgee new_mapped i.time_hr i.conc_c, family(binomial new_total) corr(independent) vce(robust) /* zeros in the new_total variable*/




*/Generalized Estimating Equations: linear regression on proportion of new reads (incident mapped and total reads) */

xtgee new_mapped i.time_hr##i.conc_cat, corr(independent) vce(robust)

testparm i.time_hr

testparm i.conc_cat

testparm i.conc_cat#i.time_hr

margins i.time_hr#i.conc_cat

marginsplot, by(conc_cat) title("") ytitle("Total new 16S mapped bases") plotopts(mcolor(black) lcolor(black))recastci(rarea) ciopts(color(green%10))

contrast ar.time_hr@conc_cat, mcompare(sidak)
contrast r.conc_cat@time_hr, mcompare(sidak)


























