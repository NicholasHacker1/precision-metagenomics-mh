import delimited "/Users/nhacker/Library/CloudStorage/OneDrive-TexasA&MUniversity/TAMU/VERO_Research/Adaptive_Sampling_M_Haemolytica/Results/Final/All_Mh_stats.csv", clear

gen pool_str = regexs(1) if regexm(sample, "Pool([0-9]+)")


gen time_str = regexs(1) if regexm(sample,"Pool[0-9]+_(0-[0-9]+)_barcode")
gen time_hr = real(substr(time_str, 3, .))


gen barcode_str = regexs(1) if regexm(sample, "barcode([0-9]+)")
egen sample_id = group(pool_str barcode_str)


generate conc = "High" 
replace conc = "Low" if inlist(pool_str, "4", "5", "6")
replace conc = "Neg" if inlist(pool_str, "7", "8", "9")

tab conc
tab pool_str

*Total Reads by concentration and overall
histogram total_reads if time_hr == 72, /// 
width(100) ///
 start(0) ///
 ylabel(0(5)30) ///
 xlabel(0(1000)12000) ///
 frequency scheme(plotplain) ///
 color(maroon) ///
 lcolor(maroon) ///
 xtitle("Total Reads") ///
 ytitle("Sample Frequency") ///
 title("Frequnecy of Total Reads")
 
histogram total_reads if time_hr == 72, by(conc) /// 
width(100) ///
 start(0) ///
 xlabel(0(2000)12000) ///
 frequency scheme(plotplain) ///
 color(maroon) ///
 lcolor(maroon) ///
 xtitle("Total Reads") ///
 ytitle("Sample Frequency") ///
 title("Frequnecy of Total Reads")

*Mapped reads by concentration and overall
histogram mapped_reads if time_hr == 72, by(conc) /// 
width(100) ///
 start(0) ///
 xlabel(0(2000)12000) ///
 frequency scheme(plotplain) ///
 color(maroon) ///
 lcolor(maroon) ///
 xtitle("Total Mapped Reads") ///
 ytitle("Sample Frequency") ///
 title("Frequnecy of Total Mapped Reads")
 
histogram mapped_reads if time_hr == 72, /// 
width(100) ///
 start(0) ///
 xlabel(0(1000)12000) ///
 frequency scheme(plotplain) ///
 color(maroon) ///
 lcolor(maroon) ///
 xtitle("Total Mapped Reads") ///
 ytitle("Sample Frequency") ///
 title("Frequnecy of Total Mapped Reads")

*histogram proportion_mapped if time_hr == 72, by(conc) width(.01) frequency

*twoway ///
    (lowess mapped_reads time_hr if inlist(real(pool_str),1,2,3), ///
        lcolor(blue) lwidth(medium) legend(label(1 "High Mh"))) ///
    (lowess mapped_reads time_hr if inlist(real(pool_str),4,5,6), ///
        lcolor(red) lwidth(medium) legend(label(2 "Low Mh"))) ///
    (lowess mapped_reads time_hr if inlist(real(pool_str),7,8,9), ///
        lcolor(green) lwidth(medium) legend(label(3 "Neg Mh"))), ///
    by(pool_str, note("") title("") legend(pos(6))) ///
    xtitle("Sequencing time") ///
    ytitle("Reads mapping to {it: Mannheimia haemolytica}")

*lowess mapped_reads time_hr, by(pool_str, note("")legend(at(4)) title("") ) xtitle("Sequencing time") mcolor(#4B0082) ytitle("Reads mapping to{it: Mannheimia haemolytica}") legend(order(1 2) label(1 "Mapped reads") label(2 "Lowess smoother"))
 
scatter mapped_reads time_hr, by(pool_str, note("")legend(at(4)) ///
title("Mapped Reads by Pool") ) ///
 xtitle("Sequencing time") ///
 ytitle("Mapped reads") ///
 xlabel(0(6)72) ///
 scheme(plotplain) ///
 msymbol(d) ///
 mcolor(maroon) ///
 mfcolor(maroon)
 
graph box mapped_reads if pool_str == "5", over(sample_id) ///
scheme(plotplain) ///
box(1, bcolor(maroon)) ///
title("Mapped Reads for Pool 5") ///
ytitle("Mapped Reads")
graph save pool5_mapped.gph, replace

preserve
keep if pool_str == "5"
gen group1 = 0
replace group1 = 1 if sample_id == 36
ttest mapped_reads, by(group1)
restore

graph box mapped_reads if pool_str == "7", over(sample_id) ///
scheme(plotplain) ///
box(1, bcolor(maroon)) ///
title("Mapped Reads for Pool 7") ///
ytitle("Mapped Reads")
graph save pool7_mapped.gph, replace

graph combine pool5_mapped.gph pool7_mapped.gph, cols(2) ///
title("Outliers in Pools 5 & 7")
 
scatter total_reads time_hr, by(pool_str, note("")legend(at(4)) ///
title("Total Reads by Concentration") ) ///
 xtitle("Sequencing time") ///
 ytitle("Total reads") ///
 scheme(plotplain) ///
 xlabel(0(6)72) ///
 msymbol(d) ///
 mcolor(maroon) ///
 mfcolor(maroon)

sort sample_id time_hr

scatter mapped_reads time_hr if pool_str=="1", by(sample_id) scheme(plotplain) mcolor(maroon) msymbol(d) ytitle("Mapped Reads") xtitle("Sequencing Time")

twoway line mapped_reads time_hr if pool_str=="2", by(sample_id)

twoway line mapped_reads time_hr if pool_str=="3", by(sample_id)

twoway line mapped_reads time_hr if pool_str=="4", by(sample_id)

twoway line mapped_reads time_hr if pool_str=="5", by(sample_id)

twoway line mapped_reads time_hr if pool_str=="6", by(sample_id)

twoway line mapped_reads time_hr if pool_str=="7", by(sample_id)

twoway line mapped_reads time_hr if pool_str=="8", by(sample_id)

twoway line mapped_reads time_hr if pool_str=="9", by(sample_id)


xtset sample_id time_hr, delta(6)

encode(conc), generate(conc_cat)



* xtgee mapped i.time_hr i.conc_cat, family(binomial total_reads) corr(unstructured) vce(robust) /* doesn't converge*/

* xtgee mapped i.time_hr i.conc_cat, family(binomial total_reads) corr(ar 1) vce(robust) /* doesn't converge*/


* xtgee mapped i.time_hr i.conc_cat, family(binomial total_reads) corr(independent) vce(robust) /*main effects model*/


xtgee mapped i.time_hr##i.conc_cat, family(binomial total_reads) corr(independent) vce(robust) 

* xtgee mapped i.time_hr##i.conc_cat if conc!="Neg", family(binomial total_reads) corr(independent) vce(robust) /* does not include the Neg group*/


testparm i.time_hr

testparm i.conc_cat

testparm i.time_hr#i.conc_cat /* we cant use with main effects model*/

testparm  i.conc_cat i.time_hr#i.conc_cat /*cant use with main effects model*/

margins time_hr#conc_cat, expression(exp(predict(xb))/(1+exp(predict(xb))))
marginsplot, yscale(range(0 1)) ytitle("Proportion of Mapped Reads") ytick(0 0.2 0.4 0.6 0.8 1) ylabel(0 "0" 0.2 "0.2" 0.4 "0.4" 0.6 "0.6" 0.8 "0.8" 1 "1") title("Proportion of {it:M. haemolytica} Reads Over Sequencing Time") xtitle("Sequencing Time")scheme(plotplain)plot1opts(mcolor(maroon%70) lcolor(maroon%70) msymbol(d))plot2opts(mcolor(green%70)lcolor(green%70) msymbol(s))plot3opts(mcolor(blue%70)lcolor(blue%70) msymbol(t))ci1opts(lcolor(maroon%70))ci2opts(lcolor(green%70))ci3opts(lcolor(blue%70))

contrast ar.time_hr@i.conc_cat, mcompare(sidak)
*The reverse adjacent contrast between 12-18 hours in the high group was significant (p = 0.0351)


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
byopts(title("Reads Mapping to {it:M. haemolytica}")) ///
 ytitle("Mapped Read Count") ///
 xtitle("Sequencing time") ///
 plot1opts(mcolor(maroon%80) lcolor(maroon%80)) ///
 recastci(rarea) ci1opts(color(maroon%10))
	
contrast ar.time_hr@i.conc_cat, mcompare(sidak)
*The only significant adjacent reverse contrast are in the high group.





/* Creation of incident reads (mapped and total): new reads from a specific time span*/

sort sample_id time_hr

by sample_id: generate new_mapped=mapped[_n] - mapped[_n-1]


by sample_id: replace new_mapped=mapped if time_hr==6



by sample_id: generate new_total=total_reads[_n] - total_reads[_n-1]


by sample_id: replace new_total=total_reads if time_hr==6


*/Generalized Estimating Equations: logistic regression on proportion of new reads (incident reads) */


* xtgee new_mapped i.time_hr i.conc_c, family(binomial new_total) corr(independent) vce(robust) /* zeros in the new_total variable*/




*/Generalized Estimating Equations: linear regression on proportion of new reads (incident mapped and total reads) */

xtgee new_mapped i.time_hr##i.conc_cat, corr(independent) vce(robust)

testparm i.time_hr

testparm i.conc_cat

testparm i.conc_cat#i.time_hr

margins i.time_hr#i.conc_cat

marginsplot, by(conc_cat) byopts(title("Non-Cumulative Reads per Timepoint")) ytitle("Mapped Read Count")  xtitle("Sequencing time") ///
 plot1opts(mcolor(maroon%80) lcolor(maroon%80)) ///
 recastci(rarea) ci1opts(color(maroon%10)) ///
 yscale(range(0(100)500))

contrast ar.time_hr@conc_cat, mcompare(sidak)
contrast r.conc_cat@time_hr, mcompare(sidak)




* =========================================
*/ 16S Models
* =========================================

import delimited "/Users/nhacker/Library/CloudStorage/OneDrive-TexasA&MUniversity/TAMU/VERO_Research/Adaptive_Sampling_M_Haemolytica/Results/Final/All_16S_stats.csv", clear

gen pool_str = regexs(1) if regexm(sample, "Pool([0-9]+)")


gen time_str = regexs(1) if regexm(sample,"Pool[0-9]+_(0-[0-9]+)_barcode")
gen time_hr = real(substr(time_str, 3, .))


gen barcode_str = regexs(1) if regexm(sample, "barcode([0-9]+)")
egen sample_id = group(pool_str barcode_str)


generate conc = "High" 
replace conc = "Low" if inlist(pool_str, "4", "5", "6")
replace conc = "Neg" if inlist(pool_str, "7", "8", "9")


histogram mapped_reads


histogram proportion_mapped, by(conc)

twoway ///
    (lowess mapped_reads time_hr if inlist(real(pool_str),1,2,3), ///
        lcolor(blue) lwidth(medium) legend(label(1 "High Mh"))) ///
    (lowess mapped_reads time_hr if inlist(real(pool_str),4,5,6), ///
        lcolor(red) lwidth(medium) legend(label(2 "Low Mh"))) ///
    (lowess mapped_reads time_hr if inlist(real(pool_str),7,8,9), ///
        lcolor(green) lwidth(medium) legend(label(3 "Neg Mh"))), ///
    by(pool_str, note("") title("") legend(pos(6))) ///
    xtitle("Sequencing time") ///
    ytitle("Reads mapping to 16S Region")

lowess mapped_reads time_hr, by(pool_str, note("")legend(at(4)) title("") ) xtitle("Sequencing time") mcolor(#4B0082) ytitle("Reads mapping to 16S Region") legend(order(1 2) label(1 "Mapped reads") label(2 "Lowess smoother"))


lowess total_reads time_hr, by(pool_str, note("")legend(at(4)) title("") ) xtitle("Sequencing time")  ytitle("Total reads") legend(order(1 2) label(1 "Total reads") label(2 "Lowess smoother"))


sort sample_id time_hr

twoway line mapped_reads time_hr if pool_str=="1", by(sample_id)

twoway line mapped_reads time_hr if pool_str=="2", by(sample_id)

twoway line mapped_reads time_hr if pool_str=="3", by(sample_id)

twoway line mapped_reads time_hr if pool_str=="4", by(sample_id)

twoway line mapped_reads time_hr if pool_str=="5", by(sample_id)

twoway line mapped_reads time_hr if pool_str=="6", by(sample_id)

twoway line mapped_reads time_hr if pool_str=="7", by(sample_id)

twoway line mapped_reads time_hr if pool_str=="8", by(sample_id)

twoway line mapped_reads time_hr if pool_str=="9", by(sample_id)


xtset sample_id time_hr, delta(6)

encode(conc), generate(conc_cat)



* xtgee mapped i.time_hr i.conc_cat, family(binomial total_reads) corr(unstructured) vce(robust) /* doesn't converge*/

* xtgee mapped i.time_hr i.conc_cat, family(binomial total_reads) corr(ar 1) vce(robust) /* doesn't converge*/


* xtgee mapped i.time_hr i.conc_cat, family(binomial total_reads) corr(independent) vce(robust) /*main effects model*/


xtgee mapped i.time_hr##i.conc_cat, family(binomial total_reads) corr(independent) vce(robust) 

* xtgee mapped i.time_hr##i.conc_cat if conc!="Neg", family(binomial total_reads) corr(independent) vce(robust) /* does not include the Neg group*/


testparm i.time_hr

testparm i.conc_cat

testparm i.time_hr#i.conc_cat /* we cant use with main effects model*/

testparm  i.conc_cat i.time_hr#i.conc_cat /*cant use with main effects model*/

margins time_hr#conc_cat, expression(exp(predict(xb))/(1+exp(predict(xb))))
marginsplot, yscale(range(0 1)) ytitle("Proportion of 16S mapped reads") ytick(0 0.2 0.4 0.6 0.8 1) ylabel(0 "0" 0.2 "0.2" 0.4 "0.4" 0.6 "0.6" 0.8 "0.8" 1 "1")

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
 ytitle("Reads mapping to 16S Region") ///
 xtitle("Sequencing time") ///
 plotopts(mcolor(#4B0082) lcolor(#4B0082)) ///
 recastci(rarea) ciopts(color(green%10))

contrast ar.time_hr@i.conc_cat, mcompare(sidak)






/* Creation of incident reads (mapped and total): new reads from a specific time span*/

sort sample_id time_hr

by sample_id: generate new_mapped=mapped[_n] - mapped[_n-1]


by sample_id: replace new_mapped=mapped if time_hr==6



by sample_id: generate new_total=total_reads[_n] - total_reads[_n-1]


by sample_id: replace new_total=total_reads if time_hr==6


*/Generalized Estimating Equations: logistic regression on proportion of new reads (incident reads) */


* xtgee new_mapped i.time_hr i.conc_c, family(binomial new_total) corr(independent) vce(robust) /* zeros in the new_total variable*/




*/Generalized Estimating Equations: linear regression on proportion of new reads (incident mapped and total reads) */

xtgee new_mapped i.time_hr##i.conc_cat, corr(independent) vce(robust)

testparm i.time_hr

testparm i.conc_cat

testparm i.conc_cat#i.time_hr

margins i.time_hr#i.conc_cat

marginsplot, by(conc_cat) title("") ytitle("Total new 16S mapped reads") plotopts(mcolor(black) lcolor(black))recastci(rarea) ciopts(color(green%10))

contrast ar.time_hr@conc_cat, mcompare(sidak)
contrast r.conc_cat@time_hr, mcompare(sidak)


























