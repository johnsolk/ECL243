Reference:
    https://github.com/mojaveazure/angsd-wrapper/wiki/Tutorial

* Estimate the site frequency spectrum for maize and teosinte
(Bonus: plot the SFS -- why does it look weird?)

The site frequency spectrum (num of sites vs. num of variants at that site) of a neutral population will have a large num of sites with a small number of variants at that site. As the number of sites increases, there will be a small number of variants. In the Maize population, there is a slight U shape with slightly higher number of sites with number of variants. This could be due to drift, introducing diversity into the population. 

* Get summary statistics diversity for maize and teosinte
Use angsd-wrapper Thetas Thetas_Config. Report the overall diversity values for
each taxa.

mean theta pi

https://github.com/ljcohen/ECL243/blob/master/Assignment3/Teosinte_Thetas/Teosinte.thetas.graph.me
https://github.com/ljcohen/ECL243/blob/master/Assignment3/Maize_Thetas/Maize.thetas.graph.me

* Explain the differences you see in nucleotide diversity and Tajima's D between
maize and teosinte

plots: http://rpubs.com/rossibarra/ecl243_diversity

Nucleotide diversity across the 10kb window of the chromosome 10 decreases around 17.5 and then plateaus. Overall, diversity is positive, with theta pi (common variants) greater than theta W (rare variants). This shoudl be zero under balancing selection. When Tajima's D is positive, this could mean that there is a deficit of rare variants and the population could be bottlenecking.
