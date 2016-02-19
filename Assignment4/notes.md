#Questions:

* How and why do your results differ from Hufford et al. 2013?

In [Hufford et al. 2013](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003477#pgen-1003477-g002), Figure 2D shows a STRUCTURE analysis with K=2 between mexicana and maize individuals. Mexicana and Maize were both used as a reference, with a chosen K=2 value, meaning two ancestral populations. Slightly asymmetric introgression was seen across individuals, with more highland maize derived from mexicana (19% versus 12% of mexicana from maize). The authors indicated higher K values showed admixture in mexicana populations but not in maize, suggesting that gene flow from mexicana into maize may have been more ancient.

According to this [post](https://www.biostars.org/p/138371/), [STRUCTURE](http://pritchardlab.stanford.edu/structure.html) ([paper](http://pritchardlab.stanford.edu/publications/pdfs/PritchardEtAl00.pdf)) and [ADMIXTURE](https://www.genetics.ucla.edu/software/admixture/) ([paper](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-246)) are both used for population structure, where ADMIXTURE uses a maximum likelihood framework and STRUCTURE uses bayesian approach. Admixture estimates individual ancestries by computing maximum likelihood estimates in a parametric model. Given K ancestral populations, the allele's "success" probability depends on the fraction inherited from those population.

For this analysis: 

    source ~/.bash_profile
    angsd-wrapper SFS ./Site_Frequency_Spectrum_Config
    angsd-wrapper Admixture ./Admixture_Config
    angsd-wrapper Abbababa ./Abbababa_Config

an admixture framework was used with several K values below: 

![image](https://raw.githubusercontent.com/ljcohen/ECL243/master/Assignment4/all_Admixture.png)

These results are slightly different from Hufford et al. 2013. In addition to Maize and Mexicana, Teosinte samples are included. Tripsacum was used as a reference. Introgression appears more variable across samples. Mexicana and Maize appear diferent, with more introgression between Teosinte and Maize, with some overlaps between Teosinte (parviglumis) and Mexicana.

(((Mays,Parv),Mex),Trip)

* Do the D statistics jive with the admixture plot? Why might these two methods give
different answers?



