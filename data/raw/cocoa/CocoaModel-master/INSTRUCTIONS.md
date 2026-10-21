# Cocoa bean fermentation model
Datasets and codes supporting the article: "*A mathematical model of cocoa bean fermentation*" by Mauricio Moreno-Zambrano, Sergio Grimbs, Matthias S. Ullrich and Marc-Thorsten Hütt

## Contents
Within the folder "simulation":
* Three R scripts, one for each collection of data, containing the reported data as well as the parameters under which the Bayesian framework is ran.
  * **camu.R**: data reported by Camu *et al.* (2007).
  * **p_box1.R**: data reported by Papalexandratou *et al.* (2011) for trial 'box 1'.
  * **p_box2.R**: data reported by Papalexandratou *et al.* (2011) for trial 'box 2'.

* One 'Stan' script (**model.stan**) file containing the model in 'Stan' language.

Within the folder "results":
* One R script named **table.R** that generates LaTeX code for a table containing the re-scaled parameters, similar to Table 3 of our article.

## R Script
Independently of the collection of data, the R scripts follows a structure of three blocs: (1) Header, (2) Data and (3) Execution.

* **Header block**: Common for all three scripts, it contains the basic settings of the script, i.e., removal of objects from R environment, setting up of a working directory and call to RStan library.

* **Data block**: Here, time-series for each dataset are stated. Values shown for microbial groups' time-series are the result of a previous step of conversion of *log*(CFU) to mass units (see section 2.2 of the article for further explanation).

* **Execution block**: Here, the parameters for the Bayesian parameter estimation framework are specified.

## General instructions
1. Download the R script(s) in which you are interested in, as well as the 'Stan' script (**model.stan**) within a same folder.
2. In the R script, change the path according to which folder you saved the Scripts in.
3. In case that you have a multicore processing machine with at **least 4 cores and plenty of RAM**, go to **step 6**. Otherwise, go to step 4.
4. Delete the line of code in the header block that says: `rstan_options(auto_write = TRUE)`.
5. Delete the line of code in the execution block that says: `cores = min(4,  parallel::detectCores())`.
6. Run the R script

## Final remarks
* Once simulations are finished, these will be given in non-scaled quantities (see section 2.4.2 of the article for further explanation). Run script **table.R** for generating re-scaled estimates to their original units.
* Note that, due to the stochastic nature of the MCMC-No-U-Turn sampler, the parameter estimates might be slightly different to those reported in our article.
* For further questions regarding these codes, you can contact Mauricio Moreno-Zambrano at m.morenozambrano@jacobs-university.de.

## References
Camu, N., De Winter, T., Verbrugghe, K., Cleenwerck, I., Vandamme, P., Takrama, J. S., Vancanneyt, M. & De Vuyst, L. 2007 Dynamics and biodiversity of populations of lactic acid bacteria and acetic acid bacteria involved in spontaneous heap fermentation of cocoa beans in Ghana. Applied and Environmental Microbiology, 73(6), 1809--1824. [(doi:10.1128/AEM.02189-06)](http://aem.asm.org/content/73/6/1809.full).

Papalexandratou, Z., Vrancken, G., Bruyne, K. D., Vandamme, P. & Vuyst, L. D. 2011 Spontaneous organic cocoa bean box fermentations in Brazil are characterized by a restricted species diversity of lactic acid bacteria and acetic acid bacteria. Food Microbiology, 28(7), 1326--1338. [(doi:10.1016/j.fm.2011.06.003)](http://www.sciencedirect.com/science/article/pii/S0740002011001407?via%3Dihub).
