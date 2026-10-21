# Sources

Trust levels:
- **A**: peer-reviewed, with data available
- **B**: peer-reviewed, data only in figures (digitize)
- **C**: practitioner or community knowledge, used only as a loose check

## Shared engine equations

| Source | Implementation reference |
|---|---|
| Rosso et al. 1995, doi:10.1128/aem.61.2.610-616.1995 ([paper](https://pbil.univ-lyon1.fr/members/lobry/articles/AEM_1995_61_610.pdf)) | Equation 6: cardinal temperature; equation 7: cardinal pH. P1.1 transcribes equations, not library code. |
| Monod 1949, doi:10.1146/annurev.mi.03.100149.002103 ([paper](https://www.annualreviews.org/content/journals/10.1146/annurev.mi.03.100149.002103)) | Substrate-limited growth in P1.3. |
| Contois 1959, doi:10.1099/00221287-21-1-40 ([paper](https://www.microbiologyresearch.org/content/journal/micro/10.1099/00221287-21-1-40)) | Biomass-dependent substrate limitation in P1.3. |

The Gänzle and Breidt sources below also supply generic equations. Numerical
constants in engine tests and examples are explicitly synthetic and do not
represent extracted organism or food parameters.

## Estimation methods

| Source | Use |
|---|---|
| [OptimizationOptimJL documentation](https://docs.sciml.ai/Optimization/stable/optimization_packages/optim/) | Bounded Fminbox/LBFGS fitting with automatic derivatives (P1.5). |
| [Turing Bayesian differential equations](https://turinglang.org/docs/tutorials/bayesian-differential-equations/) and [sampling API](https://turinglang.org/docs/usage/sampling-options/) | NUTS through hand-written ODEs, explicit warmup and retained-chain handling (P1.6). |
| Vehtari et al., rank-normalization and improved R̂ ([paper](https://arxiv.org/abs/1903.08008)) and [MCMCDiagnosticTools](https://turinglang.org/MCMCDiagnosticTools.jl/) | Rank R̂, bulk/tail ESS, and diagnostic thresholds. |

These method references do not supply food-model parameters or prior widths.

P1.7 uses the [ForwardDiff Jacobian API](https://juliadiff.org/ForwardDiff.jl/stable/user/api/)
and profile-likelihood diagnostics following [Raue et al. (2009),
doi:10.1093/bioinformatics/btp358](https://pubmed.ncbi.nlm.nih.gov/19505944/).

## Cocoa

| Source | What it gives | License / access | Trust |
|---|---|---|---|
| Moreno-Zambrano, Grimbs, Ullrich & Hütt 2018, *R. Soc. Open Sci.*, doi:10.1098/rsos.180964 ([PMC6227950](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6227950/)) | Baseline ODE model (sugars, ethanol, lactic, acetic; yeast, LAB, AAB), fitted with Bayesian estimation | Open access | A |
| Data and code: doi:10.5061/dryad.321d33v ([Zenodo 4935403](https://zenodo.org/records/4935403)) | `CocoaModel.zip`: seven files, three embedded datasets, one Stan model; [P2.1 inventory](cocoa_original.md) | CC0 | A |
| Supplement: doi:10.6084/m9.figshare.7152614.v2 ([Figshare](https://rs.figshare.com/articles/journal_contribution/A_mathematical_model_of_cocoa_bean_fermentation_Supplementary_Material/7152614)) | S3 scaling derivation (pp. 7–10); S4–S6 normalized posterior tables (pp. 11–13). Download pinned with provenance; S3 confirms the fructose-yield conversion error in the archived table script. | CC BY 4.0 | A |
| Moreno-Zambrano, Ullrich & Hütt 2022 ([PMC8847890](https://pmc.ncbi.nlm.nih.gov/articles/PMC8847890)) | 32 model variants screened on 32 trials; 23 retained. [P2.5 audit](cocoa_benchmark.md) | CC BY 4.0 | A |
| 2022 source archive: doi:10.5281/zenodo.5510370 ([Zenodo](https://zenodo.org/records/5510370)) | 78 files; embedded observations for all 32 trials; 64 main Stan files in two observation families; posterior outputs absent | Zenodo `other-open`; no explicit license file; do not assume 2018 CC0 applies | A |
| 2022 supplement: doi:10.6084/m9.figshare.19149492.v1 ([Figshare](https://doi.org/10.6084/m9.figshare.19149492.v1)) | Mechanism equations, dataset catalog, convergence summaries, BMA weights and PSIS-LOO plots; archived with checksums | CC BY 4.0 | A |
| Moreno-Zambrano PhD thesis ([Constructor University](https://opus.constructor.university/frontdoor/index/index/docId/1001), [DNB PDF](https://d-nb.info/1240902832/34)) | Temperature via Arrhenius terms; chapter 2 reproduces the 2018 article. Figures 2.3–2.5 supply the archived figure crops and 168 digitized reference points for P2.3. | Public PDF; only CC BY article figure crops redistributed, not the full thesis | A/B |

## Sourdough


| Source | What it gives | License / access | Trust |
|---|---|---|---|
| Gänzle, Ehmann & Hammes 1998, AEM 64:2616 ([PMC106434](https://pmc.ncbi.nlm.nih.gov/articles/PMC106434)) | Growth of *L. sanfranciscensis* and *C. milleri* vs. T, pH, NaCl, lactate, acetate, ethanol | Open access | A |
| Gänzle temperature formula via rec.food.sourdough FAQ ([link](https://nyx.org/~dgreenw/whatistherelationshipbetwe.html)) | Secondary summary; implementation uses the original paper's Table 1, `a·x^b·e^(-c·x)` for positive c (see `secondary_models.md`) | Public web page | C; primary equation checked in P1.1 |
| Landis et al. 2021, *eLife* 10:e61644 ([article](https://elifesciences.org/articles/61644), [Dryad](https://datadryad.org/dataset/doi:10.5061/dryad.0p2ngf1z1)) | 500-starter community survey; rise rate and aroma for a subset | Open data | A |
| Sourdough gLV rock-paper-scissors study ([PMC10714952](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10714952/)) | Interaction model for a 5-species LAB–yeast community | Open access | A |
| The Fresh Loaf rise-time table ([comment 27403](https://www.thefreshloaf.com/comment/27403)) | Doubling, bulk and proof times vs. T and inoculation | Forum post | C |

## Vegetables

| Source | What it gives | License / access | Trust |
|---|---|---|---|
| Zabat et al. 2018, *Foods* 7:77 ([PMC5977097](https://pmc.ncbi.nlm.nih.gov/articles/PMC5977097)) | Sauerkraut community over time, about 21 °C | Open access | A/B |
| Jung et al. 2011 kimchi metagenome ([MGnify MGYS00000276](https://ebi.ac.uk/metagenomics/api/v1/studies/MGYS00000276)) | 29-day succession | Open data | A |
| Kim et al. 2022, 500-day kimchi at 4 °C ([Korea Univ.](https://pure.korea.ac.kr/en/publications/long-term-population-dynamics-of-viable-microbes-in-a-closed-ecos/)) | Long-term succession and metabolites | Check access | B |
| Breidt et al., buffer models for cucumber fermentation ([USDA-ARS PDF](https://ars.usda.gov/ARSUserFiles/60702500/Fermented%20Veg%20Manuscripts/p450.pdf)) | pH ↔ acid concentration | Public (US gov) | A/B |
| Breidt et al., *E. coli* O157:H7 survival in cucumber fermentation brines ([USDA-ARS PDF](https://ars.usda.gov/ARSUserFiles/60701000/Pickle%20Pubs/p373.pdf)) | 5-log reduction times vs. pH, acid, T | Public (US gov) | A/B |
| Sauerkraut production course notes ([Ankara Univ.](https://acikders.ankara.edu.tr/mod/resource/view.php?id=165281)) | Salt/temperature effects on hetero- vs. homofermenters; final acidity | Course notes | C (find primary source) |
| Pederson & Albury sauerkraut studies (NY Agric. Exp. Station) | Classic T and salt trials | To locate | — |
| makesauerkraut.com salt guidance ([link](https://www.makesauerkraut.com/how-much-salt-use-to-make-sauerkraut)) | Home practice: 1.5–2.5 % salt, adjust for temperature | Blog | C |

## Broader context (not in current plan)

- ComBase predictive microbiology database: https://www.combase.cc
- curatedFoodMetagenomicData (Carlino et al. 2024, *Cell*): https://pubmed.ncbi.nlm.nih.gov/39214080/
- Wolfe et al. 2014, cheese rind communities, *Cell* 158:422
