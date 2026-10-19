# EpiModelHIV

EpiModelHIV is an R package for simulating HIV and other sexually transmitted infections over dynamic sexual networks. It extends [EpiModel](https://www.epimodel.org/), our general platform for network-based epidemic modeling, with modules for the populations, partnerships, behaviors, and prevention and care services in our HIV and STI research.

This repository holds the version of EpiModelHIV used in each of the published papers listed below, so that their results can be reproduced. Development happens in a private repository, and we add a snapshot here when a paper is published.

## How this repository is organized

- Each paper has a release with a tag, such as `EDP-JAIDS`. The tag points to the exact engine code behind the paper's results.
- Most releases also have a branch with the same code. Apart from notes added to their READMEs, snapshot branches are frozen: they are not updated or merged.
- Each paper's analysis code (network estimation, calibration, scenarios, and the tables and figures) is in its own project repository, linked in the table below.
- The default branch, `CombPrev`, is the snapshot for Jenness et al., AIDS 2020. Installing without a tag gets that version, so always give the tag for the paper you need.

## Reproducing a paper

1. Find the paper in the table below and clone its analysis repository.
2. If the analysis repository has an `renv.lock` file, restore the package versions used for the paper by running this in R from the project folder:

   ```r
   install.packages("renv")
   renv::restore(exclude = "ARTnetData")
   ```

   This installs EpiModelHIV from the paper's release, along with the versions of EpiModel, the statnet packages, and the other R packages recorded in the lockfile. Use the R version listed in the table if you can.

3. If there is no lockfile, install the engine by its tag:

   ```r
   remotes::install_github("EpiModel/EpiModelHIV@<tag>")
   ```

   These snapshots were built for EpiModel 1.x and depend on EpiModelHPC, which is on GitHub only ([EpiModel/EpiModelHPC](https://github.com/EpiModel/EpiModelHPC)). Current releases of EpiModel, EpiModelHPC, and tergmLite are unlikely to work with them, so install versions from around the paper's date.

ARTnetData, the data package for the ARTnet survey of men who have sex with men, is private. It is needed only to re-estimate the network models; see the [ARTnet repository](https://github.com/EpiModel/ARTnet) for how to request access.

## Papers

Newest first. The last column gives the EpiModel and R versions recorded in the analysis repository's `renv.lock`, or, where there is no lockfile, the oldest EpiModel version the snapshot accepts.

| Paper | Topic | Release | Analysis code | renv.lock |
|---|---|---|---|---|
| Chandra et al., [J Acquir Immune Defic Syndr 2026](https://doi.org/10.1097/QAI.0000000000003934) | Event-driven PrEP | [EDP-JAIDS](https://github.com/EpiModel/EpiModelHIV/releases/tag/EDP-JAIDS) | [EventDrivenPrEP](https://github.com/EpiModel/EventDrivenPrEP) | [EpiModel 2.3.2, R 4.3.2](https://github.com/EpiModel/EventDrivenPrEP/blob/main/renv.lock) |
| Le Guillou et al., [J Infect Dis 2026](https://doi.org/10.1093/infdis/jiag369) | Low-barrier PrEP access and monitoring | [LowBarrier-JID](https://github.com/EpiModel/EpiModelHIV/releases/tag/LowBarrier-JID) | [offscript](https://github.com/EpiModel/offscript) | [EpiModel 2.5.0, R 4.5.2](https://github.com/EpiModel/offscript/blob/main/renv.lock) |
| Chandra et al., [Lancet Reg Health Am 2026](https://doi.org/10.1016/j.lana.2026.101620) | PrEP reinitiation, uptake, and persistence in Boston | [PrEPReinit-LRHA](https://github.com/EpiModel/EpiModelHIV/releases/tag/PrEPReinit-LRHA) | [PrEP-Reinitiation](https://github.com/EpiModel/PrEP-Reinitiation) | [EpiModel 2.5.0, R 4.4.2](https://github.com/EpiModel/PrEP-Reinitiation/blob/main/renv.lock) |
| Onwubiko et al., [BMC Health Serv Res 2026](https://doi.org/10.1186/s12913-026-15074-0) | Mental health care in PrEP programs | [MentalHealth-BMCHSR](https://github.com/EpiModel/EpiModelHIV/releases/tag/MentalHealth-BMCHSR) | [StigmaDepressionHIV](https://github.com/EpiModel/StigmaDepressionHIV) | [EpiModel 2.4.0, R 4.3.2](https://github.com/EpiModel/StigmaDepressionHIV/blob/main/renv.lock) |
| Onwubiko et al., [Sex Transm Dis 2025](https://doi.org/10.1097/OLQ.0000000000002167) | Partner services for people previously diagnosed with HIV | [PartServDx-STD](https://github.com/EpiModel/EpiModelHIV/releases/tag/PartServDx-STD) | [PartnerServiceYr3](https://github.com/EpiModel/PartnerServiceYr3) | [EpiModel 2.3.2, R 4.2.3](https://github.com/EpiModel/PartnerServiceYr3/blob/main/renv.lock) |
| Jones et al., [Sex Transm Dis 2023](https://doi.org/10.1097/OLQ.0000000000001724) | HIV infections per gonorrhea or chlamydia infection | [GCCT-STD](https://github.com/EpiModel/EpiModelHIV/releases/tag/GCCT-STD) | [PAFGuidelines](https://github.com/EpiModel/PAFGuidelines) | [EpiModel 2.0.5, R 4.2.0](https://github.com/EpiModel/PAFGuidelines/blob/main/renv.lock) |
| Jenness et al., [Sex Transm Dis 2022](https://doi.org/10.1097/OLQ.0000000000001711) | HIV partner services | [PartServ-STD](https://github.com/EpiModel/EpiModelHIV/releases/tag/PartServ-STD) | [CombPrevNet](https://github.com/EpiModel/CombPrevNet) | [EpiModel 2.0.4, R 4.1.2](https://github.com/EpiModel/CombPrevNet/blob/master/renv.lock) |
| Jones et al., [Sex Transm Dis 2022](https://doi.org/10.1097/OLQ.0000000000001685) | Gonorrhea and chlamydia screening and HIV incidence | [GCCT-STD](https://github.com/EpiModel/EpiModelHIV/releases/tag/GCCT-STD) | [PAFGuidelines](https://github.com/EpiModel/PAFGuidelines) | [EpiModel 2.0.5, R 4.2.0](https://github.com/EpiModel/PAFGuidelines/blob/main/renv.lock) |
| Wheatley et al., [J Acquir Immune Defic Syndr 2022](https://doi.org/10.1097/QAI.0000000000002921) | Cost-effectiveness of PrEP interventions | [PrEPCEA-JAIDS](https://github.com/EpiModel/EpiModelHIV/releases/tag/PrEPCEA-JAIDS) | [PrEP-CEA](https://github.com/EpiModel/PrEP-CEA) | none (EpiModel 1.7.0 or later) |
| Le Guillou et al., [J Acquir Immune Defic Syndr 2021](https://doi.org/10.1097/QAI.0000000000002781) | PrEP and ART in San Francisco | [SFPrEP-JAIDS](https://github.com/EpiModel/EpiModelHIV/releases/tag/SFPrEP-JAIDS) | [SFO_PrEP](https://github.com/EpiModel/SFO_PrEP) | [EpiModel 1.8.0, R 4.0.2](https://github.com/EpiModel/SFO_PrEP/blob/main/renv.lock) |
| Jenness et al., [AIDS 2021](https://doi.org/10.1097/QAD.0000000000002909) | Investment in the PrEP care cascade | [Optim-AIDS](https://github.com/EpiModel/EpiModelHIV/releases/tag/Optim-AIDS) | [PrEP-Optimize](https://github.com/EpiModel/PrEP-Optimize) | [EpiModel 1.7.5, R 4.0.3](https://github.com/EpiModel/PrEP-Optimize/blob/master/renv.lock) |
| Jenness et al., [J Infect Dis 2021](https://doi.org/10.1093/infdis/jiab051) | COVID-19 sexual distancing and service interruptions | [SexDist-JID](https://github.com/EpiModel/EpiModelHIV/releases/tag/SexDist-JID) | [SexualDistancing](https://github.com/EpiModel/SexualDistancing) | [EpiModel 1.8.0, R 4.0.2](https://github.com/EpiModel/SexualDistancing/blob/master/renv.lock) |
| Maloney et al., [J Infect Dis 2021](https://doi.org/10.1093/infdis/jiaa552) | Long-acting injectable and daily oral PrEP | [LAPREP-JID](https://github.com/EpiModel/EpiModelHIV/releases/tag/LAPREP-JID) | [injectable-prep](https://github.com/EpiModel/injectable-prep) | [EpiModel 1.7.5, R 4.0.2](https://github.com/EpiModel/injectable-prep/blob/master/renv.lock) |
| Jenness et al., [AIDS 2020](https://doi.org/10.1097/QAD.0000000000002681) | Integrated HIV prevention and care for the Ending the HIV Epidemic goals | [EHE-AIDS](https://github.com/EpiModel/EpiModelHIV/releases/tag/EHE-AIDS) | [CombPrev](https://github.com/EpiModel/CombPrev) | [EpiModel 1.7.3, R 4.0.2](https://github.com/EpiModel/CombPrev/blob/master/renv.lock) |
| Weiss et al., [Sex Transm Dis 2019](https://doi.org/10.1097/OLQ.0000000000001058) | Expedited partner therapy | [ept](https://github.com/EpiModel/EpiModelHIV/releases/tag/ept) | [EPT](https://github.com/EpiModel/EPT) | none (EpiModel 1.2.7 or later) |
| Weiss et al., [Open Forum Infect Dis 2019](https://doi.org/10.1093/ofid/ofz405) | STI screening coverage and frequency | [STITest-OFID](https://github.com/EpiModel/EpiModelHIV/releases/tag/STITest-OFID) (replaces `ofid1`) | [stitestguidelines](https://github.com/EpiModel/stitestguidelines) | none (EpiModel 1.2.7 or later) |
| Jones et al., [Sex Transm Dis 2019](https://doi.org/10.1097/OLQ.0000000000000980) | Share of HIV infections attributable to gonorrhea and chlamydia | [paf1](https://github.com/EpiModel/EpiModelHIV/releases/tag/paf1) | [sti_paf](https://github.com/EpiModel/sti_paf) | none (EpiModel 1.2.7 or later) |
| Jenness et al., [Am J Epidemiol 2019](https://doi.org/10.1093/aje/kwy230) | PrEP care gaps and racial disparities | [p4](https://github.com/EpiModel/EpiModelHIV/releases/tag/p4) | [PrEPdisparities](https://github.com/EpiModel/PrEPdisparities) | none (EpiModel 1.2.7 or later) |
| Jenness et al., [Clin Infect Dis 2017](https://doi.org/10.1093/cid/cix439) | Gonorrhea and chlamydia after PrEP | [p3](https://github.com/EpiModel/EpiModelHIV/releases/tag/p3) | [stiPrEP](https://github.com/EpiModel/stiPrEP) | none (EpiModel 1.2.7 or later) |
| Jenness et al., [PLoS One 2017](https://doi.org/10.1371/journal.pone.0169484) | PrEP risk compensation | [p2](https://github.com/EpiModel/EpiModelHIV/releases/tag/p2) | [HIV-Risk-Comp](https://github.com/EpiModel/HIV-Risk-Comp) | none (EpiModel 1.2.7 or later) |
| Jenness et al., [J Infect Dis 2016](https://doi.org/10.1093/infdis/jiw223) | CDC PrEP guidelines | [p1](https://github.com/EpiModel/EpiModelHIV/releases/tag/p1) | [PrEPGuidelines](https://github.com/EpiModel/PrEPGuidelines) | none (EpiModel 1.2.0 or later) |

## Other public branches

These branches from the University of Washington group support published papers but have no release, and they keep their full development history.

| Branch | Papers |
|---|---|
| [Kenya-TM](https://github.com/EpiModel/EpiModelHIV/tree/Kenya-TM) | Hamilton et al., [J Acquir Immune Defic Syndr 2022](https://doi.org/10.1097/QAI.0000000000003013); Hamilton et al., [Epidemics 2023](https://doi.org/10.1016/j.epidem.2023.100696) |
| [Hometest](https://github.com/EpiModel/EpiModelHIV/tree/Hometest) | Hamilton et al., [Epidemics 2021](https://doi.org/10.1016/j.epidem.2021.100518) |
| [campcl](https://github.com/EpiModel/EpiModelHIV/tree/campcl) | Hamilton et al., [Am J Public Health 2018](https://doi.org/10.2105/AJPH.2018.304471); Hamilton et al., [PLoS One 2019](https://doi.org/10.1371/journal.pone.0217315); Hamilton et al., [J Adolesc Health 2021](https://doi.org/10.1016/j.jadohealth.2020.06.041); Katz et al., [Sex Transm Dis 2021](https://doi.org/10.1097/OLQ.0000000000001485) |
| [MSM-race-disparities](https://github.com/EpiModel/EpiModelHIV/tree/MSM-race-disparities) | Goodreau et al., [Lancet HIV 2017](https://doi.org/10.1016/S2352-3018(17)30067-X) |

The remaining branches (`master`, `main`, `shamp`, `shamp_params`, `shamp_summary_stats`, and `stiratetest`) and the tags `croi2017` and `p3LID` are from earlier development and are not tied to a published paper.

## About the model

EpiModelHIV simulates a population of individuals and their sexual partnerships over time. In most of the papers above, partnerships form and dissolve on separate main, casual, and one-time networks, using temporal exponential-family random graph models from [statnet](https://statnet.org/) ([Krivitsky and Handcock 2014](https://doi.org/10.1111/rssb.12014)). For the papers since 2020, these network models were estimated from the ARTnet survey. At each time step, modules for demography, sexual behavior, HIV and STI transmission, and the HIV prevention and care continuum (testing, antiretroviral therapy, and PrEP) act on the population and its networks. The modules differ across snapshots, so each paper's methods and appendix describe the version it used.

The [EpiModel website](https://www.epimodel.org/) has documentation, tutorials, and course materials for network-based epidemic modeling, including our [SISMID course](https://epimodel.github.io/sismid/).

## Citation

Please cite the paper whose results you are using. For the modeling framework, also cite:

> Jenness SM, Goodreau SM, Morris M. EpiModel: An R Package for Mathematical Modeling of Infectious Disease over Networks. _J Stat Softw._ 2018; 84(8): 1-47. https://doi.org/10.18637/jss.v084.i08

## Questions

Questions about a specific analysis are best sent to the paper's corresponding author. If a snapshot does not install as described here, please open an issue in this repository. For EpiModel itself, see the [EpiModel repository](https://github.com/EpiModel/EpiModel).

## License

GPL-3. See [LICENSE.md](LICENSE.md).
