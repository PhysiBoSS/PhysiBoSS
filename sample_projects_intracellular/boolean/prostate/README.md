# About

This repository is an update of the prostate model to version 2.2.3 of PhysiBoSS. The project does not rely on the cell rules to perform the inhibition, but on an ad-hoc cell sensitivity that depends on the cell line.

A prior version of this model was used to simulate the drug response in prostate cell lines in Figures 6, 7 and 8 of [the PhysiBoSS 2.0 paper](https://www.nature.com/articles/s41540-023-00314-4).

# Getting started

The XML in the `config` folder are: 
- PhysiCell_settings_LNCaP.xml: LNCaP cell line simulation
- PhysiCell_settings_LNCaP_Luminespib_2.xml: LNCaP cell line simulation adding Luminespib drug (anti_HSPs).
- PhysiCell_settings_LNCaP_Pictilisib_2.xml: LNCaP cell line simulation adding Pictilisib drug (anti_PI3K).
- PhysiCell_settings_LNCaP_Luminespib_Pictilisib_2.xml: LNCaP cell line simulation adding Luminespib and Pictilisib drugs.

Files in the `scripts`  folder have post processing scripts.

# Authors:

Arnau Montagud, Annika Meert, Gerard Pradas and Miguel Ponce de Leon, BSC-CNS

Contact: arnau.montagud at csic.es
