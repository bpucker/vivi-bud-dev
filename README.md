# vivi-bud-dev
Collection of scripts related to the bud development study in Vitis vinifera.

## Gene expression plot scripts

ap1_ail2_soc1a_plot.py

Generates a  gene expression plot of VviAP1 (VIT_201s0011g00100), VviAIL2 (VIT_209s0002g01370), and VviSOC1a (VIT_215s0048g01250) across the investigated time. CPMs calculated based on RNA-Seq read mappings with STAR are displayed.


dyl1_plot.py

Generates a gene expression plot of VviDRM1/VviDYL1 (VIT_210s0003g00090) across the investigated time. CPMs calculated based on RNA-Seq read mappings with STAR are displayed.


## Get heatshock genes

This script extracts potential heat shock genes based on the functional annotation.


## Correlation of heatshock gene expression with temperature

weather_plot.py

Generates a gene expression plot for potential heatshock genes to show the correlation with the temperature profile.

## Heatmap generation

mads.py

myb.py

wrky.py

These scripts generate gene expression heatmaps of the respective gene family based on RPKM values of the RNA-Seq experiment.


## References

Pucker B., Schwandner A., Becker S., Hausmann L., Viehöver P., Töpfer R., Weisshaar B., Holtgräwe D. (2020). RNA-Seq Time Series of Vitis vinifera Bud Development Reveals Correlation of Expression Patterns with the Local Temperature Profile. bioRxiv 2020.10.18.344176; doi: [10.1101/2020.10.18.344176](https://doi.org/10.1101/2020.10.18.344176).

