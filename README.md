# tht_analysis
Levine Lab
ThT Analysis Script
By Bryan Bogin
Analyze, fit, and plot ThT data in digestible manner. Designed for visualization, normalization, and customization of datasets.

Prepare ThT data in .csv or excel format. If headers are different, you must select the appropriate strings to aquire the data of interest.

Can average replicates of ThT

Do not recommend averaging curves if noticible lag time differences are observed.
In that case, turn on the lmfit function to fit the curve with a sigmoid.
Then, normalized based on lag time and t50s.

