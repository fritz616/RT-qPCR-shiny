# RT-qPCR-shiny
Simple script for ddCt analysis of qPCR data for pharmacologists. It calculates fold change based on raw CT values for reference gene and gene of interest for specific treatment groups and returns XLSX file containing summary statistics (including means, standard deviations, N, SEM, min, max and normality distributtion) and descriptive statistics from ANOVA and Tukey post-hoc analysis of between group differences.

To use the app:

1. format your data accordingly to example below
2. **do not change the column name in first column - treatment**, do not use whitespace or special characters in headers or treatment groups
3. save you file as CSV and upload it (use semicolon as column separator and comma as decimal separator)
4. Type the name of the reference group, eg. VEH
5. Calculate duplicate means or remove outliers if required
6. You have all the results below and graph in the tab in upper-right corner
7. To save the graph just right-click it and choose save as and add image type extensions, for example *.jpg or *.png
