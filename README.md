# What is zeta_combiner?
Zeta_combiner is a Python script that combines replicate electrophoretic mobility / ELS / PALS / zeta potential distributions performed by a Brookhaven NanoBrook Omni. This was created because the Brookhaven NanoBrook Omni instrument/data processing software I had access to could not export all replicate scans of a sample into one file except as a difficult-to-digitize PDF report, but could export each replicate scan as a separate .csv or .xls file.
# Requirements

- Python 3.8+
- Matplotlib 3.1.2+
- NumPy 1.17.4+
- OpenPyxl 3.0.2+
- Pandas 0.25.3+
- XlsxWriter 1.2.6+
# How to use zeta_combiner
Zeta_combiner is run by copying zeta_combiner into a folder containing one or more sets of replicate data files that were manually exported from the Brookhaven software in .csv or .xls format with identical names except for their endings (Sample-1.csv, Sample-2.csv, etc. using dashes or underscores before the replicate number), then running the script from the terminal with the following syntax

`python zeta_combiner.py C:/example/path/to/data`

For each group of samples, zeta_combiner will export

- an .xlsx workbook summarizing the entire folder of data with a worksheet for each sample containing the raw and combined zeta potential distributions of that sample's replicate data files
- an .svg vectorized image of a plot of the combined zeta potential distribution.
# How zeta_combiner works
Zeta_combiner first identifies all of the unique zeta potential (x-axis) values that appear in at least one of the replicate data files. The Brookhaven NanoBrook Omni does not use a fixed bin size on the x-axis, so every replicate file has a different set of x-axis values. Zeta_combiner then calculates the mean, sample standard deviation, and number of power (y-axis) values that were collected at each unique x-axis value. Due to the variable bin size of the x-axes from each replicate data file and the fact that each replicate data file is normalized to a mode of 1.0 power, the mean zeta potential distribution can have large local discontinuities that can be reduced with a simple 1-D box blur / moving average. This smoothing can be adjusted in the source code through `nSmooths` (default 2) and `w` (default 1% of x-axis). The low and high ends of the x-axis range can be set by the user with `minZeta` and `maxZeta` (default -/+ 90 mV).
