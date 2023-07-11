# BSISO Index based on Extended-EOF Analysis

By Sandro W. Lubis and Ziming Chen (PNNL)

This script is used to construct BSISO Index based on EEOF (lags -10, -5, and 0) following Kikuchi (2021) https://doi.org/10.2151/jmsj.2021-045

Data: Daily NOAA interpolated OLR dataset (Liebmann & Smith, 1996)
Steps (Kikuchi, 2021):
1. cal_bf_filter.ncl:  Apply 25~90-day Lanczos bandpass filter (Duchon, 1979) to OLR 
2. Perform EEOF analysis with three-time lags (–10, –, and 0 days) on the intraseasonal OLR data 

<p align="center">
  <img src="https://github.com/sandrolubis/BSISO_Index/blob/main/bsiso_plots.png" width="1000">
</p>
