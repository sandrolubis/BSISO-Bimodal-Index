# BSISO Index based on EEOF Analysis

This script is used to construct BSISO Index based on EEOF (lags -10, -5, and 0) following Kikuchi (2021) https://doi.org/10.2151/jmsj.2021-045

Data: Daily NOAA interpolated OLR dataset (Liebmann & Smith, 1996)
Steps (Kikuchi, 2020, 2021):
- Apply 25~90-day Lanczos bandpass filter (Duchon, 1979) to OLR 
- Perform EEOF analysis with three-time lags (–10, –, and 0 days) on the intraseasonal OLR data


![image](https://github.com/sandrolubis/BSISO_Index/assets/26489636/f7f9a0e6-eff5-4e5d-9692-bca65dad5a6c)

<p align="center">
  <img src="https://github.com/sandrolubis/BSISO_Index/blob/main/bsiso_plots.png" width="1000">
</p>
