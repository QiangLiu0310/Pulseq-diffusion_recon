EPI GE BOLD data
3T Siemens Trio
32 Channel Head coil
TR 3sec
TE 30msec
Echo Spacing 770usec
2.0mm isotropic resolution

96x60   matrix size
192x120 (mm) FOV
32 slices, 2mm thick

2015-05-18  MID902

Filelist:
dpg_acs1.nd     ACS data, with the first readout acquired using RO+
dpg_acs2.nd     ACS data, with the first readout acquired using RO-
dpg_C.nd        The PLACE ACS image, formed from a combination of dpg_acs1 and dpg_acs2

dpg_k.nd        one frame of R=2 EPI data
dpg_k_pcref.nd  navigator data for dpg_k, to perform traditional Nyquist ghost correction
y.nd            NGC coefficients, derived from dpg_k_pcref and used with phzapply.m
y_acs2.nd       NGC coefficients, derived from A and B and used with phzapply.m

dpg_cal.m       the DPG calibration code
dpg_recon.m     the DPG reconstruction code
phzapply.m      a script to apply standard Nyquist ghost correction coefficients
fif.m           a 2D FFT macro script
readnd.m        script to read the .nd data files

runthis.m       the main demo script
runthis2.m      an alternate demo script (w/o GESTE/PLACE acs data)

