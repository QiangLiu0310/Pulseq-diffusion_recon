function vrgf = form_vrgf( hdr )

inres  = hdr.Config.DestSamples;
outres = hdr.Config.BaseResolution;
t_ru   = hdr.Config.RampupTime;
t_ds   = hdr.Config.DelaySamplesTime;
t_ft   = hdr.Config.FlattopTime;
t_sw   = hdr.Config.ADCDuration;

vrgf = mk_vrgf_dat( inres, outres, t_ru, t_ds, t_ft, t_sw );

return;