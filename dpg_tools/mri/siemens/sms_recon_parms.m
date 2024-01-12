function [SMS,FOVshift,NSlc,Ngroup,PhaseShift] = sms_recon_parms(hdr)

SMS      = hdr.SMS;
FOVshift = hdr.FOVshift;
NSlc     = hdr.NSlc;

Ngroup   = NSlc/SMS;   % number of measured groups per volume

if ( FOVshift == 1 ),
  PhaseShift = 0;
else,
  PhaseShift = 2*pi/FOVshift;
end;
