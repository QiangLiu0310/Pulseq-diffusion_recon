# Pulseq-diffusion_recon
Pulseq-Siemens data reconstruction
using mapVBVD https://github.com/CIC-methods/FID-A/tree/master/inputOutput/mapVBVD to read the data
Pulseq toolbox is required to calculate the trajectory for EPI
# EPI recon includes: 
1. ramp sampling intepolaration (Pulseq code) 
2. Ghost correction (LPC based on Dr.William Scott Hoge's code) 
3. GRAPPA reconstruction
4. Homodyne for Partial Fourier recon
5. Apodization (based on Dr.Bilgic Berkin's toolbox @Pulseq_Mprage_Recon_Toolbox)
# by Qiang Liu et al., Brigham and Women's Hospital, Harvard Medical School
qliu30@mgh.harvard.edu
qiangliu2019@163.com
