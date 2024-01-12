dpgfullpath = mfilename('fullpath');
dpgname = mfilename;
dpgpath = regexprep(dpgfullpath,dpgname,'');
addpath([ dpgpath '/mbin/' ]);
addpath([ dpgpath '/mathlib/conjgrad/' ]);
addpath([ dpgpath '/mathlib/tensor/' ]);
addpath([ dpgpath '/mri/mtlib/' ]);
addpath([ dpgpath '/mri/pmri/' ]);
addpath([ dpgpath '/mri/sie_pcm/' ]);
addpath([ dpgpath '/mri/siemens/' ]);
addpath([ dpgpath '/mri/siemens/mapVBVD/' ]);
addpath([ dpgpath '/nifti/' ]);
addpath([ dpgpath '/yj_function/' ]);
clear dpgpath dpgname dpgfullpath