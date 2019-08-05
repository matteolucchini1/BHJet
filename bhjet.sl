
% Create a S-lang library of the function using SLIRP:

% unix%> slirp -make mffjet.f
% unix%> make 

% Move mffjet-module.so to your S-lang modules directory (for me,
% /usr/local/isis/modules).  In your ISIS .isisrc, or other S-lang
% input file, run all the code below.

import("bhjet");

define bhjet_fit(lo,hi,par)
{
   variable ear = [_A(hi),_A(lo[0])];
   variable nein = length(lo), photar=Double_Type[nein], photer=@photar;
   
   variable pars = par[[1:28]];

   % !!!  It is important to match the type for arrays that are input, !!!
   % !!!  and recieve return values.  I.e., phot & energ go in and     !!!
   % !!!  out as real*8 == Double_Type, not real*4 == Float_Type       !!!

   variable ne = 299, ifl=1, phot = Double_Type[ne], energ = @phot;
   variable ebins = [-11.00:10.0:0.07];
   
   ebins = 10.^ebins;

   % Model called on a fixed grid, then interpolated

   jetmain(ebins,ne,pars,energ,phot);

   % Convert to keV and ph/cm^2/s/keV

   energ = (10.^energ)*(6.6260755e-18/1.60217733);
   phot = (10.^phot)*(10./6.6260755/energ); 
	
   % Interpolate back to original grid size
   jetinterp(ear,energ,phot,photar,ne,nein);

   % Multiply by normalization, and return in wavelength ascending order

   return reverse(par[0]*photar);
}
%change parameters here to not have corona and have inverted spectrum
static variable jet_input = { 
"norm",              1.,     0.,    10.,
"mbh [msun]",      1.e9,     3.,  3.e10,
"incl [deg]",       20.,     2.,    80.,
"dkpc [kpc]",     16.e3,     1.,   1.e6,
"redshift",			 0.,	 0.,	 5.,
"jetrat [L_edd]", 5.e-3,  1.e-7,     1.,
"r0 [r_g]",         20.,     2.,    80.,
"hratio",	     	 2.,     1.,     5.,		
"zdiss [r_g]",    1000.,    50.,   5.e5,
"zacc [r_g]",     1000.,   3.e2,   5.e5,
"zmax [r_g]",     1.e7,   1.e5,    1.e9,
"Te [Kev]",        100.,    10.,   2.e3,
"plfrac",           0.1,   0.05,   0.95,
"pldist",            0.,     0.,   	10.,
"pspec",             2.,    1.5,     3.,
"heat",              1.,     1.,    50.,
"betaeff",	     	0.1,  0.001,     1.,
"fsc",           1.5e-7,   1e-9,    0.1,
"pbeta",          2.e-3,   1e-3,     1.,
"sigf",		   	   1e-1,  1.e-2,     1.,
"Tin [keV/-L_edd]",1.e-2,  1.e-6,    5.,
"rin [r_g]",         1.,     1.,   2.e2,
"rout [r_g]",     1000.,    10.,   1.e5,
"compar1",         1.e3,     3.,   1.e6,
"compar2",           0.,     0.,     1.,
"compar3",		 3.e-10,     0.,     1.,
"compsw",			 0.,     0.,	 3.,
"velsw", 	         1.,     1.,    25.,
"infosw",	         0.,     0.,     5.,	
};

variable npar=length(jet_input)/4;
variable jet_pars=String_Type[npar];
variable jet_def=Float_Type[npar];
variable jet_min=Float_Type[npar];
variable jet_max=Float_Type[npar];

variable i;
for(i=0;i<=npar-1;i++)
{
    jet_pars[i]=jet_input[4*i];
    jet_def[i]=jet_input[4*i+1];
    jet_min[i]=jet_input[4*i+2];
    jet_max[i]=jet_input[4*i+3];
}

static variable jet_frz = Integer_Type[npar];

% Default above is no frozen parameters, now set the frozen ones

jet_frz[[0,1,2,3,4,4,9,10,12,13,14,15,16,19,22,23,24,25,26,27,28]]=1;

add_slang_function("bhjet",jet_pars);

define bhjet_defaults(i)
{
   return (jet_def[i],jet_frz[i],jet_min[i],jet_max[i]);
}

set_param_default_hook("bhjet","bhjet_defaults");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Just making some wavelength grids and parameters to test the output
% of the above

variable par = @jet_def;
variable pars = par[[1:28]];
variable lo = _A([-11.0:10.0:0.07]);
variable hi = make_hi_grid(lo);
