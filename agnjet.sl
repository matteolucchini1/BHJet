
% Create a S-lang library of the function using SLIRP:

% unix%> slirp -make mffjet.f
% unix%> make 

% Move mffjet-module.so to your S-lang modules directory (for me,
% /usr/local/isis/modules).  In your ISIS .isisrc, or other S-lang
% input file, run all the code below.

import("agnjet");

define agnjet_fit(lo,hi,par)
{
   variable ear = [_A(hi),_A(lo[0])];
   variable nein = length(lo), photar=Double_Type[nein], photer=@photar;
   
   variable pars = par[[1:29]];

   % !!!  It is important to match the type for arrays that are input, !!!
   % !!!  and recieve return values.  I.e., phot & energ go in and     !!!
   % !!!  out as real*8 == Double_Type, not real*4 == Float_Type       !!!

   variable ne = 300, ifl=1, phot = Double_Type[ne], energ = @phot;
   variable ebins = [-11.00:10.0:0.07];
   
   ebins = 10.^ebins;

   % Model called on a fixed grid, then interpolated

   bhjet(ebins,ne,pars,energ,phot);

   % Convert to keV and ph/cm^2/s/keV

   energ = (10.^energ)*(6.6260755e-18/1.60217733);
   phot = (10.^phot)*(10./6.6260755/energ); 

   % Interpolate back to original grid size
   bhinterp(ear, energ, phot, photar, ne, nein);

   % Multiply by normalization, and return in wavelength ascending order

   return reverse(par[0]*photar);
}

static variable jet_input = { 
"norm",              1.,     0.,    10.,
"mbh [msun]",      1.e9,     3.,  3.e10,
"incl [deg]",       20.,     2.,    80.,
"dkpc [kpc]",     16.e3,     1.,   1.e6,
"jetrat [L_edd]", 5.e-3,  1.e-7,     1.,
"r0 [r_g]",         20.,     2.,    80.,
"hratio",	     	 2.,     1.,     5.,		
"zsh [r_g]",      1000.,    50.,   5.e5,
"zacc [r_g]",     1000.,   3.e2,   5.e5,
"zmax [r_g]",     1.e6,   1.e5,    1.e8,
"eltemp [gamma]",    3.,     1.,    20.,
"pspec",             2.,    1.5,     3.,
"heat",              1.,     1.,    50.,
"brk",		     	 1.,     1.,   100.,
"fsc",           1.5e-7,   1e-9,    0.1,
"gamfac",           10.,     1.,   1.e4,
"betapl",         2.e-3,   1e-3,     1.,
"sig",		   	   1e-2,  1.e-3,     1.,
"Tin [keV/-L_edd]",1.e-2,  1.e-6,    5.,
"rin [r_g]",         1.,     1.,   2.e2,
"rout [r_g]",     1000.,    10.,   1.e5,
"compar1",         1.e3,     3.,   1.e6,
"compar2",        1.e-9, 1.e-12,  1.e-1,
"compar3",		   3.e5,   3.e1,   3.e7,
"plfrac",           0.1,   0.05,   0.95,
"mxsw",              1.,     0.,     1.,
"velsw", 	         1.,     1.,    25.,
"compsw",			 0.,     0.,	 1.,
"plotsw",            1.,     0.,     1.,
"infosw",	         0.,     0.,     1.,	
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

jet_frz[[0,1,2,3,6,15,16,20,21,22,23,24,25,26,27,28,29]]=1;

add_slang_function("agnjet",jet_pars);

define agnjet_defaults(i)
{
   return (jet_def[i],jet_frz[i],jet_min[i],jet_max[i]);
}

set_param_default_hook("agnjet","agnjet_defaults");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Just making some wavelength grids and parameters to test the output
% of the above

variable par = @jet_def;
variable pars = par[[1:29]];
variable lo = _A([-11.0:10.0:0.07]);
variable hi = make_hi_grid(lo);
