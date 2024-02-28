# BHJET

The free parameters of the model are:

 - mbh         : mass of the black hole. Always frozen before the fit, bhjet should NOT be used to estimate BH masses from SED modelling.
 - Incl        : viewing angle of the jet. Sets Doppler factor for the various regions. Typically frozen before the fit.
 - dkpc        : distance from source in kpc. Always frozen before the fit.
 - redshfit    : self explanatory. Only used for modelling AGN and set before the fit.
 - jetrat      : amount of power injected at the base of the jet. Increasing values increase the normalisation of the synchrotron flux, and the importance of Comptonisation. Measured in Eddington units. 
 - r\_0        : radius of the nozzle/corona, described by an outflowing, magnetized cylinder of radius r\_0 and height 2r\_0. Decreasing values increase the optical depth of the nozzle. Measured in units of r\_g.
 - z\_diss     : location of non-thermal particle injection region. Sets optically thick to thin break, overall normalization and Compton dominance of non-thermal component, and for high accretion rate AGN sets the EC target field (BLR or BLR+torus). Decreasing values move the synchrotron thick-to-thin break to higher frequencies. Measured in units of r\_g.
 - z\_acc      : if velsw > 1, sets location of the jet acceleration. Sets the jet speed and dependency of magnetic field with distance; smaller values correspond to faster dissipation before z\_acc (see velsw below). Measured in r\_g like zdiss. 
 - z\_max      : maximum length over which jet calculations are carried out, in units of Rg. Always frozen to large values in order to get a flat radio spectrum.
 - t\_e        : temperature of relativistic electrons in the nozzle/corona, expressed in keV.
 - f\_nth      : Fraction of thermal particles accelerated into power-law tail. Typical values range between 0.01 and 1. Sets the relative importance of thermal to non-thermal emission in the compact jet. Typically frozen before the fit.
 - f\_pl       : reduces particle temperature and percentage of accelerated particles along the jet after z\_diss, resulting in an inverted radio spectrum. Set to 0 for a standard flat spectrum, increase for a more inverted spectrum.
 - pspec       : slope of non-thermal particle distribution. Sets the slope of the non-thermal synchrotron and inverse Compton components.
 - f\_heat     : imitates shock heating; at z=z\_diss bumps Te by a fixed factor f\_heat, increasing values increase the radiative efficiency of the jet after z\_diss. Set to 1 and frozen unless necessary.
 - f\_b        : sets the effective adiabatic cooling timescale, defined as t_ad = r/f\_b*c. This in turn sets the location of the cooling break in the radiating particle distribution. Kept free if the SED shows a cooling break, otherwise set to 0.1 and frozen.
 - f\_sc       : If < 0.1, sets the maximum energy of non-thermal particles by parametrizing the acceleration timescale (and therefore acceleration efficiency). If > 10, sets the maximum Lorenz factor of the non-thermal particles. Either way, sets the energy up to which the non-thermal synchrotron and inverse Compton spectra extend to. Values between 0.1 and 10 are non-physical and should be disregarded. 
 - p\_beta     : plasma beta (ratio between lepton and magnetic field energy density) at the base. If velsw = 0 or 1 this sets (and freezes) the equipartition value throughout the jet as well as the pair content, if velsw > 1 it only sets the pair content and has to be frozen to a value high enough to have pair content ~approx unity. This is because the assumptions going into the calculation of the magnetic field when velsw>1 only hold if lepton energy density <<< (cold) proton energy density. Additionally, if velsw > 1, this can be set to 0 to always enforce one proton per electron. Can be used to tune the optical depth of the jet nozzle for a given jet power. _See the main paper for a discussion on valid values of p\_beta for all model flavours._
 - sig\_acc    : only used if velsw > 1, sets the value for magnetization. Higher values decrease the Compton dominance and increase the peak frequencies of the non-thermal synchrotron and inverse Compton components. Typically frozen to ~0.01-1 unless required by the data.
 - l\_disk     : sets the luminosity of the Shakura-Sunyaev disk in Eddington units, and the corresponding temperature is computed from knowing this luminosity and Rin. If set to a negative value, the code still includes the disk as a seed photon field for IC scattering, but it is not included in the output of the model. This is in case users want to use a more complex disk model (e.g. diskir). 
 - r\_in       : inner radius of the disk.
 - r\_out      : outer radius of the disk. Only has a minimal contribution to the SED, typically frozen. rout<=rin disables the disk.
 - compar1     : first parameter of the external photon field, depending on the value of compsw
 - compar2     : second parameter of the external photon field, depending on the value of compsw
 - compar3     : third parameter of the external photon field, depending on the value of compsw
 - compsw      : sets the external photon fields to be used in the IC calculation. 0 only includes SSC (+disk IC, if present) emission, 1 adds a uniform external black body contribution (e.g. a host galaxy), 2 adds the Broad Line Region and Torus of an AGN and ties their luminosity to that of the disk. If compsw = 1, compar1 is the temperature in Kelvin, compar2 is the total black-body luminosity, and compar3 is the BB energy density. If compsw = 2, compar1 is the fraction of disk photons reprocessed by the BLR, compar2 the fraction of disk photons reprocessed by the torus, and compar3 is not used.
 - velsw       : This sets the jet velocity profile used by the code. If velsw = 0 or 1 the jet is pressure-driven and mildly relativistic, otherwise the jet base is highly magnetized and the jet is magnetically driven. In this case, the value of velsw is also the final Lorentz factor of the jet, achieved at a distance z\_acc (see above).
 - infosw      : information switch; the higher the value, the more info the code returns. 0 simply stores a total flux value in the photspec array; 1 also prints the total emission to, as well as contribution from each radiative component (e.g. thermal synchrotron, non-thermal synchrotron, thermal Comptonization, non-thermal IC emission, external fields, disk) to a file; 2 also prints the emission and particle distribution in each zone, which will be automatically plotted by Plot.py; 3, 4 and 5 print increasing amounts of information to the terminal.
 - EBLsw	   : switch to activate extragalactic background light (EBL) attenuation. If EBLsw = 0 does not use EBL attenuation, and 1 uses the fiducial model by Gilmore et al. (2012).
 
Generally, no more than ~9 parameters should be fitted at the same time due to model degeneracies; depending on the application, many parameters can be frozen to a reasonable educated guess.

---------------------------------------------------------------------------------------------------------------------------------------

There are two important notes for running the code with acceptable performance:
1) The jet is arbitrarily divided into 100 zones (controlled by the integer variable nz) which do not interact with each other in any way. As the calculations move from the base of the jet to zmax, the grid which defines each zone changes. Up to 1000 Rg, the grid advances in steps of 2*r, where r is the radius of the jet at that particular distance. Outwards, the grid uses logarithmic steps instead. The reason for this choice (discussed in Connors et al. 2018) is that it prevents the final spectrum from being resolution-dependent. On top of that, the code runtime increases linearly with the zone count, so anything above ~70 will simply slow the code down for no gain. Note that particularly large values of zmax and/or small values of r0 can result in the number of zones being insufficient to cover the entire jet length; this is why nz is set slightly higher than the 70 segments necessary to optimize performance in most applications.
2) The bulk of the code runtime is caused by the inverse Compton calculation, in particular for cases of moderate to high optical depth (>~0.05, when multiple scatters are considered) and/or when many zones in the jet produce bright IC emission. The code uses two adjustments to improve the efficiency of the radiation calculations: a) the frequency grid over which the emission of each zone (both inverse Compton and synchrotron) is updated dynamically after calculating the appropriate scale frequencies, and b) the inverse Compton itself is only calculated when it's expected to be bright enough to contribute meaningfully to the SED through the Compton\_check function. Depending on the source and model parameters, this reduces the code run time by a factor of ~3-50. _It is extremely important to double check whether Compton\_check is being too aggressive in neglecting zones to compute IC or not! You can do this simply by forcing Compton\_check to always return true, thus calculating the IC emission for every zone, and comparing the result with the standard prescriptions in the function._

---------------------------------------------------------------------------------------------------------------------------------------

For running outside as a stand-alone code: 

Compile using:

./MakeBHJet

Input parameters can be changed in ip.dat in the Input directory. If you want to use a different file, change the name in jetwrap.cc and recompile. 

Run using:

./bhwrap.x

If you want to quickly plot the output, a simple python script is provided (Plot.py) and ran automatically by bhwrap.x. If you want to use your own plotting tools, this can be disabled by commenting out line 79 in the wrapper: 

system("python Plot.py");


---------------------------------------------------------------------------------------------------------------------------------------

For running inside ISIS:

./slirpAgnjet
make
make test

add in your .isisrc the following lines:
append_to_isis_load_path(path+"path/to/your/agnjet/folder");
append_to_isis_module_path(path+"path/to/your/agnjet/folder");

(on some systems if gsl libs aren't on obvious path you may need to
       explicitly include with -I and -L, but try this first!)

---------------------------------------------------------------------------------------------------------------------------------------
