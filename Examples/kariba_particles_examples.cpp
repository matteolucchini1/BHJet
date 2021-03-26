#include "kariba_examples.hh"

//This example compares all the non-thermal particle distrubtions in Kariba, for a common set of input parameters. 

int main(){
    
    int nel = 200;                                  //array size for particle distributions

    //Set the radius of the emitting regions in cm, the electron temperature, the maximum Lorenz factor, magnetic field in the
    //emitting region, effective expansion speed to define the adiabatic timescale, and non-thermal index. These parameter values
    //roughly mimic the corona of a BHB
    double R = 3.e7;
    double Te = 511.;
    double gmax = 1.e3;
    double bfield = 1.e5;
    double beta_exp = 0.1;
    double s = 2.;
    
    //These calls remove the output of previous runs from the output files
    clean_file("Output/lowfrac.dat",0);
    clean_file("Output/highfrac.dat",0);
    clean_file("Output/kdist.dat",0);
    clean_file("Output/bkndist.dat",0);
    clean_file("Output/lowfrac_cool.dat",0);
    clean_file("Output/highfrac_cool.dat",0);
    clean_file("Output/kdist_cool.dat",0);
    clean_file("Output/bkndist_cool.dat",0);
    
    //Declare all the objects for the different particle distributions; the only parameter needed is the size of the arrays to be
    //insantiated
    Mixed lowfrac(nel);
    Mixed highfrac(nel);
    Kappa kdist(nel);
    Bknpower bkndist(nel);
    Thermal placeholder(nel);
    
    //For the mixed particle distribution, you need to set the temperature of the thermal poo before setting up the Lorenz factor
    //and momentum arrays with set_p. Here it is called to a fixed maximum Lorenz factor. 
    //After that, specify the fraction of particles in non-thermal particles, the index of the non-thermal particles, and then
    //calculate the normalisation of the distributions for a given number density, which here is just set to unity.
    //The full distribution before cooling can only be calculated correctly with set_ndens() after all these setters have been
    //called correctly
    //Finally, after all this is done, the steady-state distribution including cooling can be computed by passing the energy 
    //density of external photons (here set to 0), the total number density of the final distribution, the magnetic field, the
    //radius of the emitting region and its expansion speed
    lowfrac.set_temp_kev(Te);
    lowfrac.set_p(gmax);
    lowfrac.set_plfrac(0.1);
    lowfrac.set_pspec(s);
    lowfrac.set_norm(1.);
    lowfrac.set_ndens();
    plot_write(nel,lowfrac.get_p(),lowfrac.get_gamma(),lowfrac.get_pdens(),lowfrac.get_gdens(),"Output/lowfrac.dat");
    lowfrac.cooling_steadystate(0.,1.,bfield,R,beta_exp);
    plot_write(nel,lowfrac.get_p(),lowfrac.get_gamma(),lowfrac.get_pdens(),lowfrac.get_gdens(),"Output/lowfrac_cool.dat");
    
    highfrac.set_temp_kev(Te);
    highfrac.set_p(gmax);
    highfrac.set_plfrac(0.9);
    highfrac.set_pspec(s);
    highfrac.set_norm(1.);
    highfrac.set_ndens();
    plot_write(nel,highfrac.get_p(),highfrac.get_gamma(),highfrac.get_pdens(),highfrac.get_gdens(),"Output/highfrac.dat");
    highfrac.cooling_steadystate(0.,1.,bfield,R,beta_exp);
    plot_write(nel,highfrac.get_p(),highfrac.get_gamma(),highfrac.get_pdens(),highfrac.get_gdens(),"Output/highfrac_cool.dat");
    
    //The kappa distribution is identical to the mixed distribution, except a) you don't need to specify the non-thermal fraction
    //and b) by definition the kappa index needs to be increased by 1 to match the same slope
    kdist.set_temp_kev(Te);
    kdist.set_p(gmax);
    kdist.set_kappa(s+1.);
    kdist.set_norm(1.);
    kdist.set_ndens();
    plot_write(nel,kdist.get_p(),kdist.get_gamma(),kdist.get_pdens(),kdist.get_gdens(),"Output/kdist.dat");
    kdist.cooling_steadystate(0.,1.,bfield,R,beta_exp);
    plot_write(nel,kdist.get_p(),kdist.get_gamma(),kdist.get_pdens(),kdist.get_gdens(),"Output/kdist_cool.dat");
    
    //For the broken power-law distribution, first instatiate a dummy thermal distribution. For the thermal distribution, before
    //calculating the momentum and Lorenz factor arrays you need to set the temperature, and before calculating the number density
    //arrays you need to specify the total number density. After that, it is possible to calculate the average Thermal momentum to
    //set the break momentum in the broken power-law.
    placeholder.set_temp_kev(Te);
    placeholder.set_p();
    placeholder.set_norm(1.);
    placeholder.set_ndens();
    double pbrk = placeholder.av_p();
    
    //For the broken power-law, first set up the slopes before and after the break, then set the momentum and Lorenz factor arrays
    //setting the minimum of the distribution to be one tenth of the thermal break, and again hard-coding the maximum.
    //After that, as usual, you need to specify the number density before calculating the number density of the uncooled 
    //distribution, and you need to pass photon energy density/particle number density/magnetic field/radius/expansion speed for 
    //the cooled distribution
    bkndist.set_pspec1(-2.);
    bkndist.set_pspec2(s);
    bkndist.set_p(0.1*pbrk,pbrk,gmax);
    bkndist.set_norm(1.);
    bkndist.set_ndens();
    plot_write(nel,bkndist.get_p(),bkndist.get_gamma(),bkndist.get_pdens(),bkndist.get_gdens(),"Output/bkndist.dat");
    bkndist.cooling_steadystate(0.,1.,bfield,R,beta_exp);
    plot_write(nel,bkndist.get_p(),bkndist.get_gamma(),bkndist.get_pdens(),bkndist.get_gdens(),"Output/bkndist_cool.dat");

    system("python3 Particles.py");

	return 0;
}
