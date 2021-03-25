#include "kariba_examples.hh"

int main(){
    
    int nel = 200;

    double R = 3.e7;
    double Te = 511.;
    double gmax = 1.e3;
    double bfield = 1.e5;
    double beta_exp = 0.1;
    double s = 2.;
    
    clean_file("Output/lowfrac.dat",0);
    clean_file("Output/highfrac.dat",0);
    clean_file("Output/kdist.dat",0);
    clean_file("Output/bkndist.dat",0);
    clean_file("Output/lowfrac_cool.dat",0);
    clean_file("Output/highfrac_cool.dat",0);
    clean_file("Output/kdist_cool.dat",0);
    clean_file("Output/bkndist_cool.dat",0);
    
    Mixed lowfrac(nel);
    Mixed highfrac(nel);
    Kappa kdist(nel);
    Bknpower bkndist(nel);
    Thermal placeholder(nel);
    
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
    
    kdist.set_temp_kev(Te);
    kdist.set_p(gmax);
    kdist.set_kappa(s+1.);
    kdist.set_norm(1.);
    kdist.set_ndens();
    plot_write(nel,kdist.get_p(),kdist.get_gamma(),kdist.get_pdens(),kdist.get_gdens(),"Output/kdist.dat");
    kdist.cooling_steadystate(0.,1.,bfield,R,beta_exp);
    plot_write(nel,kdist.get_p(),kdist.get_gamma(),kdist.get_pdens(),kdist.get_gdens(),"Output/kdist_cool.dat");
    
    placeholder.set_temp_kev(Te);
    placeholder.set_p();
    placeholder.set_norm(1.);
    placeholder.set_ndens();
    double pbrk = placeholder.av_p();
    
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
