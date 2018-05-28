/* Boost libs/numeric/odeint/examples/simple1d.cpp

 Copyright 2012-2013 Mario Mulansky
 Copyright 2012 Karsten Ahnert

 example for a simple one-dimensional 1st order ODE

 Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#include <iostream>
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#define zeta 1.0

using namespace std;
using namespace boost::numeric::odeint;
// using namespace boost::multiprecision;

//typedef boost::multiprecision::cpp_dec_float_50 value_type;
//typedef boost::array< value_type , 1 > state_type;
//state_type = double;
typedef runge_kutta_dopri5< double > stepper_type;

//const int Steps = 4;
//typedef double value_type;
//typedef boost::array< double , 2 > state_type;
//typedef runge_kutta_fehlberg78<state_type> initializing_stepper_type;
//typedef adams_bashforth_moulton< Steps , state_type > stepper_type;

struct push_back_z_and_gb{
    std::vector< double >& tt;
    std::vector< double >& yy;
    
    push_back_z_and_gb(std::vector< double > &states,std::vector<double>&times):yy(states),tt(times){}
    
    void operator()(const double &y,double t){
        tt.push_back( t );
        yy.push_back( y );
    }
};


/* Eq (2) from Falcke 1996. What is used in AGNJET */
void rhs0(const double y,double &dydt,const double t){
    double Gam = 4.0/3.0;
    double xi = pow((y/(Gam*(Gam-1)/(Gam+1))),(1-Gam));
    dydt = 2/t;
//    cout << dydt << endl;
    dydt *= ( y/( (Gam+xi)/(Gam-1)*pow(y,2)-Gam));
//    cout << t << " " << y << " " << dydt << " " << xi << endl;
}

/* Eq (9) & (10) from my writeup. Basically what you would get 
 * from the assumptions in Falcke 1996 minus the algebra mistakes i.e.
 * quasi-isothermal, d/dz(U/n) = 0
 */
void rhs1(const double y,double &dydt,const double t){
    double Gam = 4.0/3.0;
    double delta;
    delta = pow((y / sqrt(Gam*(Gam-1)/(1+2*Gam-Gam*Gam))),(Gam-1));
    dydt = 2/t;
    dydt *= y/( (Gam+delta)/(Gam-1)*y*y-Gam);
}
/* Full Euler assuming quasi-isothermal like Falcke 1996 but no longer neglecting d/dz(U/n) */
void rhs2(const double y,double &dydt,const double t){
    double Gam = 4.0/3.0;
    double alpha, beta;
//    double zeta;
    double y0;
    alpha = Gam;
    beta = 1.0;
//    zeta = 10.0;
//    y0 = sqrt(Gam*(Gam-1.0)/(1.0+2.0*Gam-Gam*Gam));
    y0 = sqrt(zeta*Gam*(Gam-1.0)/(1.0+2.0*Gam*zeta-Gam*Gam*zeta));
    
    dydt = pow((Gam*y/(Gam-1.0) * (1./(Gam*zeta) *pow((y/y0),(-beta+alpha)) + beta - alpha + 1.0) - alpha/y),-1.0);
    
//    cout << y << " " << y/(Gam-1) << endl;
//    cout << (Gam*(y/(Gam-1.))) * ((1./Gam/zeta) *pow((y/y0),-beta+alpha)+beta-alpha+1.) << endl;
//    cout << alpha/y << endl;
    
//    cout << dydt << endl;
    dydt *= 2./t;
//    cout << dydt << endl;
    
//    exit(0);
}

/* Full Euler equation assuming adiabatic jet. For whatever reason
 * I went with the Bernoulli equation instead of the Euler equation 
 * here but Euler eq would give the same result. 
 */
void rhs3(const double y,double &dydt,const double t){
    double Gam = 4.0/3.0;
    double t0, n_normed,xi;
    double y0, beta_s;
//    double zeta;
    
    t0 = 1.0;
//    zeta = 10.0;
    beta_s = sqrt(zeta*Gam*(Gam-1)/(zeta*Gam+1)); // initial sound speed assuming U_0 = n_0 m_p c^2
    y0 = 1/sqrt(pow(beta_s,-2)-1); //\gamma_0\beta_0 if  n_0 m_p c^2
    
    n_normed = pow((y/y0),-1)*pow((t/t0),-2);
    xi = (1+zeta*Gam)*(Gam-1)*zeta*Gam*pow(n_normed,(Gam-1));
    xi *= pow((1+zeta*Gam*pow(n_normed,(Gam-1))),-2);
    xi *= sqrt(y0*y0 +1);
    dydt = y/sqrt(y*y+1);
    dydt -= xi/y;
    dydt = pow(dydt,-1);
    dydt *=xi*2/t;
    
//    cout << dydt << endl;
//    exit(0);
}

void write_cout(const double &y,double t){
    cout << t << '\t' << y << endl;
}

//int main()
void VelProf(int velsw,double Gam,double gbx[],double gby[],int &sizegb){
    /* 
     * velsw = 0 AGNJET;
     * velsw = 1 CORRECTED AGNJET;
     * velsw = 2 full 1D quasi-isothermal Bernoulli Eq.;
     * velsw = 3 full 1D adiabatic Bernoulli Eq.;
     */
//    int velsw = 3;
//    cout.precision( 20 );
    double beta_s,y0;
//    double zeta = 1.;
//    Gam = 4./3.;
    beta_s = sqrt(zeta*Gam*(Gam-1.)/(zeta*Gam+1.)); // initial sound speed assuming U_0 = n_0 m_p c^2
    y0 = 1./sqrt(pow(beta_s,-2.)-1.); // \gamma_0\beta_0 if  n_0 m_p c^2
    
    size_t steps;
    vector<double> yy;
    vector<double> tt;
    
    double y = y0;
    
    if (velsw == 0){
        steps = integrate_adaptive(make_controlled(1E-3,1E-3,stepper_type()),rhs0,y,pow(10,0.0),pow(10,20.0),
        pow(10,0.2),push_back_z_and_gb(yy,tt));
    } else if (velsw == 1){
        steps = integrate_adaptive(make_controlled(1E-3,1E-3,stepper_type()),rhs1,y,pow(10,0.0),pow(10,20.0),
        pow(10,0.2),push_back_z_and_gb(yy,tt));
    } else if (velsw == 2){
        steps = integrate_adaptive(make_controlled(1E-3,1E-3,stepper_type()),rhs2,y,pow(10,0.0),pow(10,20.0),
        pow(10,0.2),push_back_z_and_gb(yy,tt));
    } else {
        steps = integrate_adaptive(make_controlled(1E-3,1E-3,stepper_type()),rhs3,y,pow(10,0.0),pow(10,20.0),
        pow(10,0.2),push_back_z_and_gb(yy,tt));
    }
    /* output */
    sizegb = steps+2;
//    double gbx[sizegb],gby[sizegb];
    for( size_t i=0; i<=steps+1; i++ )
    {
        gbx[i] = tt[i];
        gby[i] = yy[i];
        if (i==steps+1) {
            gbx[i] = steps+2;
            gby[i] = steps+2;
        }
//        cout << gbx[i] << '\t' << gby[i] << '\n';
    }
    tt.clear();
    yy.clear();
}
