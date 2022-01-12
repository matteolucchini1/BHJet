//
// The XSPEC interface is slightly different to bhjet so we need a
// wrapper routine. This is also used by the Sherpa interface.
//
// Note that this has not yet been tested with XSPEC!
//

#include <cmath>

// The bhjet interface, but we use the interpolation scheme, as used
// by ISIS, and hard-code the grid the model is evaluated on here.
// This is not quite as flexible as in ISIS (since users can tweak
// the underlying grid the model is evaluated on) but it's not obvious
// how to do better for XSPEC.
//
extern void jetmain(double *ear,int ne,double *param,double *photeng,double *photspec);
extern void jetinterp(double *ear,double *energ,double *phot,double *photar,int ne,int newne);

// Match the grid used in bhjet.sl - this C++ code could take advantage of
// std library facitlities but has been kept pretty basic for now.
//
// I forget my S-Lang but
//     -11 + 0.07 * 298 =  9.86
//     -11 + 0.07 * 299 =  9.93
//     -11 + 0.07 * 298 = 10.00
// so I'm not quite sure what grid is meant to be being used.
//
const double LOG10_EMIN = -11.0;
const double LOG10_EMAX = 10.0;
const double LOG10_ESTEP = 0.07;
const int NBINS = 299;

#include <iostream>
using namespace std;

extern "C" {

  void xspec_jetinterp(const double *energy, int Nflux, const double *parameter,
		       int spectrum, double *flux, double *fluxError,
		       const char* init) {

    std::cout << "# DBG: user energy grid .." << std::endl;
    for (int i = 0; i < std::min(Nflux, 5); i++) {
      std::cout << "  energy[" << i << "]" << " = " << energy[i] << std::endl;
    }
    if (Nflux > 5) {
      std::cout << " ... " << std::endl;
      std::cout << "  energy[" << (Nflux - 1) << "]" << " = " << energy[Nflux - 1] << std::endl;
      std::cout << "  energy[" << Nflux << "]" << " = " << energy[Nflux] << std::endl;
    }

    // Can the definition of egrid be moved to compile time?
    //
    double *egrid = new double [NBINS];
    double *energ = new double [NBINS];
    double *phot = new double [NBINS];

    std::cout << "# DBG: egrid for jetmain .." << std::endl;
    double ebin = LOG10_EMIN;
    for (int i = 0; i < NBINS; i++) {
      egrid[i] = std::pow(10, ebin);
      ebin += LOG10_ESTEP;
      if ((i < 5) || (i > NBINS-5))
	std::cout << "  egrid[" << i << "] = " << egrid[i] << std::endl;
      if (i == 5)
	std::cout << " .. " << std::endl;
    }

    // Evaluate the model on the fixed grid
    jetmain(egrid, NBINS, const_cast<double *>(parameter), energ, phot);

    std::cout << "# Have called jetmain" << std::endl;

    // Interpolate onto the user grid after converting the values
    // to match XSPEC requirements.
    //
    for (int i = 0; i < NBINS; i++) {
      energ[i] = std::pow(10, energ[i]) * (6.6260755e-18 / 1.60217733);
      phot[i] = std::pow(10, phot[i]) * (10. / 6.6260755 / energ[i]);
    }

    jetinterp(const_cast<double *>(energy), energ, phot, flux, NBINS, Nflux);

    std::cout << "# Have called jetinterp" << std::endl;

    delete[] phot;
    delete[] energ;
    delete[] egrid;

  }

}
