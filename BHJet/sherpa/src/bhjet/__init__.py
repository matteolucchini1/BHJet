"""Sherpa interface to the BHJet model."""

from sherpa.models.parameter import Parameter, hugeval
from sherpa.astro.xspec import XSAdditiveModel, XSParameter

from . import _models


class BHJet(XSAdditiveModel):
    """BHJet model from Lucchini et al. 2021, arxiv: 2108.12011

    See [1]_.

    Parameters
    ----------
    mbh
        Mass of the Black Hole.
    incl
        Viewing angle of the jet.
    dkpc
        Distance to the source.
    redshift
        Only used for modelling the AGN.
    jetrat
        Ampunt of power injected at the base of the jet.
    r_0
        Rasius of the nozzle/corona.
    z_diss
        Location of the non-thermal particle injection region.
    z_acc
        Sets jet speed and dependency of magnetic field or the
        location of the jet accelaration. See [1]_.
    z_max
        Maximum length over which the jet calculations are made.
    t_e
        Temperature of relativistic electrons in the nozzle/corona.
    f_nth
        Fraction of thernal particles accelerated into the power-law
        tail.
    f_pl
        Reduces particle temperature and percentage of accelerated
        particles along the jet after z_diss.
    pspec
        Slope of non-thermal particle distribution.
    f_heat
        Imitates shock heating. Leave at 1 unless necessary.
    f_b
        The effective adiabatic cooling timescale.
    f_sc
        Complicated - see [1]_.
    p_beta
        The plasma beta at the base. See [1]_.
    sig_acc
        Sets the value for magnetization when velsw > 1.
    l_disk
        The luminosity of the Shakura-Sunyaev disk.
    r_in
        The inner radius of the disk.
    r_out
        The outer radius of the disk.
    compar1
        The first parameter of the external photon field. See compsw.
    compar2
        The second parameter of the external photon field. See compsw.
    compar3
        The third parameter of the external photon field. See compsw.
    compsw
        Sets the external photon fields to be  used in the IC calculation.
        See [1]_.
    velsw
        The jet-velocity profile. See [1]_.
    infosw
        What information to create along with the model. See [1]_.
    norm
        Unlike most XSPEC-like models, the normalization is frozen by default.

    References
    ----------

    .. [1] https://github.com/matteolucchini1/BHJet

    """
    _calc = _models.xspec_jetinterp

    def __init__(self, name='bhjet'):
        self.mbh = XSParameter(name, 'mbh', 1000000.0, min=3.0, max=30000000000.0, hard_min=3.0, hard_max=30000000000.0, frozen=True, units='msun')
        self.incl = XSParameter(name, 'incl', 40.0, min=2.0, max=80.0, hard_min=2.0, hard_max=80.0, frozen=True, units='deg')
        self.dkpc = XSParameter(name, 'dkpc', 3000.0, min=1.0, max=1000000.0, hard_min=1.0, hard_max=1000000.0, frozen=True, units='kpc')
        self.redshift = XSParameter(name, 'redshift', 0.0, min=0.0, max=7.0, hard_min=0.0, hard_max=7.0, frozen=True)
        self.jetrat = XSParameter(name, 'jetrat', 0.005, min=1e-07, max=1.0, hard_min=1e-07, hard_max=1.0, units='L_edd')
        self.r_0 = XSParameter(name, 'r_0', 20.0, min=2.0, max=80.0, hard_min=2.0, hard_max=80.0, units='r_g')
        self.z_diss = XSParameter(name, 'z_diss', 1000.0, min=50.0, max=500000.0, hard_min=50.0, hard_max=500000.0, units='r_g')
        self.z_acc = XSParameter(name, 'z_acc', 1000.0, min=300.0, max=500000.0, hard_min=300.0, hard_max=500000.0, frozen=True, units='r_g')
        self.z_max = XSParameter(name, 'z_max', 10000000.0, min=100000.0, max=1000000000.0, hard_min=100000.0, hard_max=1000000000.0, frozen=True, units='r_g')
        self.t_e = XSParameter(name, 't_e', 100.0, min=10.0, max=2000.0, hard_min=10.0, hard_max=2000.0, units='keV')
        self.f_nth = XSParameter(name, 'f_nth', 0.1, min=0.05, max=0.95, hard_min=0.05, hard_max=0.95, frozen=True)
        self.f_pl = XSParameter(name, 'f_pl', 0.0, min=0.0, max=10.0, hard_min=0.0, hard_max=10.0, frozen=True)
        self.pspec = XSParameter(name, 'pspec', 2.0, min=1.5, max=3.0, hard_min=1.5, hard_max=3.0)
        self.f_heat = XSParameter(name, 'f_heat', 1.0, min=1.0, max=50.0, hard_min=1.0, hard_max=50.0, frozen=True)
        self.f_b = XSParameter(name, 'f_b', 0.1, min=0.001, max=1.0, hard_min=0.001, hard_max=1.0, frozen=True)
        self.f_sc = XSParameter(name, 'f_sc', 1.5e-07, min=1e-09, max=0.1, hard_min=1e-09, hard_max=0.1)
        self.p_beta = XSParameter(name, 'p_beta', 0.05, min=0.001, max=1.0, hard_min=0.001, hard_max=1.0)
        self.sig_acc = XSParameter(name, 'sig_acc', 0.1, min=0.01, max=1.0, hard_min=0.01, hard_max=1.0, frozen=True)
        self.l_disk = XSParameter(name, 'l_disk', 0.01, min=0.0001, max=1.0, hard_min=0.0001, hard_max=1.0, units='L_edd')
        self.r_in = XSParameter(name, 'r_in', 1.0, min=1.0, max=200.0, hard_min=1.0, hard_max=200.0, units='r_g')
        self.r_out = XSParameter(name, 'r_out', 1000.0, min=10.0, max=100000.0, hard_min=10.0, hard_max=100000.0, frozen=True, units='r_g')
        self.compar1 = XSParameter(name, 'compar1', 1000.0, min=3.0, max=1000000.0, hard_min=3.0, hard_max=1000000.0, frozen=True)
        self.compar2 = XSParameter(name, 'compar2', 0.0, min=0.0, max=1.0, hard_min=0.0, hard_max=1.0, frozen=True)
        self.compar3 = XSParameter(name, 'compar3', 3e-10, min=0.0, max=1.0, hard_min=0.0, hard_max=1.0, frozen=True)
        self.compsw = XSParameter(name, 'compsw', 0, alwaysfrozen=True)
        self.velsw = XSParameter(name, 'velsw', 1.0, min=1.0, max=25.0, hard_min=1.0, hard_max=25.0, frozen=True)
        self.infosw = XSParameter(name, 'infosw', 0, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, min=0.0, max=1e+24, hard_min=0.0, hard_max=1e+24, frozen=True)
        XSAdditiveModel.__init__(self, name, (self.mbh,self.incl,self.dkpc,self.redshift,self.jetrat,self.r_0,self.z_diss,self.z_acc,self.z_max,self.t_e,self.f_nth,self.f_pl,self.pspec,self.f_heat,self.f_b,self.f_sc,self.p_beta,self.sig_acc,self.l_disk,self.r_in,self.r_out,self.compar1,self.compar2,self.compar3,self.compsw,self.velsw,self.infosw,self.norm))
