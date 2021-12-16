# BHJET

BHJet is a semi-analytical, multi-zone jet model designed for modelling steady-state SEDs of jets launched from accreting black holes. The key features of the model are:
1) It is applicable across the entire black hole mass scale, from black hole X-ray binaries (both low and high mass) to active galactic nuclei of any class (from low-luminosity AGN to flat spectrum radio quasars),
2) It is designed to be more comparable than other codes to GRMHD simulations and/or RMHD semi-analytical solutions.

The model is fairly complex and fitting it to data is not always straightforward. As such, it is highly recommended to read this file carefully before running the code. It takes little time and will save you a lot of headaches later on. All the physics of the model, the assumptions going into it, as well as some applications, are discussed in depth in Lucchini et al. 2021, arxiv: [2108.12011](https://arxiv.org/abs/2108.12011), DOI: WIP , and references therein.

If you have any questions, or find issues or bugs, please feel free to open a ticket in the Issues section of the repository.

---------------------------------------------------------------------------------------------------------------------------------------

## Code structure

The BHJET model, contained in the /BHJet/ subfolder, is built on a simple C++ library called Kariba. The library is contained in the /Kariba/ folder; /Examples/ includes a set of simple codes to highlight features of Kariba and to replicate the plots included in Lucchini et al, 2021. The classes in Kariba handle the particle distribution and radiation calculations through their methods. More documentation is included in the repository sub-folders.

The codes in /BHJet/ and /Examples/ can be compiled with the ./MakeBHJet and ./MakeCode commands in their respective folders, in standard C++ fashion. Python plotting scripts are also included and automatically called to visualize the output. Furthermore, /BHJet/ can be imported in spectral fitting packages ISIS and Sherpa (with experimental support for XSPEC), allowing users to combine it with the traditional Heasoft library of spectral models.

---------------------------------------------------------------------------------------------------------------------------------------

## Citing BHJet

Please cite (Lucchini et al. 2021a,b -- placeholder) if you find Kariba useful in your research. If you find BHJet useful, please also cite the original agnjet paper, Markoff et al. 2001, 2005. The bibtex entries are:

    @ARTICLE{2005ApJ...635.1203M,
           author = {{Markoff}, Sera and {Nowak}, Michael A. and {Wilms}, J{\"o}rn},
            title = "{Going with the Flow: Can the Base of Jets Subsume the Role of Compact Accretion Disk Coronae?}",
          journal = {\apj},
         keywords = {Accretion, Accretion Disks, Black Hole Physics, Radiation Mechanisms: Nonthermal, X-Rays: Binaries, X-Rays: General, Astrophysics},
             year = 2005,
            month = dec,
           volume = {635},
           number = {2},
            pages = {1203-1216},
              doi = {10.1086/497628},
    archivePrefix = {arXiv},
           eprint = {astro-ph/0509028},
     primaryClass = {astro-ph},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2005ApJ...635.1203M},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }


    @ARTICLE{2001A&A...372L..25M,
           author = {{Markoff}, S. and {Falcke}, H. and {Fender}, R.},
            title = "{A jet model for the broadband spectrum of XTE J1118+480. Synchrotron emission from radio to X-rays in the Low/Hard spectral state}",
          journal = {\aap},
         keywords = {X-RAYS: BINARIES, X-RAYS: INDIVIDUAL: XTE J1118+480, RADIATION MECHANISMS: NON-THERMAL, STARS: WINDS, OUTFLOWS -BLACK HOLE PHYSICS, ACCRETION, ACCRETION DISKS, Astrophysics},
             year = 2001,
            month = jun,
           volume = {372},
            pages = {L25-L28},
              doi = {10.1051/0004-6361:20010420},
    archivePrefix = {arXiv},
           eprint = {astro-ph/0010560},
     primaryClass = {astro-ph},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2001A&A...372L..25M},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
