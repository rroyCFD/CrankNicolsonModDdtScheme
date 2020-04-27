### Author:
- Rajib Roy
- University of Wyoming
- rroy@uwyo.edu, roy.rajib@live.com

### Preamble:

The OpenFOAM-v6 Crank-Nicolson scheme has a potential bug, which described in the cfd-online forum post [Crank Nicolson scheme implemented wrong?](https://www.cfd-online.com/Forums/openfoam-programming-development/226125-crank-nicolson-scheme-implemented-wrong.html)

The major bug included an inaccurate calculation of the temporal discretization coefficients. This implementation attempts to fix the coefficients. Another minor bug was to use previous temporal derivative instead of old-old time value.

The scheme is implemented following the explanation in the book by Moukalled et al. (section 13.3.5, page 512) 

	Moukalled, F., Mangani, L., & Darwish, M. (2015). 
	The Finite Volume Method in Computational Fluid Dynamics: An Advanced Introduction with OpenFOAM® and Matlab.
	https://doi.org/10.1007/978-3-319-16874-6

### Description:

Second-order CrankNicolsonMod-differencing implicit ddt using the current and two previous time-step values. 

The Crank-Nicolson scheme is often unstable for complex flows in complex
geometries and it is necessary to "off-centre" the scheme (i.e. blend with Euler-implicit scheme) to stabilize it while retaining greater temporal accuracy than the first-order Euler-implicit scheme. Off-centering is specified via the mandatory coefficient \c ocCoeff in the range [0,1] following the scheme name e.g.

    ddtSchemes
    {
        default         CrankNicolsonMod 0.9;
    }

    
or with an optional "damp" function to dynamically blend (basded on maximum Courant number (maxCo)) from the pure/blended Crank-Nicolson schme towards Euler-implicit scheme. Blending begins when the maxCo exceeds the threshold value and finaly drives the ocCoeff when maxCo is equal or higher than maximum permitted Corant number (maxValue). The nature of blend parmeter depends on the damping function (linear, quadratic, cosine, exponential, etc.), e.g.

    ocCoeff
    {
        type            scale;

        scale
        {
            type            linearDampCFL; 
            // type            exp10DampCFL; 
            // type            quadraticDampCFL;
            // type            halfCosineDampCFL;
            start           0;
            thresholdValue  0.50;
            maxValue        1.00;
        }

        value           1.0; // constValue
    };
    
Need **Function1Mesh** class to utilize the dammping functions. https://github.com/rroyCFD/Function1Mesh.git


### Make

#### files

	ddtSchemes = finiteVolume/ddtSchemes
	$(ddtSchemes)/CrankNicolsonModDdtScheme/CrankNicolsonModDdtSchemes.C

	LIB = $(FOAM_USER_LIBBIN)/libfiniteVolume_$(USER)

#### options

	EXE_INC = \
	    -I$(LIB_SRC)/triSurface/lnInclude \
	    -I$(LIB_SRC)/meshTools/lnInclude \
	    -I$(LIB_SRC)/finiteVolume/lnInclude \
	    -I$(WM_PROJECT_USER_DIR)/src/OpenFOAM/lnInclude \

	LIB_LIBS = \
	    -ltriSurface \
	    -lmeshTools \
	    -lfiniteVolume \
	    -L$(FOAM_USER_LIBBIN) \
	    -lOpenFOAM_$(USER)
