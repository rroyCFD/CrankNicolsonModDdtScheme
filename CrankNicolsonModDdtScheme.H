/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Author:
- Rajib Roy
- University of Wyoming
- rroy@uwyo.edu, roy.rajib@live.com

License
    This file is compatible with OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::fv::CrankNicolsonModDdtScheme

Description
    Second-order CrankNicolsonMod-differencing implicit ddt using the current and two previous time-step values.

    The Crank-Nicolson scheme is often unstable for complex flows in complex
    geometries and it is necessary to "off-centre" the scheme (i.e. blend with Euler-implicit scheme) to stabilize it while retaining greater temporal accuracy than the first-order Euler-implicit scheme.  

    Off-centering is specified via the mandatory coefficient \c ocCoeff in the range [0,1] following the scheme name e.g.
    \verbatim
    ddtSchemes
    {
        default         CrankNicolsonMod 0.9;
    }
    \endverbatim
    
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
    
    Need Function1Mesh class to utilize the dammping functions.

SourceFiles
    CrankNicolsonModDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef CrankNicolsonModDdtScheme_H
#define CrankNicolsonModDdtScheme_H

#include "ddtScheme.H"
#include "Function1Mesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class CrankNicolsonModDdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CrankNicolsonModDdtScheme
:
    public fv::ddtScheme<Type>
{
    // Private Member Functions

        //- Return the current time-step
        scalar deltaT_() const;

        //- Return the previous time-step
        scalar deltaT0_() const;

        //- Return the previous time-step or the old timestep, if field
        //  wasn't available in which case Euler ddt is used
        template<class GeoField>
        scalar deltaT0_(const GeoField&) const;

        //- Off-centering coefficient function
        //  1 -> CN, less than one blends with EI
        autoPtr<Function1Mesh<scalar>> ocCoeff_;

        //- Disallow default bitwise copy construct
        CrankNicolsonModDdtScheme(const CrankNicolsonModDdtScheme&);

        //- Disallow default bitwise assignment
        void operator=(const CrankNicolsonModDdtScheme&);


public:

    //- Runtime type information
    TypeName("CrankNicolsonMod");


    // Constructors
        //- Construct from mesh
        CrankNicolsonModDdtScheme(const fvMesh& mesh);

        //- Construct from mesh and Istream
        CrankNicolsonModDdtScheme(const fvMesh& mesh, Istream& is);
        

    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return fv::ddtScheme<Type>::mesh();
        }

        //- Return the current off-centreing coefficient
        scalar ocCoeff() const
        {
            return ocCoeff_->value(mesh().time().value());
        }

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const dimensioned<Type>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );

        typedef typename ddtScheme<Type>::fluxFieldType fluxFieldType;

        tmp<fluxFieldType> fvcDdtUfCorr
        (
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
        );

        tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const fluxFieldType& phi
        );

        tmp<fluxFieldType> fvcDdtUfCorr
        (
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
        );

        tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const fluxFieldType& phi
        );

        tmp<surfaceScalarField> meshPhi
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );
};


template<>
tmp<surfaceScalarField> CrankNicolsonModDdtScheme<scalar>::fvcDdtUfCorr
(
    const GeometricField<scalar, fvPatchField, volMesh>& U,
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& Uf
);

template<>
tmp<surfaceScalarField> CrankNicolsonModDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& U,
    const surfaceScalarField& phi
);

template<>
tmp<surfaceScalarField> CrankNicolsonModDdtScheme<scalar>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& Uf
);

template<>
tmp<surfaceScalarField> CrankNicolsonModDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& phi
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CrankNicolsonModDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
