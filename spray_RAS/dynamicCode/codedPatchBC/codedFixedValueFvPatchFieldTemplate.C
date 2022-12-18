/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "codedFixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = bc14de22b75262eefd017eccf6696ae94d9539bc
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void codedPatchBC_bc14de22b75262eefd017eccf6696ae94d9539bc(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchScalarField,
    codedPatchBCFixedValueFvPatchScalarField
);


const char* const codedPatchBCFixedValueFvPatchScalarField::SHA1sum =
    "bc14de22b75262eefd017eccf6696ae94d9539bc";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

codedPatchBCFixedValueFvPatchScalarField::
codedPatchBCFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF)
{
    if (false)
    {
        Info<<"construct codedPatchBC sha1: bc14de22b75262eefd017eccf6696ae94d9539bc"
            " from patch/DimensionedField\n";
    }
}


codedPatchBCFixedValueFvPatchScalarField::
codedPatchBCFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct codedPatchBC sha1: bc14de22b75262eefd017eccf6696ae94d9539bc"
            " from patch/dictionary\n";
    }
}


codedPatchBCFixedValueFvPatchScalarField::
codedPatchBCFixedValueFvPatchScalarField
(
    const codedPatchBCFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct codedPatchBC sha1: bc14de22b75262eefd017eccf6696ae94d9539bc"
            " from patch/DimensionedField/mapper\n";
    }
}


codedPatchBCFixedValueFvPatchScalarField::
codedPatchBCFixedValueFvPatchScalarField
(
    const codedPatchBCFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct codedPatchBC sha1: bc14de22b75262eefd017eccf6696ae94d9539bc "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

codedPatchBCFixedValueFvPatchScalarField::
~codedPatchBCFixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy codedPatchBC sha1: bc14de22b75262eefd017eccf6696ae94d9539bc\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void codedPatchBCFixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs codedPatchBC sha1: bc14de22b75262eefd017eccf6696ae94d9539bc\n";
    }

//{{{ begin code
    #line 30 "/home/andrea/spray_RAS3/0/tracer/boundaryField/in_air"
const fvPatch& patch = this->patch();
            const vectorField& cf = patch.Cf();
            
            scalarField& field = *this;
            forAll(cf, i)
            {
                field[i] = ((pow(cf[i].y(),2) + pow(cf[i].x(),2)) <= 0.75) ? 1. : 0.;
            }
//}}} end code

    this->fixedValueFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

