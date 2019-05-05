/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "lognormal.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributionModels
{
    defineTypeNameAndDebug(lognormal, 0);
    addToRunTimeSelectionTable(distributionModel, lognormal, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionModels::lognormal::lognormal
(
    const dictionary& dict,
    Random& rndGen
)
:
    distributionModel(typeName, dict, rndGen),
    minValue_(readScalar(distributionModelDict_.lookup("minValue"))),
    maxValue_(readScalar(distributionModelDict_.lookup("maxValue"))),
    expectation_(m_to_mu(readScalar(distributionModelDict_.lookup("expectation")),
            readScalar(distributionModelDict_.lookup("variance")))),
    variance_(v_to_var(readScalar(distributionModelDict_.lookup("expectation")),
            readScalar(distributionModelDict_.lookup("variance")))),
    a_(0.147)
{
    if (maxValue_ < minValue_)
    {
        FatalErrorInFunction
            << "Maximum value is smaller than the minimum value:"
            << "    maxValue = " << maxValue_ << ", minValue = " << minValue_
            << abort(FatalError);
    }
}


Foam::distributionModels::lognormal::lognormal(const lognormal& p)
:
    distributionModel(p),
    minValue_(p.minValue_),
    maxValue_(p.maxValue_),
    expectation_(p.expectation_),
    variance_(p.variance_),
    a_(p.a_)
{}


// * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::distributionModels::lognormal::m_to_mu
(
        const scalar m,
        const scalar v
)
{
    return log(m / sqrt(1.0 + (v / m*m)));
}


Foam::scalar Foam::distributionModels::lognormal::v_to_var
(
        const scalar m,
        const scalar v
)
{
    return log(1.0 + (v / m*m));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributionModels::lognormal::~lognormal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributionModels::lognormal::sample() const
{

    scalar a = erf((minValue_ - expectation_)/variance_);
    scalar b = erf((maxValue_ - expectation_)/variance_);

    scalar y = rndGen_.sample01<scalar>();
    scalar x = erfInv(y*(b - a) + a)*variance_ + expectation_;

    // Note: numerical approximation of the inverse function yields slight
    //       inaccuracies

    x = exp(x);

//    x = min(max(x, minValue_), maxValue_);

    return x;

//    return exp(-pow((lnr - log(mean_r)), 2) / 2.0 / pow(log(stdev),2))
//           / log(stdev)
//           / sqrt(2*constant::mathematical::pi);

}


Foam::scalar Foam::distributionModels::lognormal::minValue() const
{
    return minValue_;
}


Foam::scalar Foam::distributionModels::lognormal::maxValue() const
{
    return maxValue_;
}


Foam::scalar Foam::distributionModels::lognormal::meanValue() const
{
    return expectation_;
}


Foam::scalar Foam::distributionModels::lognormal::erfInv(const scalar y) const
{
    scalar k = 2.0/(constant::mathematical::pi*a_) + 0.5*log(1.0 - y*y);
    scalar h = log(1.0 - y*y)/a_;
    scalar x = sqrt(-k + sqrt(k*k - h));
    if (y < 0.0)
    {
        x *= -1.0;
    }
    return x;
}


// ************************************************************************* //
