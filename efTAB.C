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
    ANY WARRANTY; without even the implied warranty of MERCHANefTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "efTAB.H"
#include "electromagneticConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::efTAB<CloudType>::efTAB
(
    const dictionary& dict,
    CloudType& owner
)
:
    BreakupModel<CloudType>(dict, owner, typeName, true),
    SMDCalcMethod_(this->coeffDict().lookup("SMDCalculationMethod"))
{
    // calculate the inverse function of the Rossin-Rammler Distribution
    const scalar xx0 = 12.0;
    const scalar rrd100 =
        1.0/(1.0 - exp(-xx0)*(1.0 + xx0 + sqr(xx0)/2.0 + pow3(xx0)/6.0));

    forAll(rrd_, n)
    {
        scalar xx = 0.12*(n + 1);
        rrd_[n] =
            (1.0 - exp(-xx)*(1.0 + xx + sqr(xx)/2.0 + pow3(xx)/6.0))*rrd100;
    }

    if (SMDCalcMethod_ == "method1")
    {
        SMDMethod_ = method1;
    }
    else if (SMDCalcMethod_ == "method2")
    {
        SMDMethod_ = method2;
    }
    else
    {
        SMDMethod_ = method2;
        WarningInFunction
            << "Unknown SMDCalculationMethod. Valid options are "
            << "(method1 | method2). Using method2" << endl;
    }
}


template<class CloudType>
Foam::efTAB<CloudType>::efTAB(const efTAB<CloudType>& bum)
:
    BreakupModel<CloudType>(bum),
    SMDCalcMethod_(bum.SMDCalcMethod_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::efTAB<CloudType>::~efTAB()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::efTAB<CloudType>::update
(
    const scalar dt,
    const vector& g,
    scalar& d,
    scalar& tc,
    scalar& ms,
    scalar& nParticle,
    scalar& KHindex,
    scalar& y,
    scalar& yDot,
    const scalar d0,
    const scalar rho,
    scalar& Qdens,
    const scalar mu,
    const scalar sigma,
    scalar& sigma_s,
    const scalar Qdens0,
    scalar& netChargeChild,
    const vector& U,
    const scalar rhoc,
    const scalar muc,
    const vector& Urel,
    const scalar Urmag,
    const scalar tMom,
    scalar& dChild,
    scalar& massChild
)
{
    Random& rndGen = this->owner().rndGen();

    scalar r = 0.5*d;
    scalar r2 = r*r;
    scalar r3 = r*r2;

    scalar semiMass = nParticle*pow3(d);

    //scalar Qnet0 = Qdens0*pow(d0,3.0)*pi/6.0;  //need access to Qdens0 !
    //Qdens= Qnet0*6/(pi*pow(d,3.0)); //update of Qdens
    Info << "Qdens_efTAB_test1:" << Qdens << endl;

    scalar Qdens_old = Qdens;
    scalar dParent_old = d;
    scalar nParticle_old = nParticle;
    Info << "nParticle_old_efTAB:" << nParticle_old << endl;
    Info << "dParent_old_efTAB:" << dParent_old << endl;

    sigma_s = sigma - pow((Qdens_old*pi/6.0*pow(d0, 3.0)),2)/(8.0*Foam::constant::electromagnetic::epsilon0.value()*pi*pi*pow(dParent_old,3.0));

    Info << "Sigma_s_efTAB:" << sigma_s << endl;

    scalar netChargeParent = Qdens_old*pow(d0,3.0)*pi/6.0;
    Info << "netChargeparent_efTAB_test:" << netChargeParent << endl;
    Info << "Qdens0_efTAB_test:" << Qdens0 << endl;

    // inverse of characteristic viscous damping time
    scalar rtd = 0.5*this->TABCmu_*mu/(rho*r2);

    // oscillation frequency (squared)
    scalar omega2 = this->TABComega_*sigma_s/(rho*r3) - rtd*rtd;

    Info << "omega2_efTAB:" << omega2 << endl;

    //scotty we may have a Coulomb break-up
    if (sigma_s <= 0)
    {
        scalar We = rhoc*sqr(Urmag)*r/sigma_s;
        Info << "We_s_efTAB:" << We << endl;
        //update droplet size
        scalar rs = r;
        scalar rNew = rs/10;
        d = 2*rNew;
    }

    if (omega2 > 0)
    {
        scalar omega = sqrt(omega2);
        scalar We = rhoc*sqr(Urmag)*r/sigma_s;
        Info << "We_efTAB:" << We << endl;
        scalar Wetmp = We/this->TABtwoWeCrit_;

        scalar y1 = y - Wetmp;
        scalar y2 = yDot/omega;

        scalar a = sqrt(y1*y1 + y2*y2);

        // scotty we may have break-up
        if (a+Wetmp > 1.0)
        {

            Info << "BREAKUP !!!" << endl;
            scalar phic = y1/a;

            // constrain phic within -1 to 1
            phic = max(min(phic, 1), -1);

            scalar phit = acos(phic);
            scalar phi = phit;
            scalar quad = -y2/a;
            if (quad < 0)
            {
                phi = constant::mathematical::twoPi - phit;
            }

            scalar tb = 0;

            if (mag(y) < 1.0)
            {
                scalar coste = 1.0;
                if ((Wetmp - a < -1) && (yDot < 0))
                {
                    coste = -1.0;
                }

                scalar theta = acos((coste-Wetmp)/a);

                if (theta < phi)
                {
                    if (constant::mathematical::twoPi - theta >= phi)
                    {
                        theta = -theta;
                    }
                    theta += constant::mathematical::twoPi;
                }
                tb = (theta-phi)/omega;

		Info << "tb and dt " << tb << "; " << dt << endl;

                // breakup occurs
                if (dt > tb)
                {
                    y = 1.0;
                    yDot = -a*omega*sin(omega*tb + phi);

		    Info << "It's really a BREAKUP " << endl;
                }

            }

            // update droplet size
            if (dt > tb)
            {
                scalar rs =
                    r/(1.0 + (4.0/3.0)*sqr(y) + rho*r3/(8*sigma_s)*sqr(yDot));

                label n = 0;
                scalar rNew = 0.0;
                switch (SMDMethod_)
                {
                    case method1:
                    {
                        #include "efTABSMDCalcMethod1.H"
                        break;
                    }
                    case method2:
                    {
                        #include "efTABSMDCalcMethod2.H"
                        break;
                    }
                }

                if (rNew < r)
                {
                    d = 2*rNew;
                    y = 0;
                    yDot = 0;
                }
            } 
        }
    }
    else
    {
        // reset droplet distortion parameters
        y = 0;
        yDot = 0;
    }

    // update the nParticle count to conserve mass
    nParticle = semiMass/pow3(d);

    Info << "nParticle_efTAB:" << nParticle << endl;
    Info << "d_efTAB:" << d << endl;

    // update droplet charge density to conserve charge
    netChargeParent = Qdens_old*pow(d0,3.0)*pi/6.0;
    Info << "netChargeparent_efTAB:" << netChargeParent << endl;
    netChargeChild = (nParticle_old/nParticle)*netChargeParent;
    //Qnet0= netChargeChild;
    //Info << "Qnet0_efTAB:" << Qnet0 << endl;
    Info << "netChargeChild_efTAB:" << netChargeChild << endl;
    //QdensChild = netChargeChild/(pow(d,3.0)*pi/6.0);
    Qdens = netChargeChild/(pow(d0,3.0)*pi/6.0);
    Info << "Qdens_efTAB:" << Qdens << endl;



    // Do not add child parcel
    return false;
}


// ************************************************************************* //
