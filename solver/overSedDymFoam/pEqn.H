/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
{
    word schemeA("div(alpha,nu)");

    surfaceScalarField faceMask(localMin<scalar>(mesh).interpolate(cellMask));

    surfaceScalarField alphaf(fvc::interpolate(alpha));
    surfaceScalarField betaf(scalar(1.0) - alphaf);

    // Compute inverse matrix for momentum balance of phase a and b
    //     at cell centers
    volScalarField rUaA(1.0/UaEqn.A());
    volScalarField rUbA(1.0/UbEqn.A());
    mesh.interpolate(rUaA);
    mesh.interpolate(rUbA);

    // Interpolate inverse matrix for momentum balance of phase a and b
    //     at cell faces
    rUaAf = fvc::interpolate(rUaA);
    surfaceScalarField rUbAf(fvc::interpolate(rUbA));

    // limiter for the explicit stress terms to avoid numerical instability and
    // too strong limitation on the time step
    volScalarField nuFraExp(nuFra);
    nuFraExp.min(nuMaxExp);
    nuFraExp.correctBoundaryConditions();

    // interpolate 1/(alpha+alphaSmall) on the faces
    surfaceScalarField oneOverAlpha
    (
        fvc::interpolate
        (
            scalar(1.0)/(alpha+alphaSmall),
            "interpolate(oneoveralpha)"
        )
    );

    // local Schmidt number
    if (SUSlocal)
    {
        volScalarField ustar
        (
            sqrt(mag(beta*nuEffb*(fvc::grad(Ub)+gradUbT)))
        );
        SUStilde = scalar(1.0) + scalar(2.0)*pow(Ufall/(ustar+Usmall), 2);
        SUStilde.min(scalar(3.0));
        SUS = SUStilde;

        Info<<"ustarmax"<<max(ustar).value()<<endl;
    }

    // PIMPLE algorithm loop
    while (pimple.correct())
    {
        // Compute predicted velocity of phase a and b at cell centers
        volVectorField HabyA(Ua);
        HabyA = rUaA*UaEqn.H();
        volVectorField HbbyA(Ub);
        HbbyA = rUbA*UbEqn.H();

//if (massFluxInterpolation)
//{
//    #include "interpolatedFaces.H"
//}

        // Momentum fluxes computed on the faces for the particle phase phiDraga
        // and the fluid phase phiDragb
        surfaceScalarField phiDraga(
        faceMask*(
        // turbulent suspension term
           -(fvc::interpolate
            (
                SUS*(K/rhoa)*turbulenceb->nut()*(beta/(alpha + alphaSmall))
            )
           *fvc::snGrad(alpha))*rUaAf*mesh.magSf()
            // grad(pa)
          - (oneOverAlpha/rhoa*fvc::snGrad(pa))*rUaAf*mesh.magSf()
            // grad(Pf)
             - (oneOverAlpha/rhoa*fvc::snGrad(pff))*rUaAf*mesh.magSf()
            // diagonal part of stress
          + (oneOverAlpha/rhoa*fvc::snGrad(lambdaUa*tr(gradUaT)))
           *rUaAf*mesh.magSf()
            // diagonal part of stress 2
          - (
                scalar(2.0/3.0)*oneOverAlpha
               *fvc::snGrad((alpha+alphaSmall)*(sqr(Ct)*turbulenceb->k()
              + nuEffa*tr(gradUaT)))
            )
           *rUaAf*mesh.magSf()
          - (scalar(2.0/3.0)*oneOverAlpha*fvc::snGrad(nuFraExp*tr(gradUaT)))
           *rUaAf*mesh.magSf()
           // stress part 3
          + rUaAf*(oneOverAlpha
           *fvc::interpolate
            (
                fvc::div((alpha+alphaSmall)*nuEffa*gradUaT, schemeA)
            )
          & mesh.Sf())
          + rUaAf*(oneOverAlpha
           *fvc::interpolate(fvc::div(nuFraExp*gradUaT, schemeA)) & mesh.Sf())
            // gravity
          + rUaAf
           *(
               (1.0-rhob/rhoa)*(gravity & mesh.Sf())
            )
          // mean pressure gradient
          + rUaAf*((gradPOSC/(tilt*rhob +(1.0-tilt)*rhoa)) & mesh.Sf())
        ));

        surfaceScalarField phiDragb(
        faceMask*(
        // turbulent suspension term
            fvc::interpolate(SUS*(K/rhob)*turbulenceb->nut())*fvc::snGrad(alpha)
            *rUbAf*mesh.magSf()
            // diagnoal part of stress
          - (fvc::interpolate(scalar(2.0)/((beta+alphaSmall)*scalar(3.0)))
           *fvc::snGrad(((nuEffb*tr(gradUbT)))*beta))*rUbAf*mesh.magSf()
           // stress part 3
          + rUbAf*(fvc::interpolate(scalar(1.0)/(beta+alphaSmall))
           *fvc::interpolate
            (
                fvc::div(beta*nuEffb*gradUbT, schemeA)) & mesh.Sf()
            )
          // driving force
          + rUbAf*((gradPOSC/rhob) & mesh.Sf())
          // Convective acceleration for 1D non-uniform simulation
          - rUbAf*fvc::interpolate(neg(alpha-0.3))*((gradConvAcc/rhob)
            & mesh.Sf())
        ));

        //    fvOptions.makeRelative(phiHabyA);
        //    fvOptions.makeRelative(phiHbbyA);


        // Fix for gravity on outlet boundary.
        forAll(p_rbgh.boundaryField(), patchi)
    {
      if (isA<zeroGradientFvPatchScalarField>(p_rbgh.boundaryField()[patchi]))
      {
          phiDraga.boundaryFieldRef()[patchi] = 0.0;
          phiDragb.boundaryFieldRef()[patchi] = 0.0;
      }
        if (isA<symmetryFvPatchScalarField>(p_rbgh.boundaryField()[patchi]))
        {
            phiDraga.boundaryFieldRef()[patchi] = 0.0;
            phiDragb.boundaryFieldRef()[patchi] = 0.0;
        }
        if (p_rbgh.boundaryField().types()[patchi] == "symmetryPlane")
        {
            phiDraga.boundaryFieldRef()[patchi] = 0.0;
            phiDragb.boundaryFieldRef()[patchi] = 0.0;
        }
        if (p_rbgh.boundaryField().types()[patchi] == "groovyBC")
        {
            phiDraga.boundaryFieldRef()[patchi] = 0.0;
            phiDragb.boundaryFieldRef()[patchi] = 0.0;
        }
        if (p_rbgh.boundaryField().types()[patchi] == "outletInlet")
        {
            phiDraga.boundaryFieldRef()[patchi] = 0.0;
            phiDragb.boundaryFieldRef()[patchi] = 0.0;
        }
        if (p_rbgh.boundaryField().types()[patchi] == "fixedFluxPressure")
        {
            phiDraga.boundaryFieldRef()[patchi] = 0.0;
            phiDragb.boundaryFieldRef()[patchi] = 0.0;
        }
    }
    surfaceScalarField faceMaskOld
    (
        localMin<scalar>(mesh).interpolate(cellMask.oldTime())
    );

        surfaceScalarField phiUa
        (
            (fvc::interpolate(HabyA) & mesh.Sf())
           + faceMaskOld*rUaAf*fvc::ddtCorr(Ua, Ufa)
        );

        surfaceScalarField phiUb
        (
           (fvc::interpolate(HbbyA) & mesh.Sf())
          + faceMaskOld*rUbAf*fvc::ddtCorr(Ub, Ufb)
        );

        surfaceScalarField phiaHbyA(phiUa + phiDraga);
        surfaceScalarField phibHbyA(phiUb + phiDragb);

        // Face-drag coefficients
        surfaceScalarField rAUKda(fvc::interpolate(beta/rhoa*K*rUaA));
        surfaceScalarField rAUKdb(fvc::interpolate(alpha/rhob*K*rUbA));
        // Compute predicted velocity fluxes at cell faces including
        //       including semi-implicit drag contribution
        MRF.makeRelative(phiaHbyA);
        MRF.makeRelative(phibHbyA);
        surfaceScalarField phiHbyA
        (
            alphaf*(phiaHbyA + faceMask*rAUKda*phib)
          + betaf*(phibHbyA + faceMask*rAUKdb*phia)
        );

    //MRF.makeRelative(phiHbyA);
    //adjustPhi(phiHbyA, U, p_rbgh);
    //surfaceScalarField phiHbyA = alphaf*phiaHbyA + betaf*phibHbyA;

    //    adjustPhi(phiHbyA, U, p_rbgh);

        // Poisson equation for the pressure
        surfaceScalarField Dp
        (
            "(rho*(1|A(U)))",
            alphaf*rUaAf/rhoa + betaf*rUbAf/rhob
        );

        // Update the fixedFluxPressure BCs to ensure flux consistency
        setSnGrad<fixedFluxPressureFvPatchScalarField>
        (
            p_rbgh.boundaryFieldRef(),
            (
                phiHbyA.boundaryField()
              - alphaf.boundaryField()*phia.boundaryField()
               //*MRF.relative(mesh.Sf().boundaryField() & Ua.boundaryField())
              + betaf.boundaryField()*phib.boundaryField()
               //*MRF.relative(mesh.Sf().boundaryField() & Ub.boundaryField())
             )
            /(mesh.magSf().boundaryField()*Dp.boundaryField())
        );

//    constrainPressure(p_rbgh, U, phiHbyA, Dp, MRF);


        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(Dp, p_rbgh) == fvc::div(phiHbyA)
            );

            pEqn.setReference(pRefCell, getRefCellValue(p_rbgh, pRefCell));

#if OFVERSION >= 700
            pEqn.solve();
#else
            pEqn.solve(mesh.solver(p_rbgh.select(pimple.finalInnerIter())));
#endif

            if (pimple.finalNonOrthogonalIter())
            {
                // Compute pressure gradient at cell faces
                surfaceScalarField SfGradp(pEqn.flux()/Dp);
                // Compute corrected velocity fluxes at cell faces
                //       due to pressure gradient
                surfaceScalarField phias(phiaHbyA - rUaAf*SfGradp/rhoa);
                surfaceScalarField phibs(phibHbyA - rUbAf*SfGradp/rhob);
                // Compute mixture velocity flux at cell faces
                phi = phiHbyA - pEqn.flux();
                // Compute corrected relative velocity fluxes due to drag force
                //       at cell faces
                surfaceScalarField phir
                (
                    ((phias + rAUKda*phibs) - (phibs + rAUKdb*phias))
                   /(1 - rAUKda*rAUKdb)
                );
                // Compute corrected phase velocity fluxes at cell faces
                //       using updated semi-implicit drag
                phia = phi + betaf*(phir);
                phib = phi - alphaf*(phir);

                #include "continuityErrs.H"
                p_rbgh.relax();

                // Velocity correction step at cell centers
                volVectorField Ubs
                (
                    HbbyA + rUbA*fvc::reconstruct(phiDragb/rUbAf-SfGradp/rhob)
                );
                volVectorField Uas
                (
                    HabyA + rUaA*fvc::reconstruct(phiDraga/rUaAf-SfGradp/rhoa)
                );
                //Uas.correctBoundaryConditions();
                //Ubs.correctBoundaryConditions();
                volScalarField Da(beta/rhoa*K*rUaA);
                volScalarField Db(alpha/rhob*K*rUbA);
                U = cellMask*(alpha*(Uas + Da*Ub) + beta*(Ubs + Db*Ua));
                volVectorField Ur(cellMask*((1 - Db)*Uas - (1 - Da)*Ubs)
                    /(1 - Da*Db));

                Ua = U + beta*Ur;
                Ua.correctBoundaryConditions();
                //fvOptions.correct(Ua);

                Ub = U - alpha*Ur;
                Ub.correctBoundaryConditions();
                //fvOptions.correct(Ub);

                p = p_rbgh + rhob*gh;
                if (p_rbgh.needReference())
                {
                    p += dimensionedScalar
                    (
                        "p",
                        p.dimensions(),
                        pRefValue - getRefCellValue(p, pRefCell)
                    );
                    p_rbgh = p - rhob*gh;
                    p_rbgh.correctBoundaryConditions();
                }
                #include "UaUbLimiter.H"
                Ua.correctBoundaryConditions();
                // Update mixture velocity
                U = alpha*Ua + beta*Ub;
                U.correctBoundaryConditions();
{
    Ufa = fvc::interpolate(Ua);
    Ufb = fvc::interpolate(Ub);
    surfaceVectorField n(mesh.Sf()/mesh.magSf());
    Ufa += n*(phia/mesh.magSf() - (n & Ufa));
    Ufb += n*(phib/mesh.magSf() - (n & Ufb));
}
    //fvc::makeRelative(phi, U);
    fvc::makeRelative(phia, Ua);
    fvc::makeRelative(phib, Ub);

    phia *= faceMask;
    phib *= faceMask;
    phi = alphaf*phia + betaf*phib;

            }
        }
    }
}
