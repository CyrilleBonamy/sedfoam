/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Application
    SedFoam

Description
    Solver for a system of 2 phases with one phase dispersed,
    e.g.  solid particles in a fluid.

    Reference:
    \verbatim
        Version 2.0:
        Chauchat J., Cheng Z., Nagel T., Bonamy C., & Hsu T-J. (2017).
        SedFoam-2.0: a 3D two-phase flow numerical model for sediment transport
        Geosci. Model Dev. Discuss.
        http://dx.doi.org/10.5194/gmd-2017-101
    \endverbatim


\*---------------------------------------------------------------------------*/
/**
 * \file sedFoam.C
 * \brief 2 phases Solver
 * \author Julien Chauchat, Cyrille Bonamy, Antoine Mathieu, Tim Nagel,
           Zhen Cheng, Tian-Jian Hsu and Eduard Puig Montella.
 * \version 3.1
 * \date September 16, 2019
 *
 * Solver for a system of 2 phases with one phase dispersed
 *
 */
/*
* Changelog [Higuera]
* December 10, 2018 - Adapted to compile automatically with OpenFOAM 6.
*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "PhaseIncompressibleTurbulenceModel.H"

#include "symmetryFvPatchFields.H"
#include "fixedFluxPressureFvPatchScalarField.H"

#include "dragModel.H"
#include "phaseModel.H"
#include "ppModel.H"

#include "kineticTheoryModel.H"
#include "granularRheologyModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "fixedValueFvsPatchFields.H"

#include "cellCellStencilObject.H"
#include "zeroGradientFvPatchFields.H"
#include "localMin.H"
#include "interpolationCellPoint.H"
#include "transform.H"
#include "fvMeshSubset.H"
#include "oversetAdjustPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow"
        " on a moving mesh."
    );

//    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"

    pimpleControl pimple(mesh);

    #include "readGravity.H"
    #include "createFields.H"
    #include "createRASTurbulence.H"
    #include "createUf.H"
    #include "createMRF.H"
    #include "createFvOptions.H"
    #include "createControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    if (SUSlocal)
    {
        Info<< "\nLocal Schmidt number activated" << endl;
    }
    else
    {
       if (max(SUS).value() == 0)
       {
           Info<< "Turbulence suspension term is neglected" << endl;
       }
       else if (max(SUS).value() > 0)
       {
           Info<< "Turbulence suspension term is included" << endl;
       }
       else
       {
           Info<< "Turbulence suspension coefficient SUS can't be negative"
               << endl;
       }
    }

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTwoPhaseEulerFoamControls.H"
        #include "CourantNos.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "gravityRamp.H"

        bool changed = mesh.update();

        if (changed)
        {            
            Info << "MESH CHANGED" << endl;
            #include "setCellMask.H"
            #include "setInterpolatedCells.H"

            surfaceScalarField faceMaskOld
            (
                localMin<scalar>(mesh).interpolate(cellMask.oldTime())
            );

            // Zero Uf on old faceMask (H-I)
            Ufa *= faceMaskOld;
            Ufb *= faceMaskOld;
            // Update Uf and phi on new C-I faces

            const surfaceVectorField Uinta(fvc::interpolate(Ua));
            const surfaceVectorField Uintb(fvc::interpolate(Ub));
            // Update Uf and phi on new C-I faces
            Ufa += (1-faceMaskOld)*Uinta;
            Ufb += (1-faceMaskOld)*Uintb;

            // Update Uf boundary
            forAll(Ufa.boundaryField(), patchI)
            {
                Ufa.boundaryFieldRef()[patchI] =
                    Uinta.boundaryField()[patchI];
            }
            forAll(Ufb.boundaryField(), patchI)
            {
                Ufb.boundaryFieldRef()[patchI] =
                    Uintb.boundaryField()[patchI];
            }



            phia = mesh.Sf() & Ufa;
            phib = mesh.Sf() & Ufb;



            // Zero phi on current H-I
            surfaceScalarField faceMask
            (
                localMin<scalar>(mesh).interpolate(cellMask)
            );
            phia *= faceMask;
            phib *= faceMask;
            Ua *= cellMask;
            Ub *= cellMask;
        // Make the flux relative to the mesh motion
        fvc::makeRelative(phia, Ua);
        fvc::makeRelative(phib, Ub);
        surfaceScalarField alphaf = fvc::interpolate(alpha);
        surfaceScalarField betaf = scalar(1.0) - alphaf;
        phi = alphaf*phia + betaf*phib;
                    phi *= faceMask;
                    U   *= cellMask;
                    alpha   *= cellMask;
                    // Make the flux relative to the mesh motion
                  //  fvc::makeRelative(phi, U);
        }


        if (mesh.changing() && correctPhi)
        {
            // Calculate absolute flux from the mapped surface velocity
            #include "correctPhi.H"
        }


        if (mesh.changing() && checkMeshCourantNo)
        {
            #include "meshCourantNo.H"
        }


//      Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "alphaEqn.H"
            #include "liftDragCoeffs.H"

//          Compute the Kinetic Theory parameters: nuEffa and lambdaUa from the
//          solution of the Granular temperature equation
            #include "callKineticTheory.H"

//          Compute the contact pressure pff and the frictional viscosity nuFra
//          from a Coulomb model if using the kinetic theory
//          and from the mu(I) rheology if using the granular rheology
            #include "callFrictionStress.H"
Info << "pass callfric" << endl;
//          Create the momentum balance equations for both phases a and b
            #include "UEqns.H"

Info << "pass UEqn" << endl;
            #include "pEqn.H"

            #include "DDtU.H"

            if (pimple.turbCorr())
            {
                #include "updateTwoPhaseRASTurbulence.H"
                turbulenceb->correct();
                if (debugInfo)
                {
                    Info << " max(nutb) = "
                         << max(turbulenceb->nut()).value() << endl;
                }
            }
        }
        if (debugInfo)
        {
            Info<< "min(Ua) = " << gMin(Ua)
                << "max(Ua) = " << gMax(Ua) << endl;
            Info<< "min(Ub) = " << gMin(Ub)
                << "max(Ub) = " << gMax(Ub) << nl << endl;
        }
        #include "OutputGradPOSC.H"
        #include "writeOutput.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}
// ************************************************************************* //
