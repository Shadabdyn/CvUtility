/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "CvCalc.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(CvCalc, 0);
    addToRunTimeSelectionTable(functionObject, CvCalc, dictionary);
}
}

// * * * * * * * * * * * * * * * * Protected members  * * * * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::CvCalc::createFileNames
(
    const dictionary& dict
) const
{
    DynamicList<word> names(1);

    // use type of the utility as specified in the dict as the top-level dir name
    const word objectType(dict.lookup("type"));

    // Name for file(MAIN_FILE=0)
    names.append(objectType);

    return names;
}


void Foam::functionObjects::CvCalc::writeFileHeader(const label i)
{
    // Find the correct file to write to from the enum defined in the header.
    switch (fileID(i))
    {
        case MAIN_FILE:
        {
            writeHeader(file(i), "Flow coefficient of the valve is measured at 2D upstream and 6D downstream");
            writeCommented(file(i), "Time [s] | Flow Coefficient CV [gpm/psi[-1/2]]");
            file() << endl;
            break; // exit the case structure
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled file index: " << i
                << abort(FatalError);
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::functionObjects::CvCalc::CvCalc
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    // NOTE: call the base class constructor
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),

    name_(name),
    active_(true),
    UName_("U"),
    pName_("p"),
    phiName_("phi"),
    inletName_("Inlet"),
    // NOTE: Read the face zone to integrate over. Get its name from the dict, find
    // it in the mesh, and get a reference to the list of its faces.
    twoDFaceZoneName_(dict.lookup("twoDFaceZoneName")),
    sixDFaceZoneName_(dict.lookup("sixDFaceZoneName")),

    twoDFaceZoneLabel_( mesh_.faceZones().findZoneID(twoDFaceZoneName_) ),
    sixDFaceZoneLabel_( mesh_.faceZones().findZoneID(sixDFaceZoneName_) ),

    twoDfaces_( mesh_.faceZones()[twoDFaceZoneLabel_] ),
    sixDfaces_( mesh_.faceZones()[sixDFaceZoneLabel_] )

    //Cv Convergence criteria parameters
    //interval_("interval"),
    //startTime_("startTime"),
    //CoIterations_("CoIterations"),
    //tolPercentage_("tolPercentage"),
    //tolDegree_("tolDegree")


{

    read(dict);

    // built-in logFiles method for creating file streams.
    resetNames(createFileNames(dict));

    if (active_)
    {
        Info << "Finished initialising " << type() << ": " << name_ << nl << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::CvCalc::~CvCalc()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::CvCalc::read(const dictionary& dict)
{
    if (active_)
    {
        UName_ = dict.lookupOrDefault<word>("UName", "U");
        pName_ = dict.lookupOrDefault<word>("pName", "p");
        phiName_ = dict.lookupOrDefault<word>("phiName", "phi");
        inletName_ = dict.lookupOrDefault<word>("inletName", "Inlet");

        interval_ = dict.lookupOrDefault<label>("interval", 50);
        startTime_ = dict.lookupOrDefault<label>("startTime", 300);
        CoIterations_ = dict.lookupOrDefault<label>("CoIterations", 50);
        tolPercentage_ = dict.lookupOrDefault<scalar>("tolPercentage", 0.5);
        tolDegree_ = dict.lookupOrDefault<scalar>("tolDegree", 5.0);

    }
    return true;
}

bool Foam::functionObjects::CvCalc::execute()
{
    if (active_)
    {

    }
    return true;
}

bool Foam::functionObjects::CvCalc::end()
{
    if (active_)
    {
        execute();
    }
    return true;
}

void Foam::functionObjects::CvCalc::timeSet()
{}

bool Foam::functionObjects::CvCalc::write()
{
    if (active_)
    {

        // Retrieve a reference to the velocity field
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

        surfaceScalarField phi = obr_.lookupObject<surfaceScalarField>(phiName_);


        // itnerpolate onto the faces
        surfaceVectorField Uface = fvc::interpolate(U);

        // If operating on mesh faces (faceZone, patch)

        //- Local list of face IDs
        labelList twoDfaceId_;
        labelList sixDfaceId_;

        //- Local list of patch ID per face
        labelList twoDfacePatchId_;
        labelList sixDfacePatchId_;

        //- List of +1/-1 representing face flip map
        //  (1 use as is, -1 negate)
        labelList twoDfaceSign_;
        labelList sixDfaceSign_;

        //- Global number of faces
        label nTwoDFaces_;
        label nSixDFaces_;

        //Actual 2D face creation
        DynamicList<label> faceIds(twoDfaces_.size());
        DynamicList<label> facePatchIds(twoDfaces_.size());
        DynamicList<label> faceSigns(twoDfaces_.size());

        forAll(twoDfaces_, i)
        {
            label facei = twoDfaces_[i];

            label faceId = -1;
            label facePatchId = -1;
            if (mesh_.isInternalFace(facei))
            {
                faceId = facei;
                facePatchId = -1;
            }
            else
            {
                facePatchId = mesh_.boundaryMesh().whichPatch(facei);
                const polyPatch& pp = mesh_.boundaryMesh()[facePatchId];
                if (isA<coupledPolyPatch>(pp))
                {
                    if (refCast<const coupledPolyPatch>(pp).owner())
                    {
                        faceId = pp.whichFace(facei);
                    }
                    else
                    {
                        faceId = -1;
                    }
                }
                else if (!isA<emptyPolyPatch>(pp))
                {
                    faceId = facei - pp.start();
                }
                else
                {
                    faceId = -1;
                    facePatchId = -1;
                }
            }

            if (faceId >= 0)
            {
                //if (twoDfaces_.flipMap()[i])
                //{
                    //faceSigns.append(-1);
                //}
                //else
                {
                    faceSigns.append(1);
                }
                faceIds.append(faceId);
                facePatchIds.append(facePatchId);
            }
        }

        twoDfaceId_.transfer(faceIds);
        twoDfacePatchId_.transfer(facePatchIds);
        twoDfaceSign_.transfer(faceSigns);
        nTwoDFaces_ = returnReduce(twoDfaceId_.size(), sumOp<label>());
        //Info << "The total number of faces in 2D face zone is" << nTwoDFaces_<<endl;

        debug = false;

        if (debug)
        {
            Pout<< "Original face zone size = " << twoDfaces_.size()
                << ", new size = " << nTwoDFaces_<< endl;
        }


        //Actual 6D face creation
        DynamicList<label> sixDFaceIds(sixDfaces_.size());
        DynamicList<label> sixDFacePatchIds(sixDfaces_.size());
        DynamicList<label> sixDFaceSigns(sixDfaces_.size());

        forAll(sixDfaces_, i)
        {
            label facei = sixDfaces_[i];

            label faceId = -1;
            label facePatchId = -1;
            if (mesh_.isInternalFace(facei))
            {
                faceId = facei;
                facePatchId = -1;
            }
            else
            {
                facePatchId = mesh_.boundaryMesh().whichPatch(facei);
                const polyPatch& pp = mesh_.boundaryMesh()[facePatchId];
                if (isA<coupledPolyPatch>(pp))
                {
                    if (refCast<const coupledPolyPatch>(pp).owner())
                    {
                        faceId = pp.whichFace(facei);
                    }
                    else
                    {
                        faceId = -1;
                    }
                }
                else if (!isA<emptyPolyPatch>(pp))
                {
                    faceId = facei - pp.start();
                }
                else
                {
                    faceId = -1;
                    facePatchId = -1;
                }
            }

            if (faceId >= 0)
            {
                //if (twoDfaces_.flipMap()[i])
                //{
                    //faceSigns.append(-1);
                //}
                //else
                {
                    faceSigns.append(1);
                }
                sixDFaceIds.append(faceId);
                sixDFacePatchIds.append(facePatchId);
            }
        }

        sixDfaceId_.transfer(sixDFaceIds);
        sixDfacePatchId_.transfer(sixDFacePatchIds);
        sixDfaceSign_.transfer(sixDFaceSigns);
        nSixDFaces_ = returnReduce(sixDfaceId_.size(), sumOp<label>());
        //Info << "The total number of faces in 6D face zone is" << nSixDFaces_<<endl;

        debug = false;

        if (debug)
        {
            Pout<< "Original face zone size = " << sixDfaces_.size()
                << ", new size = " << nSixDFaces_<< endl;
        }
        /*----------------------------------------------------------------------------*/


        //Flow rate calculation at 2D location
        scalar flowRate(0.0);

        const label patchid = mesh_.boundaryMesh().findPatchID(inletName_);

        //if ((patchid >= 0) && (phi.boundaryField()[patchid].size()))
        //label nInletFace = returnReduce(phi.boundaryField()[patchid].size(), sumOp<label>());

        if (patchid >= 0)
        {
            flowRate += gSum(phi.boundaryField()[patchid]);
        }

        //Info << "The patch id of Inlet is: " << patchid << " The size of patch is " << nInletFace<< endl;


        // reduce for parallel running
        // reduce(flowRate, sumOp<scalar>());

        flowRate =-1*flowRate;

        Info << "Total flow rate at Inlet is [m3/s]: " << flowRate << nl << endl;

        /*

        forAll(twoDfaceId_, faceI)
        {
            // Flow rate = dot product of velocity and surface area vector; in Latex terms,
            // Q = \mathbf{U} \cdot \mathbf{\hat{n}} A
        flowRate += Uface[twoDfaceId_[faceI]] & mesh_.Sf()[twoDfaceId_[faceI]];
        }

        // reduce for parallel running
        reduce(flowRate, sumOp<scalar>());

        Info << "Total flow rate " << flowRate << " through "
             << returnReduce(twoDfaceId_.size(), sumOp<label>()) << " faces" << nl << endl;*/

        /*---------------------------------------------------------------------------------*/
        //Average pressure calculation at 2D location
        // Retrieve a reference to the pressure field


        const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);

        // itnerpolate onto the faces
        surfaceScalarField pface = fvc::interpolate(p);

        scalar twoDPressure(0.0);
        scalar twoDArea(0.0);

        forAll(twoDfaceId_, faceI){

        twoDPressure += pface[twoDfaceId_[faceI]]* mesh_.magSf()[twoDfaceId_[faceI]];
        twoDArea += mesh_.magSf()[twoDfaceId_[faceI]];

        }


        // reduce for parallel running
        reduce(twoDPressure, sumOp<scalar>());


        reduce(twoDArea, sumOp<scalar>());


        twoDPressure =twoDPressure/twoDArea;


        Info << "2D average pressure on the "
             << returnReduce(twoDfaceId_.size(), sumOp<label>()) << " faces is [m2/s2]: " << twoDPressure << nl << endl;

        /*---------------------------------------------------------------------------------*/

        //Average pressure calculation at 6D location
        scalar sixDPressure(0.0);
        scalar sixDArea(0.0);

        forAll(sixDfaceId_, faceI){

        sixDPressure += pface[sixDfaceId_[faceI]]*mesh_.magSf()[sixDfaceId_[faceI]];
        sixDArea += mesh_.magSf()[sixDfaceId_[faceI]];

        }

        // reduce for parallel running
        reduce(sixDPressure, sumOp<scalar>());
        reduce(sixDArea, sumOp<scalar>());

        sixDPressure =sixDPressure/sixDArea;

        Info << "6D average pressure on the "
             << returnReduce(sixDfaceId_.size(), sumOp<label>()) << " faces is [m2/s2]: " << sixDPressure << nl << endl;


        /*---------------------------------------------------------------------------------*/

        scalar Cv(0.0);

        Cv = (flowRate*15850.32)/sqrt(abs(twoDPressure -sixDPressure)*998.98*0.000145038);

        Info << "The Flow Coefficent Cv is : " << Cv <<nl<<endl;

        if (Pstream::master())
        {

            logFiles::write();

            // Add the entry for this time step that has just been computed.
            file(MAIN_FILE) << obr_.time().value() << tab << Cv << endl;
        }

        //Cv Convergence Criteria Monitoring

        //CvList.append(Cv);
        //Info << " The current CvList is: " << CvList  <<nl<< endl;

        if(CvList.size() > startTime_ +interval_)
        {
            label n= CvList.size();
            label division = (n - (startTime_))/interval_;

            List<label> meanList;
            //meanList.append(startTime_ - interval_ );
            label i,j,k;


            for(i=0; i<= division; i++)
            {
                meanList.append(startTime_ + i*interval_);

            }

            //Info << "The mean List is :  " << meanList <<nl<< endl;

            //calculate mean of CV

            label m =meanList .size();

            //Info << " The size of mean list is " << m << nl <<endl;

            List<scalar> meanCvList;

            for (j=0; j<(m-1); j++ )
            {
                scalar sum = 0;
                //Info << "The current meantList is between " << meanList[j] <<" and " << meanList[j + 1] << nl << endl;

                for(k = meanList[j]; k < (meanList[j+1]); k++)
                {
                    //Info << "The current Cv value is " <<  CvList[k] << " and its index " << findIndex(CvList, CvList[k])<< nl << endl;
                    sum =sum + CvList[k];

                }

                meanCvList.append(sum/interval_);

            }

            //Info << " The mean Cv list is : " <<meanCvList << nl << endl;

            List<scalar> meanDiffCvList;

            for (j=1; j<(m-1); j++ )
            {


                for(k = meanList[j]; k < (meanList[j+1]); k++)
                {

                    meanDiffCvList.append((CvList[k] - meanCvList[j]/meanCvList[j]) * 100);
                    //Info << "The current Cv value is " <<  CvList[k] << " and its index " << findIndex(CvList, CvList[k])<< nl << endl;
                    //Info << "The current mean Cv is " << meanCvList[j] << nl << endl;

                }

            }

            //Info << " The complete meanDiffCvList is : " << meanDiffCvList << nl << endl;

            label percentageFlag=0;
            label degreeFlag=0;
            label sizeDiff = meanDiffCvList.size();

            for (i=0; i<sizeDiff; i++)
            {
                if(abs(meanDiffCvList[i]< tolPercentage_))
                {
                    percentageFlag = percentageFlag + 1;

                    if(percentageFlag> CoIterations_)
                    {
                        scalar slope = (CvList[startTime_ + i] - CvList[startTime_ + (i - CoIterations_)])/CoIterations_;
                        scalar degreeSlope = abs(atan(slope)*57.58);

                        if(degreeSlope < tolDegree_)
                        {
                            degreeFlag = degreeFlag + 1;

                            if(degreeFlag > CoIterations_)
                            {
                                Info <<" The solution is converged at iterations at " <<findIndex(CvList, CvList[startTime_ + i]) <<" iteration "
                                    <<" The final converged Cv is " << CvList[startTime_ + i]<<nl<< endl;
                                //fileName outputDir = mesh_.time().path()/"ConvergedCv";
                                //mkDir(outputDir);
                                autoPtr<OFstream> outputFilePtr;
                                outputFilePtr.reset(new OFstream("finalCv.dat"));
                                outputFilePtr() <<" The solution is converged at iterations at " <<findIndex(CvList, CvList[startTime_ + i]) <<" iteration "
                                                <<" The final converged Cv is " << CvList[startTime_ + i]<< endl;

                            }
                        }

                        else
                        {
                            degreeFlag =0;
                        }


                    }
                }
                else
                {
                    percentageFlag =0;
                }
            }





        }


    }
    return true;
}

// ************************************************************************* //
