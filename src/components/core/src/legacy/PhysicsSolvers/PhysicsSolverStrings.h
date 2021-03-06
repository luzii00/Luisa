//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  GEOS Computational Framework - Core Package, Version 3.0.0
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)
//  Stuart Walsh(walsh24@llnl.gov)
//  Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//  Chandrasekhar Annavarapu Srinivas
//  Eric Herbold
//  Michael Homel
//
//
//  All rights reserved.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
// SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
// TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S.
// Department of Energy (DOE). This work was produced at Lawrence Livermore
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National
// Security, LLC nor any of their employees, makes any warranty, express or
//     implied, or assumes any liability or responsibility for the accuracy,
// completeness, or usefulness of any information, apparatus, product, or
//     process disclosed, or represents that its use would not infringe
// privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or
// services by trade name, trademark, manufacturer or otherwise does not
//     necessarily constitute or imply its endorsement, recommendation, or
// favoring by the United States Government or Lawrence Livermore National
// Security,
//     LLC. The views and opinions of authors expressed herein do not
// necessarily state or reflect those of the United States Government or
// Lawrence
//     Livermore National Security, LLC, and shall not be used for advertising
// or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The
// BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef PHYSICSSOLVERSTRINGS_H_
#define PHYSICSSOLVERSTRINGS_H_

#include <string>
/*
 * PhysicsSolverStrings.h
 *
 *  Created on: Feb 18, 2011
 *      Author: walsh24
 */

/// Field and Model Parameter Strings
namespace PS_STR
{

const static std::string HeadStr = "Head";
const static std::string PressureStr = "Pressure";
const static std::string DarcyFluxStr = "DarcyFlux";
const static std::string FluidVelocityStr = "FluidVelocity";

const static std::string PermeabilityStr = "Permeability";
const static std::string DiffusivityStr = "Diffusivity";

const static std::string ConcentrationStr = "Concentration";
const static std::string ConcentrationFluxStr = "ConcentrationFlux";
const static std::string ReactionRateStr = "ReactionRate";

const static std::string PorosityStr = "Porosity";
const static std::string ApertureStr = "Aperture";
const static std::string MinimumApertureStr = "MinimumAperture";
const static std::string MaximumApertureStr = "MaximumAperture";

const static std::string BulkModulusStr = "BulkModulus";

const static std::string VolumetricFluxStr =  "VolumetricFlux";

const static std::string ElementCenterStr =  "ElementCenter";

const static std::string FaceAreaStr =  "FaceArea";
const static std::string FaceCenterStr =  "FaceCenter";
const static std::string FaceNormalStr =  "FaceNormal";

const static std::string EdgeCenterStr =  "EdgeCenter";
const static std::string EdgeLengthStr =  "EdgeLength";

const static std::string VerboseStr = "Verbose";

const static std::string ViscosityStr = "Viscosity";

const static std::string ProppantVolumeFractionStr = "ProppantVolumeFraction";
const static std::string ProppantPackVolumeFractionStr = "ProppantPackVolumeFraction";
}

#endif /*PHYSICSSOLVERSTRINGS_H_*/
