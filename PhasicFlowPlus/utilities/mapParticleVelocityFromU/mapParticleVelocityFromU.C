/*------------------------------- PhasicFlow ---------------------------------
  Copyright (C): Antonio Di Mauro
  email: antoniodimauro03@gmail.com
------------------------------------------------------------------------------
Licence:
  This file is part of PhasicFlow, a CFD-DEM stack built on phasicFlow and
  PhasicFlowPlus (www.cemf.ir). It is free software: you can redistribute it
  and/or modify it under the terms of the GNU General Public License v3 or
  any later version.

  It is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.
-----------------------------------------------------------------------------*/

#include <cstdint>

#include "fvCFD.H"
#include "interpolationCellPoint.H"

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <iterator>
#include <limits>
#include <regex>
#include <string>
#include <vector>

namespace
{

Foam::word detectTargetTime(Foam::Time& runTime)
{
    const Foam::instantList times = runTime.times();

    Foam::label chosenIndex = -1;
    Foam::scalar chosenValue = -std::numeric_limits<Foam::scalar>::max();

    for (Foam::label i = 0; i < times.size(); ++i)
    {
        const Foam::word tName = times[i].name();
        if (tName == "constant")
        {
            continue;
        }

        Foam::scalar tVal = 0;
        try
        {
            tVal = Foam::readScalar(Foam::IStringStream(tName)());
        }
        catch (...)
        {
            continue;
        }
        if (tVal > chosenValue)
        {
            chosenValue = tVal;
            chosenIndex = i;
        }
    }

    if (chosenIndex < 0)
    {
        FatalErrorInFunction
            << "Could not find any numeric time directory." << Foam::nl
            << Foam::exit(Foam::FatalError);
    }

    runTime.setTime(times[chosenIndex], chosenIndex);
    return runTime.timeName();
}

std::vector<Foam::point> readParticlePoints(const Foam::fileName& pStructurePath)
{
    std::ifstream input(pStructurePath.c_str(), std::ios::binary);
    if (!input)
    {
        FatalErrorInFunction
            << "Cannot open particle structure file: " << pStructurePath << Foam::nl
            << Foam::exit(Foam::FatalError);
    }

    const std::string text
    (
        (std::istreambuf_iterator<char>(input)),
        std::istreambuf_iterator<char>()
    );

    bool isBinary = false;
    {
        const std::size_t key = text.find("fileFormat");
        if (key != std::string::npos)
        {
            std::string value = text.substr(key, text.find('\n', key) - key);
            std::transform
            (
                value.begin(), value.end(), value.begin(),
                [](unsigned char c){ return std::tolower(c); }
            );
            isBinary = value.find("binary") != std::string::npos;
        }
    }

    const std::size_t tag = text.find("internalPoints");
    if (tag == std::string::npos)
    {
        FatalErrorInFunction
            << "No internalPoints section in: " << pStructurePath << Foam::nl
            << Foam::exit(Foam::FatalError);
    }

    std::size_t pos = tag + std::string("internalPoints").size();
    while (pos < text.size() && !std::isdigit(static_cast<unsigned char>(text[pos])))
    {
        ++pos;
    }

    if (pos >= text.size())
    {
        FatalErrorInFunction
            << "Could not parse number of points from: " << pStructurePath << Foam::nl
            << Foam::exit(Foam::FatalError);
    }

    std::size_t digitsEnd = pos;
    while (digitsEnd < text.size()
        && std::isdigit(static_cast<unsigned char>(text[digitsEnd])))
    {
        ++digitsEnd;
    }

    const Foam::label expectedPoints = Foam::readLabel
    (
        Foam::IStringStream(text.substr(pos, digitsEnd - pos))()
    );

    const std::size_t open = text.find('(', digitsEnd);
    if (open == std::string::npos)
    {
        FatalErrorInFunction
            << "No opening bracket after the point count in: " << pStructurePath
            << Foam::nl << Foam::exit(Foam::FatalError);
    }

    std::vector<Foam::point> points;
    points.reserve(expectedPoints);

    if (isBinary)
    {
        const std::size_t need =
            std::size_t(expectedPoints)*3*sizeof(double);

        if (text.size() < open + 1 + need)
        {
            FatalErrorInFunction
                << "Binary " << pStructurePath << " is short: " << expectedPoints
                << " points need " << need << " bytes after the bracket, only "
                << (text.size() - open - 1) << " are there." << Foam::nl
                << Foam::exit(Foam::FatalError);
        }

        for (Foam::label i = 0; i < expectedPoints; ++i)
        {
            double xyz[3];
            std::memcpy
            (
                xyz, text.data() + open + 1 + std::size_t(i)*3*sizeof(double),
                3*sizeof(double)
            );
            points.emplace_back(xyz[0], xyz[1], xyz[2]);
        }
    }
    else
    {
        const std::regex pointRegex(
            R"(\(([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\))"
        );

        const std::string body = text.substr(open + 1);

        for (std::sregex_iterator it(body.begin(), body.end(), pointRegex), end;
             it != end && Foam::label(points.size()) < expectedPoints; ++it)
        {
            const std::smatch& match = *it;

            points.emplace_back
            (
                Foam::readScalar(Foam::IStringStream(match[1].str())()),
                Foam::readScalar(Foam::IStringStream(match[2].str())()),
                Foam::readScalar(Foam::IStringStream(match[3].str())())
            );
        }
    }

    if (Foam::label(points.size()) != expectedPoints)
    {
        FatalErrorInFunction
            << "Parsed " << points.size() << " points from " << pStructurePath
            << " but expected " << expectedPoints
            << " (read as " << (isBinary ? "binary" : "ascii")
            << ", per that file's fileFormat header)." << Foam::nl
            << Foam::exit(Foam::FatalError);
    }

    return points;
}

void writePhasicFlowVelocity
(
    const Foam::fileName& velocityPath,
    const std::vector<Foam::vector>& velocities,
    const char* objectName = "velocity"
)
{
    std::ofstream output(velocityPath.c_str());
    if (!output)
    {
        FatalErrorInFunction
            << "Cannot write " << objectName << " file: " << velocityPath << Foam::nl
            << Foam::exit(Foam::FatalError);
    }

    output.precision(17);

    output
        << "/* -------------------------------*- C++ -*---------------------------------- *\\\n"
        << "|  phasicFlow File                                                            | \n"
        << "|  copyright: www.cemf.ir                                                     | \n"
        << "\\* ------------------------------------------------------------------------- */  \n"
        << " \n"
        << "objectType      pointField<realx3,Host>;\n"
        << "objectName      " << objectName << ";\n"
        << "fileFormat      ASCII;\n"
        << "\n"
        << "// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - // \n"
        << " \n"
        << "ineternalField  \n"
        << velocities.size() << "\n"
        << "(";

    for (const Foam::vector& velocity : velocities)
    {
        output
            << "(" << velocity.x() << " " << velocity.y() << " " << velocity.z() << ")\n";
    }

    output << ");\n";
}

}

int main(int argc, char* argv[])
{
    Foam::argList::addOption("time", "timeName", "Target time folder (default: latest numeric time)");
    Foam::argList::addBoolOption
    (
        "nospin",
        "Leave the particle spin (rVelocity) as it is. By default it is set to "
        "half the fluid vorticity at the particle, the torque-free spin of a "
        "sphere in a linear flow at small Reynolds number"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Foam::word targetTime;
    if (args.found("time"))
    {
        targetTime = args.get<Foam::word>("time");
        runTime.setTime(Foam::instant(targetTime), 0);
    }
    else
    {
        targetTime = detectTargetTime(runTime);
    }

    Foam::Info
        << "Mapping particle velocity from U at time " << targetTime << Foam::nl;

    const Foam::volVectorField U
    (
        Foam::IOobject
        (
            "U",
            targetTime,
            mesh,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::NO_WRITE
        ),
        mesh
    );

    const Foam::fileName pStructurePath = runTime.path()/targetTime/"pStructure";
    const std::vector<Foam::point> points = readParticlePoints(pStructurePath);

    std::vector<Foam::vector> velocities(points.size(), Foam::vector::zero);

    const bool setSpin = !args.found("nospin");
    std::vector<Foam::vector> spins(setSpin ? points.size() : 0, Foam::vector::zero);
    Foam::tmp<Foam::volVectorField> tCurlU;
    Foam::autoPtr<Foam::interpolationCellPoint<Foam::vector>> curlInterp;
    if (setSpin)
    {
        tCurlU = Foam::fvc::curl(U);
        curlInterp.reset(new Foam::interpolationCellPoint<Foam::vector>(tCurlU()));
    }

    Foam::label outsideCount = 0;
    for (Foam::label i = 0; i < static_cast<Foam::label>(points.size()); ++i)
    {
        const Foam::label cellId = mesh.findCell(points[i]);

        if (cellId >= 0)
        {
            velocities[i] = U[cellId];

            if (setSpin)
            {
                spins[i] = 0.5*curlInterp->interpolate(points[i], cellId);
            }
        }
        else
        {
            ++outsideCount;
        }
    }

    if (outsideCount > 0)
    {
        WarningInFunction
            << outsideCount << " particle(s) were outside the mesh and got zero velocity." << Foam::nl;
    }

    const Foam::fileName velocityPath = runTime.path()/targetTime/"velocity";
    writePhasicFlowVelocity(velocityPath, velocities);

    Foam::Info
        << "Wrote " << velocities.size() << " particle velocities to " << velocityPath << Foam::nl
        << Foam::endl;

    if (setSpin)
    {
        const Foam::fileName spinPath = runTime.path()/targetTime/"rVelocity";
        writePhasicFlowVelocity(spinPath, spins, "rVelocity");

        Foam::Info
            << "Wrote " << spins.size() << " particle spins (half the fluid vorticity) to "
            << spinPath << Foam::nl << Foam::endl;
    }

    return 0;
}
