#include "SinogramConverter.hh"
#include "PETParameters.hh"
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>

using namespace std;

SinogramConverter::SinogramConverter(char* inputFileName, PETParameters* parameters)
{
    std::ifstream inputfile(inputFileName);
    Parameters = parameters;

    // Initialization of 3D sinogram to zeros
    Sinogram3D = new unsigned short**[Parameters->GetTanBinNum()]; // TODO: check if it necessarily has to be unsigned short.
    for(int i = 0; i < Parameters->GetTanBinNum(); i++)
    {
        Sinogram3D[i] = new unsigned short*[Parameters->GetAngBinNum()];
        for(int j = 0; j < Parameters->GetAngBinNum(); j++)
        {
            Sinogram3D[i][j] = new unsigned short[Parameters->GetTotAxialPlanes()];
            for(int k = 0; k < Parameters->GetTotAxialPlanes(); k++)
            {
                Sinogram3D[i][j][k] = 0;
            }
        }
    }

    // Read LOR line by line from GATE ASCII output file for coincidences and add it to 3D array
    std::string line, dump;
    int base, rsector, module, submodule, crystal, layer;
    while (std::getline(inputfile, line))
    {
        LOR newLOR;

        std::istringstream iss(line);
        int n = 11;
        while(n-- > 0)
        {
            iss >> dump;
        }

        iss >> base;
        iss >> rsector;
        iss >> module;
        iss >> submodule;
        iss >> crystal;
        iss >> layer;

        ID point1;
        point1.BaseID = base;
        point1.RSectorID = rsector;
        point1.ModuleID = module;
        point1.SubmoduleID = submodule;
        point1.CrystalID = crystal;
        point1.LayerID = layer;

        n = 17;
        while(n-- > 0)
        {
            iss >> dump;
        }

        iss >> base;
        iss >> rsector;
        iss >> module;
        iss >> submodule;
        iss >> crystal;
        iss >> layer;

        ID point2;
        point2.BaseID = base;
        point2.RSectorID = rsector;
        point2.ModuleID = module;
        point2.SubmoduleID = submodule;
        point2.CrystalID = crystal;
        point2.LayerID = layer;

        newLOR.point1 = point1;
        newLOR.point2 = point2;

        AddLORtoSinogram3D(&newLOR);
    }
}

void SinogramConverter::AddLORtoSinogram3D(LOR* lor)
{
    int ring1 = (int)(lor->point1.CrystalID/Parameters->GetTanCrystalNum())
            + (int)(lor->point1.SubmoduleID/Parameters->GetTanSubmoduleNum())*Parameters->GetAxialCrystalNum()
            + (int)(lor->point1.ModuleID/Parameters->GetTanModuleNum())*Parameters->GetAxialSubmoduleNum()*Parameters->GetAxialCrystalNum();

    int ring2 = (int)(lor->point2.CrystalID/Parameters->GetTanCrystalNum())
            + (int)(lor->point2.SubmoduleID/Parameters->GetTanSubmoduleNum())*Parameters->GetAxialCrystalNum()
            + (int)(lor->point2.ModuleID/Parameters->GetTanModuleNum())*Parameters->GetAxialSubmoduleNum()*Parameters->GetAxialCrystalNum();

    if ((ring1 < 0)||(ring2 < 0)||
            (ring1 >= Parameters->GetRingNum())||
            (ring2 >= Parameters->GetRingNum()))
    {
        cout << "ERROR: Wrong ring number!" << endl;
        return;
    }

    if (abs(ring2-ring1) > Parameters->GetMaxRingDiff())
    {
        cout << "In this LOR the ring difference is bigger than the set maximum. Discard." << endl;
        return;
    }

    int crystal1 = lor->point1.RSectorID*Parameters->GetTanModuleNum()*Parameters->GetTanSubmoduleNum()*Parameters->GetTanCrystalNum()
            + (lor->point1.ModuleID%Parameters->GetTanModuleNum())*Parameters->GetTanSubmoduleNum()*Parameters->GetTanCrystalNum()
            + (lor->point1.SubmoduleID%Parameters->GetTanSubmoduleNum())*Parameters->GetTanCrystalNum()
            + (lor->point1.CrystalID%Parameters->GetTanCrystalNum());

    int crystal2 = lor->point2.RSectorID*Parameters->GetTanModuleNum()*Parameters->GetTanSubmoduleNum()*Parameters->GetTanCrystalNum()
            + (lor->point2.ModuleID%Parameters->GetTanModuleNum())*Parameters->GetTanSubmoduleNum()*Parameters->GetTanCrystalNum()
            + (lor->point2.SubmoduleID%Parameters->GetTanSubmoduleNum())*Parameters->GetTanCrystalNum()
            + (lor->point2.CrystalID%Parameters->GetTanCrystalNum());

    int phi = ((crystal1 + crystal2 + Parameters->GetTotCrystalNum()/2)%Parameters->GetTotCrystalNum())/2;

    int dis = 0;
    if (((crystal1 + crystal2) < (3*Parameters->GetTotCrystalNum()/2)) && ((crystal1 + crystal2) >= (Parameters->GetTotCrystalNum()/2)))
        dis    =  (abs(crystal1 - crystal2) -  Parameters->GetTotCrystalNum()/2)*0.5;
    else dis = (-abs(crystal1 - crystal2) +  Parameters->GetTotCrystalNum()/2)*0.5;

    int z = (ring1*Parameters->GetRingNum() + ring2);

    // Adding the coincidence to sinogram array
    Sinogram3D[Parameters->GetTanBinNum()/2+dis][phi][z]++;
}

