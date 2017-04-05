#include "PETParameters.hh"

PETParameters::PETParameters()
{
    RingNum = 32; // Number of rings
    MaxRingDiff = 31; // Maximum ring difference
    TanBinNum = 288; // Number of tangential bins in sinogram
    AngBinNum = 288; // Number of angular bins in sinogram
    TanCrystalNum = 8;
    AxialCrystalNum = 8;
    TanSubmoduleNum = 1;
    AxialSubmoduleNum = 1;
    TanModuleNum = 1;
    AxialModuleNum = 4;
    TotAxialPlanes = RingNum*RingNum;
    TotPlanes = RingNum*2-1;
    TotCrystalNum = 576;
}
