#include <vector>

using namespace std;

class PETParameters
{

public:
    PETParameters();
    ~PETParameters(){}

    const int GetMaxRingDiff() { return MaxRingDiff; }
    const int GetRingNum() { return RingNum; }
    const int GetTanBinNum() { return TanBinNum; }
    const int GetAngBinNum() { return AngBinNum; }
    const int GetTotAxialPlanes() { return TotAxialPlanes; }
    const int GetTotPlanes() { return TotPlanes; }
    const int GetTanCrystalNum() { return TanCrystalNum; }
    const int GetAxialCrystalNum() { return AxialCrystalNum; }
    const int GetTanSubmoduleNum() { return TanSubmoduleNum; }
    const int GetAxialSubmoduleNum() { return AxialSubmoduleNum; }
    const int GetTanModuleNum() { return TanModuleNum; }
    const int GetAxialModuleNum() { return AxialModuleNum; }
    const int GetTotCrystalNum() { return TotCrystalNum; }

private:
    int MaxRingDiff; // Maximum ring difference
    int RingNum; // Number of rings
    int TanBinNum; // Number of tangential bins in sinogram
    int AngBinNum; // Number of angular bins in sinogram
    int TanCrystalNum;
    int AxialCrystalNum;
    int TanSubmoduleNum;
    int AxialSubmoduleNum;
    int TanModuleNum;
    int AxialModuleNum;
    int TotAxialPlanes; // Total number of axial planes (Total number of cells in michelogram)
    int TotPlanes; // Total number of 2D planes (Total number of diagonals in michelogram)
    int TotCrystalNum; // Total number of crystals in ring
};

struct ID
{
    int BaseID;
    int RSectorID;
    int ModuleID;
    int SubmoduleID;
    int CrystalID;
    int LayerID;
};

struct LOR
{
   ID point1, point2;
};


