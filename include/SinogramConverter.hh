#ifndef SINOGRAMCONVERTER_H
#define SINOGRAMCONVERTER_H

class PETParameters;
struct LOR;

class SinogramConverter
{
public:
    SinogramConverter(char* inputFileName, PETParameters* parameters);
    unsigned short*** GetSinogramArray() { return Sinogram3D; }

private:
    void AddLORtoSinogram3D(LOR* lor);

    unsigned short*** Sinogram3D;
    PETParameters* Parameters;
};

#endif // SINOGRAMCONVERTER_H
