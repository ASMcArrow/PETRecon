#ifndef FILTEREDBACKPROJECTION_HH
#define FILTEREDBACKPROJECTION_HH

class PETParameters;

class FilteredBackProjection
{

public:
    FilteredBackProjection(unsigned short*** sinogram3D, PETParameters* parameters);
    ~FilteredBackProjection();

    double*** GetReconstructedImage() { return Image; }

private:
    void SortToSinogram2D(unsigned short ***sinogram3D);
    void CalculateFilter();
    void FFT();
    void BackProject(int slice);


    double*** Sinogram2D;
    double*** Image;
    PETParameters* Parameters;
    int inputFFTSize;
    int outputFFTSize;
    double* Filter;

    double Offset;
    double* Costheta;
    double* Sintheta;
    double** XArray;
    double** YArray;

    double MaxBrightness;
};

#endif // FILTEREDBACKPROJECTION_HH
