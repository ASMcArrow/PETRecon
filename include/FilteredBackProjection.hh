#ifndef FILTEREDBACKPROJECTION_HH
#define FILTEREDBACKPROJECTION_HH

class PETParameters;

class FilteredBackProjection
{

public:
    FilteredBackProjection(unsigned short*** sinogram3D, PETParameters* parameters);
    ~FilteredBackProjection();

    float*** GetReconstructedImage() { return Image; }

private:
    void SortToSinogram2D(unsigned short ***sinogram3D);
    void ApplyFilter();
    void Reconstruct();

    float*** Sinogram2D;
    float*** Image;
    PETParameters* Parameters;
};

#endif // FILTEREDBACKPROJECTION_HH
