#include "PETParameters.hh"
#include "SinogramConverter.hh"
#include "FilteredBackProjection.hh"

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkImageMapper.h>
#include <vtkActor2D.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkImageThreshold.h>
#include <vtkCamera.h>

#include <stdio.h>

int main(int argc,char *argv[])
{
    char* inputFileName = argv[1];

    PETParameters* parameters = new PETParameters; // Pay attention to parameters inside this file. They are hardware specific.
    SinogramConverter* converter = new SinogramConverter(inputFileName, parameters);
    FilteredBackProjection* projector = new FilteredBackProjection(converter->GetSinogramArray(), parameters);

    double*** image = projector->GetReconstructedImage();

    // VTK visualization
//    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
//    imageData->SetExtent(0, parameters->GetTanBinNum()-1, 0, parameters->GetTanBinNum()-1, 0, 0);
//    imageData->SetSpacing(1.0, 1.0, 1.0);
//    imageData->SetOrigin(0.0, 0.0, 0.0);
//    imageData->AllocateScalars(VTK_double, 1);

//    int* dims = imageData->GetDimensions();
//    for (int y = 0; y < dims[1]; y++)
//    {
//        for (int x = 0; x < dims[0]; x++)
//        {
//            double* pixel = static_cast<double*>(imageData->GetScalarPointer(x,y,0));
//            pixel[0] = (double)image[x][y][60];

//        }
//    }

//    vtkSmartPointer<vtkImageMapper> mapper = vtkSmartPointer<vtkImageMapper>::New();
//    mapper->SetInputData(imageData);
//    mapper->SetColorWindow(0);
//    mapper->SetColorLevel(0);

//    vtkSmartPointer<vtkActor2D> actor = vtkSmartPointer<vtkActor2D>::New();
//    actor->SetMapper(mapper);

//    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
//    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
//    renderWindow->AddRenderer(renderer);
//    renderWindow->SetSize(dims[0], dims[1]);
//    renderer->AddActor(actor);
//    renderer->ResetCamera();

//    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
//    renderWindowInteractor->SetRenderWindow(renderWindow);
//    renderWindow->Render();
//    renderWindowInteractor->Start();
}
