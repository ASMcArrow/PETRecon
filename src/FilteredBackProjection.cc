#include "FilteredBackProjection.hh"
#include "PETParameters.hh"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <iostream>

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

using namespace std;

FilteredBackProjection::FilteredBackProjection(unsigned short*** sinogram3D, PETParameters* parameters)
{
    Parameters = parameters;
    SortToSinogram2D(sinogram3D);
    // ApplyFilter();
    Reconstruct();
}

void FilteredBackProjection::SortToSinogram2D(unsigned short*** sinogram3D)
{
    // We convert 3D sinogram to 2D sinogram with slices along the diagonals of michelogram
    int numOfPlanes = Parameters->GetTotPlanes();
    int numOfAxialPlanes = Parameters->GetTotAxialPlanes();
    int numOfSlices = Parameters->GetRingNum();
    int* planesNorm = new int[numOfPlanes];

    Sinogram2D = new float**[Parameters->GetTanBinNum()];
    for(int i = 0; i < Parameters->GetTanBinNum(); i++)
    {
        Sinogram2D[i] = new float*[Parameters->GetAngBinNum()];
        for(int j = 0; j < Parameters->GetAngBinNum(); j++)
        {
            Sinogram2D[i][j] = new float[numOfPlanes];

            for(int k = 0; k < numOfPlanes; k++)
            {
                Sinogram2D[i][j][k] = 0;
                planesNorm[k] = 0;
            }

            for(int k = 0; k < numOfAxialPlanes; k++)
            {
                int plane2D = (int)(k%numOfSlices + k/numOfSlices);

                if (abs((k/numOfSlices - k%numOfSlices)) <= Parameters->GetMaxRingDiff())
                {
                    planesNorm[plane2D]++;
                    Sinogram2D[i][j][plane2D] = (int)sinogram3D[i][j][k] + Sinogram2D[i][j][plane2D];
                }
            }

            // Normalization
            for(int k = 0; k < numOfPlanes; k++)
            {
                if (planesNorm[k] > 1)
                {
                    Sinogram2D[i][j][k] = (float)(Sinogram2D[i][j][k]/(float)planesNorm[k]);
                }
            }
        }
    }
}

void FilteredBackProjection::Reconstruct()
{
    int numOfPlanes = Parameters->GetTotPlanes();

    // Initialization of the image = TanBinNum x TanBinNum x Slices matrix
    Image = new float**[Parameters->GetTanBinNum()];
    for(int i = 0; i < Parameters->GetTanBinNum(); i++)
    {
        Image[i] = new float*[Parameters->GetTanBinNum()];
        for(int j = 0;j < Parameters->GetTanBinNum(); j++)
        {
            Image[i][j] = new float[numOfPlanes];
            for(int k = 0; k < Parameters->GetRingNum(); k++)
            {
                Image[i][j][k] = 0.0;
            }
        }
    }

    double offset = -((double)Parameters->GetTanBinNum()-1)/2;  // Position of the center of first bin
    double* costheta = new double[Parameters->GetAngBinNum()];
    double* sintheta = new double[Parameters->GetAngBinNum()];

    double** xArray = new double*[Parameters->GetTanBinNum()];
    for(int i = 0; i < Parameters->GetTanBinNum(); i++)
    {
        xArray[i] = new double[Parameters->GetAngBinNum()];
    }

    double** yArray = new double*[Parameters->GetTanBinNum()];
    for(int i = 0; i < Parameters->GetTanBinNum(); i++)
    {
        yArray[i] = new double[Parameters->GetAngBinNum()];
    }

    for(int j = 0; j  < Parameters->GetAngBinNum(); j++)
    {
        double theta = j*M_PI/(double)Parameters->GetAngBinNum();
        costheta[j] = cos(theta);
        sintheta[j] = sin(theta);
    }

    // Compute the coordinates of pixels of the image
    for(int i = 0; i < Parameters->GetTanBinNum(); i++)
    {
        double step = (i-(Parameters->GetTanBinNum()-1)/2);
        for(int j = 0; j < Parameters->GetAngBinNum(); j++)
        {
            xArray[i][j]=step*costheta[j];
            yArray[i][j]=step*sintheta[j]-offset;
        }
    }

    // Back projection
    for (int slice = 0; slice < numOfPlanes; slice++)
    {

        for(int x = 0; x < Parameters->GetTanBinNum(); x++)
        {
            for(int y = 0; y < Parameters->GetTanBinNum(); y++)
            {
                float sum = 0.0;
                for(int angle = 0; angle < Parameters->GetAngBinNum(); angle++)
                {
                    int dist = xArray[x][angle]+yArray[y][angle];
                    float w = dist - floor(dist);

                    // Using linear interpolation
                    if (dist >= 0 && dist < Parameters->GetTanBinNum()-1)
                    {
                        sum = sum + (1-w)*Sinogram2D[dist][angle][slice]+w*Sinogram2D[dist+1][angle][slice];
                    } else if (dist == Parameters->GetTanBinNum()-1)
                    {
                        sum = sum + (1-w)*Sinogram2D[dist][angle][slice];
                    } else if (dist == -1)
                    {
                        sum = sum + w*Sinogram2D[0][angle][slice];
                    }
                }

                Image[x][y][slice] = (float)(sum*M_PI/(double)Parameters->GetAngBinNum());
            }
        }
    }

    // VTK visualization
    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData->SetExtent(0, Parameters->GetTanBinNum()-1, 0, Parameters->GetTanBinNum()-1, 0, 0);
    imageData->SetSpacing(1.0, 1.0, 1.0);
    imageData->SetOrigin(0.0, 0.0, 0.0);
    imageData->AllocateScalars(VTK_FLOAT, 1);

    int* dims = imageData->GetDimensions();
    for (int y = 0; y < dims[1]; y++)
    {
        for (int x = 0; x < dims[0]; x++)
        {
            float* pixel = static_cast<float*>(imageData->GetScalarPointer(x,y,0));
            pixel[0] = (float)Image[x][y][32];
        }
    }

    vtkSmartPointer<vtkImageMapper> mapper = vtkSmartPointer<vtkImageMapper>::New();
    mapper->SetInputData(imageData);
    mapper->SetColorWindow(0);
    mapper->SetColorLevel(0);

    vtkSmartPointer<vtkActor2D> actor = vtkSmartPointer<vtkActor2D>::New();
    actor->SetMapper(mapper);

    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderWindow->AddRenderer(renderer);
    renderWindow->SetSize(dims[0], dims[1]);
    renderer->AddActor(actor);
    renderer->ResetCamera();

    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);
    renderWindow->Render();
    renderWindowInteractor->Start();
}
