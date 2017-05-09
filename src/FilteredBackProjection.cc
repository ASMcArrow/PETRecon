#include "FilteredBackProjection.hh"
#include "PETParameters.hh"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <iostream>

#include <fftw3.h>

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
    MaxBrightness = 0;
    Parameters = parameters;

    // Initialization of the image = TanBinNum x TanBinNum x Slices matrix
    Image = new double**[Parameters->GetTanBinNum()];
    for(int i = 0; i < Parameters->GetTanBinNum(); i++)
    {
        Image[i] = new double*[Parameters->GetTanBinNum()];
        for(int j = 0;j < Parameters->GetTanBinNum(); j++)
        {
            Image[i][j] = new double[Parameters->GetTotPlanes()];
            for(int k = 0; k < Parameters->GetRingNum(); k++)
            {
                Image[i][j][k] = 0.0;
            }
        }
    }

    Offset = ((double)Parameters->GetTanBinNum()-1.0)/2.0;  // Position of the center of first bin
    Costheta = new double[Parameters->GetAngBinNum()];
    Sintheta = new double[Parameters->GetAngBinNum()];

    XArray = new double*[Parameters->GetTanBinNum()];
    for(int i = 0; i < Parameters->GetTanBinNum(); i++)
    {
        XArray[i] = new double[Parameters->GetAngBinNum()];
    }

    YArray = new double*[Parameters->GetTanBinNum()];
    for(int i = 0; i < Parameters->GetTanBinNum(); i++)
    {
        YArray[i] = new double[Parameters->GetAngBinNum()];
    }

    for(int j = 0; j < Parameters->GetAngBinNum(); j++)
    {
        double theta = j*M_PI/(double)Parameters->GetAngBinNum();
        Costheta[j] = cos(theta);
        Sintheta[j] = sin(theta);
    }

    // Compute the coordinates of pixels of the image
    for(int i = 0; i < Parameters->GetTanBinNum(); i++)
    {
        double step = ((double)i-(((double)Parameters->GetTanBinNum()-1.0)/2.0));
        for(int j = 0; j < Parameters->GetAngBinNum(); j++)
        {
            XArray[i][j] = step*Costheta[j];
            YArray[i][j] = step*Sintheta[j];
        }
    }

    SortToSinogram2D(sinogram3D);
    CalculateFilter();
    FFT();

    // VTK visualization
    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData->SetExtent(0, Parameters->GetTanBinNum()-1, 0, Parameters->GetTanBinNum()-1, 0, 0);
    imageData->SetSpacing(1.0, 1.0, 1.0);
    imageData->SetOrigin(0.0, 0.0, 0.0);
    imageData->AllocateScalars(VTK_DOUBLE, 1);

    int* dims = imageData->GetDimensions();
    for (int y = 0; y < dims[1]; y++)
    {
        for (int x = 0; x < dims[0]; x++)
        {
            double* pixel = static_cast<double*>(imageData->GetScalarPointer(x,y,0));
            pixel[0] = (double)Image[x][y][31];
        }
    }

    vtkSmartPointer<vtkImageMapper> mapper = vtkSmartPointer<vtkImageMapper>::New();
    mapper->SetInputData(imageData);
    mapper->SetColorWindow(MaxBrightness);
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

void FilteredBackProjection::SortToSinogram2D(unsigned short*** sinogram3D)
{
    // We convert 3D sinogram to 2D sinogram with slices along the diagonals of michelogram
    int numOfPlanes = Parameters->GetTotPlanes();
    int numOfAxialPlanes = Parameters->GetTotAxialPlanes();
    int numOfSlices = Parameters->GetRingNum();
    int* planesNorm = new int[numOfPlanes];

    Sinogram2D = new double**[Parameters->GetTanBinNum()];
    for(int i = 0; i < Parameters->GetTanBinNum(); i++)
    {
        Sinogram2D[i] = new double*[Parameters->GetAngBinNum()];
        for(int j = 0; j < Parameters->GetAngBinNum(); j++)
        {
            Sinogram2D[i][j] = new double[numOfPlanes];

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
                    Sinogram2D[i][j][k] = (double)(Sinogram2D[i][j][k]/(double)planesNorm[k]);
                }
            }
        }
    }
}

void FilteredBackProjection::CalculateFilter()
{
    // Zero padding for FFT: searching for P=2^n with P the first integer that is >= nukber of tangential bins
    int n = (int)ceil(log(Parameters->GetTanBinNum())/log(2));
    inputFFTSize = (int)pow(2,n);
    // Input and output arrays are of different sizes and types: the input is n real numbers, while the output is n/2+1 complex numbers
    outputFFTSize = inputFFTSize/2+1;

    // Calculate Hamming filter
    Filter = new double[outputFFTSize];
    double alpha = 0.75;
    double cutNFrequency = 0.75;
    double pixelSize = 1;
    // The cut-off frequency is expressed as a fraction of the Nq frequency.
    // The cut-off frequency varies typically from 0.2 to 1.0 times the Nq frequency.
    // It defines the frequency above which the noise is eliminated.
    double cutFrequency = cutNFrequency/(2*pixelSize);  // Nyquist frecuency = 1 / (2*pixelSize)
    for(int sampleNum = 1; sampleNum < outputFFTSize; sampleNum++) // We calculate frequency for each sample in the frequency domain
    {
        double harmonicFrequency = (double)sampleNum/(pixelSize*(double)inputFFTSize);
        if(harmonicFrequency < cutFrequency)
        {
            Filter[sampleNum] = 1.0; //harmonicFrequency*(alpha + (1.0-alpha)*cos(M_PI*sampleNum/cutFrequency));
        } else
        {
            Filter[sampleNum] = 0.0;
        }
    }
}

void FilteredBackProjection::FFT()
{
    double* inputData = (double*) fftw_malloc(sizeof(double)*inputFFTSize);
    fftw_complex* outputData = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*outputFFTSize);

    fftw_plan planToFourier = fftw_plan_dft_r2c_1d(inputFFTSize, inputData, outputData, FFTW_ESTIMATE);
    fftw_plan planFromFourier = fftw_plan_dft_c2r_1d(inputFFTSize, outputData, inputData, FFTW_ESTIMATE);

    for(int k = 0; k < Parameters->GetTotPlanes(); k++)
    {
        for(int j = 0; j < Parameters->GetAngBinNum(); j++)
        {
            for(int i = 0; i < Parameters->GetTanBinNum(); i++)
            {
                inputData[i] = (double)Sinogram2D[i][j][k];
            }

            for(int i =  Parameters->GetTanBinNum(); i < inputFFTSize; i++)
            {
                inputData[i] = 0.0;
            }

            fftw_execute(planToFourier);

            outputData[0][0] = 0.0;            // Handle zero freq. specially
            outputData[0][1] = 0.0;
            //outputData[inputFFTSize/2][0] = outputData[inputFFTSize/2][0] * Filter[inputFFTSize/2]; // Handle half sampling freq.
            //outputData[inputFFTSize/2][1] = outputData[inputFFTSize/2][1] * Filter[inputFFTSize/2];

            for(int i = 1; i < outputFFTSize; i++)
            {
                // cout << "i = " << i << " " << (double)outputData[i][0] << " + i*" << (double)outputData[i][1] << endl;
                outputData[i][0] = outputData[i][0]*Filter[i];
                outputData[i][1] = outputData[i][1]*Filter[i];
            }


            // Compute Inverse FFT, fftw inverse values are unnormalized (scaled by n = ouputFFTSize)
            fftw_execute(planFromFourier);

            for(int i = 0; i < Parameters->GetTanBinNum(); i++)
            {
                Sinogram2D[i][j][k] = inputData[i]/inputFFTSize;
            }
        }

        BackProject(k);
    }
}

void FilteredBackProjection::BackProject(int slice)
{
    for(int x = 0; x < Parameters->GetTanBinNum(); x++)
    {
        for(int y = 0; y < Parameters->GetTanBinNum(); y++)
        {
            double sum = 0.0;
            for(int angle = 0; angle < Parameters->GetAngBinNum(); angle++)
            {
                double dist = XArray[x][angle] + YArray[y][angle];
                int bin = (int)(floor(dist) + Parameters->GetTanBinNum()/2);
                double w = dist - floor(dist);

                // Using linear interpolation
                if (dist >= -Parameters->GetTanBinNum()/2 && dist < Parameters->GetTanBinNum()/2)
                {
                    if (bin >= 0 && bin < Parameters->GetTanBinNum()-1)
                    {
                        sum = sum + (1-w)*Sinogram2D[bin][angle][slice] + w*Sinogram2D[bin+1][angle][slice];
                    }
                    else if (bin == Parameters->GetTanBinNum()-1)
                        sum = sum + (1-w)*Sinogram2D[bin][angle][slice];
                }
            }

            Image[x][y][slice] = (double)(sum*M_PI/(double)Parameters->GetAngBinNum());
            if (Image[x][y][slice] >= MaxBrightness)
                MaxBrightness = Image[x][y][slice];
        }
    }
}
