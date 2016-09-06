#ifndef __itkSigmoidParametersOptimizationImageCalculator_hxx
#define __itkSigmoidParametersOptimizationImageCalculator_hxx
#include "itkSigmoidParametersOptimizationImageCalculator.h"

#include <itkConstNeighborhoodIterator.h>
#include <itkNeighborhoodIterator.h>
#include <itkMaskImageFilter.h>
#include <itkSigmoidImageFilter.h>
#include <itkImageToHistogramFilter.h>

//Threshold methods
#include <itkMaximumEntropyThresholdImageFilter.h>


#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkStatisticsImageFilter.h>
#include <math.h>

namespace itk
{
template< typename TImage >
SigmoidParametersOptimizationImageCalculator< TImage >
::SigmoidParametersOptimizationImageCalculator()
{
    this->SetNumberOfRequiredInputs(1);
    this->m_KernelRadius=1;
    this->m_MaximumBeta=1.0;
    this->m_MinimumBeta=0.0;
    this->m_MaximumAlpha=1.0;
    this->m_MinimumAlpha=0.0;
    this->m_NumberOfBins=128;
    this->m_NumberOfTrials=20;
}

template< typename TImage >
void
SigmoidParametersOptimizationImageCalculator< TImage >
::SetInput(const TImage *input){
    this->SetNthInput(0, const_cast<TImage*>(input));
}

template< typename TImage >
void
SigmoidParametersOptimizationImageCalculator< TImage >
::GenerateData()
{
    //Input image
    //    typename TImage::ConstPointer input = static_cast<const TImage*>(this->ProcessObject::GetInput(0));
    typename TImage::ConstPointer input = this->GetInput();

    //Iterators
    typedef itk::ConstNeighborhoodIterator<TImage>      ConstNeighborIteratorType;
    typedef itk::NeighborhoodIterator<TImage>           NeighborIteratorType;

    //Energy image
    typename TImage::Pointer energy = TImage::New();
    energy->SetRegions(input->GetBufferedRegion());
    energy->Allocate();

    itk::ImageRegionConstIterator<TImage> inputIterator(input, input->GetBufferedRegion());
    itk::ImageRegionIterator<TImage> energyIterator(energy, energy->GetBufferedRegion());
    while(!inputIterator.IsAtEnd())
    {
        energyIterator.Set(inputIterator.Get());
        ++inputIterator;
        ++energyIterator;
    }

    //Set kernel radius
    typename TImage::SizeType radius;
    for (int imgDim = 0; imgDim < TImageDimension; ++imgDim) {
        radius[imgDim] = m_KernelRadius;
    }

    NeighborIteratorType energyIt(radius, energy, energy->GetBufferedRegion());
    ConstNeighborIteratorType inputIt(radius, input, input->GetBufferedRegion());

    //Calculating Energy image
    double x2 = 0.0;
    energyIt.GoToBegin();
    inputIt.GoToBegin();
    while (!energyIt.IsAtEnd()) {
        for (int idx = 0; idx < pow(inputIt.GetSize()[0],TImageDimension); ++idx) {
            x2=x2+pow(static_cast<double>(inputIt.GetPixel(idx)),2.0);
        }
//                energyIt.SetCenterPixel(10.0*log10(x2)/pow(inputIt.GetSize()[0],TImageDimension));
        energyIt.SetCenterPixel(x2/pow(inputIt.GetSize()[0],TImageDimension));
        x2=0.0;
        ++energyIt;
        ++inputIt;
    }


    //Split the En and Esn images
    typedef itk::MaximumEntropyThresholdImageFilter<TImage, TImage> MaximumEntropyThresholdType;
    typename MaximumEntropyThresholdType::Pointer energyNImage = MaximumEntropyThresholdType::New();
    typename MaximumEntropyThresholdType::Pointer energySNImage = MaximumEntropyThresholdType::New();
    energyNImage->SetInput(energy);
    energyNImage->SetInsideValue(1);
    energyNImage->SetOutsideValue(0);
    energyNImage->SetNumberOfHistogramBins(128);

    energySNImage->SetInput(energy);
    energySNImage->SetInsideValue(0);
    energySNImage->SetOutsideValue(1);
    energySNImage->SetNumberOfHistogramBins(128);

    typedef itk::MaskImageFilter<TImage, TImage>        MaskFilterType;
    typename MaskFilterType::Pointer energyN = MaskFilterType::New();
    energyN->SetInput(energy);
    energyN->SetMaskImage(energyNImage->GetOutput());

    typedef itk::MaskImageFilter<TImage, TImage>        MaskFilterType;
    typename MaskFilterType::Pointer energySN = MaskFilterType::New();
    energySN->SetInput(energy);
    energySN->SetMaskImage(energySNImage->GetOutput());


    //Energy N and SN histograms
    const unsigned int MeasurementVectorSize = 1; // Grayscale
    const unsigned int binsPerDimension = m_NumberOfBins;

    typedef itk::Statistics::ImageToHistogramFilter< TImage >        ImageToHistogramFilterType;
    typename ImageToHistogramFilterType::HistogramType::MeasurementVectorType     lowerBound(binsPerDimension);
    lowerBound.Fill(0);

    typename ImageToHistogramFilterType::HistogramType::MeasurementVectorType     upperBound(binsPerDimension);
    upperBound.Fill(255);

    typename ImageToHistogramFilterType::HistogramType::SizeType                  size(MeasurementVectorSize);
    size.Fill(binsPerDimension);

    typename ImageToHistogramFilterType::Pointer imageToHistEn =  ImageToHistogramFilterType::New();
    imageToHistEn->SetInput( energyN->GetOutput() );
    imageToHistEn->SetHistogramBinMinimum( lowerBound );
    imageToHistEn->SetHistogramBinMaximum( upperBound );
    imageToHistEn->SetHistogramSize( size );
    typename ImageToHistogramFilterType::Pointer imageToHistEsn =  ImageToHistogramFilterType::New();
    imageToHistEsn->SetInput( energySN->GetOutput() );
    imageToHistEsn->SetHistogramBinMinimum( lowerBound );
    imageToHistEsn->SetHistogramBinMaximum( upperBound );
    imageToHistEsn->SetHistogramSize( size );

    imageToHistEn->Update();
    imageToHistEsn->Update();

    typename ImageToHistogramFilterType::HistogramType* histogramEn = imageToHistEn->GetOutput();
    typename ImageToHistogramFilterType::HistogramType* histogramEsn = imageToHistEsn->GetOutput();

    double J[m_NumberOfTrials][m_NumberOfTrials];
    for (int beta = 0; beta < m_NumberOfTrials; ++beta) {
        for (int alpha = 0; alpha < m_NumberOfTrials; ++alpha) {
            //P(alpha,beta)= mean(sigmoid(En))
            double P = noiseFunction(histogramEn,alpha*((m_MaximumAlpha-m_MinimumAlpha)/m_NumberOfTrials) + m_MinimumAlpha, beta*((m_MaximumBeta-m_MinimumBeta)/m_NumberOfTrials) + m_MinimumBeta);

            //Dnonlinear(alpha,beta)=mean( (Aj*Esn+Bj - sigmoid(Esn))^2 )/mean(Esn^2)
            double Dnon = nonLinearFunction(histogramEsn,alpha*((m_MaximumAlpha-m_MinimumAlpha)/m_NumberOfTrials) + m_MinimumAlpha, beta*((m_MaximumBeta-m_MinimumBeta)/m_NumberOfTrials) + m_MinimumBeta);

            //V(alpha,beta)= V1 - V2
            double V = varianceFunction(histogramEsn,alpha*((m_MaximumAlpha-m_MinimumAlpha)/m_NumberOfTrials) + m_MinimumAlpha, beta*((m_MaximumBeta-m_MinimumBeta)/m_NumberOfTrials) + m_MinimumBeta);

            J[alpha][beta]=Dnon+P-V;
        }
    }

    //Minimum of J(alpha,beta) matrix
    double minJ=J[1][1];
    double optAlpha=0;
    double optBeta=0;
    for (int beta = 1; beta < m_NumberOfTrials; ++beta) {
        for (int alpha = 1; alpha < m_NumberOfTrials; ++alpha) {
            if (J[alpha][beta]<minJ) {
                optAlpha=alpha;
                optBeta=beta;
                minJ=J[alpha][beta];
            }
        }
    }

    //Set optimum Alpha and Beta
    m_OptimumAlpha=optAlpha*((m_MaximumAlpha - m_MinimumAlpha)/m_NumberOfTrials);
    m_OptimumBeta=optBeta*((m_MaximumBeta - m_MinimumBeta)/m_NumberOfTrials);
}

template<typename TImage>
double
SigmoidParametersOptimizationImageCalculator<TImage>
::nonLinearFunction(ImageToHistogramType histogram, double alpha, double beta) {
    //Dnon(alpha,beta)=D1/D2;
    //D1=mean( (Aj*Esn+Bj - sigmoid(Esn))^2 )
    //D2=mean(Esn^2)
    double D1=0.0,D2=0.0;
    for (int idx = 0; idx < histogram->GetSize()[0]; ++idx) {
        D1=D1+((Aj(histogram,alpha,beta)*histogram->GetFrequency(idx)+Bj(histogram, alpha, beta)) - sigmoid(histogram->GetFrequency(idx),alpha,beta));
        D2=D2+pow(histogram->GetFrequency(idx),2.0);
    }

    return D1/D2;
}

template<typename TImage>
double
SigmoidParametersOptimizationImageCalculator<TImage>
::noiseFunction(ImageToHistogramType histogram, double alpha, double beta) {
    double P=0.0;
    for (int idx = 0; idx < histogram->GetSize()[0]; ++idx) {
        P=P+pow(sigmoid(histogram->GetFrequency(idx),alpha, beta),2.0);
    }

    return P/m_NumberOfBins;
}

template<typename TImage>
double
SigmoidParametersOptimizationImageCalculator<TImage>
::varianceFunction(ImageToHistogramType histogram, double alpha, double beta) {
    double V=0.0;
    double V1 = 0.0;
    for (int idx = 0; idx < histogram->GetSize()[0]; ++idx) {
        V1=V1+pow(sigmoid(histogram->GetFrequency(idx),alpha, beta),2.0);
    }
    V1=V1/m_NumberOfBins;
    double V2 = 0.0;
    for (int idx = 0; idx < histogram->GetSize()[0]; ++idx) {
        V2=V2+sigmoid(histogram->GetFrequency(idx),alpha, beta);
    }
    V2=pow(V2/m_NumberOfBins,2.0);
    V=V1-V2;

    return V;
}

template<typename TImage>
double
SigmoidParametersOptimizationImageCalculator<TImage>
::Aj(ImageToHistogramType histogram, double alpha, double beta){
    double mu=0.0,sigma=0.0, Aj=0.0;
    double A1,A2;

    //mu
    for (int idx = 0; idx < histogram->GetSize()[0]; ++idx) {
        mu=mu+histogram->GetFrequency(idx);
    }
    mu=mu/m_NumberOfBins;

    //sigma
    for (int idx = 0; idx < histogram->GetSize()[0]; ++idx) {
        sigma=sigma+pow(histogram->GetFrequency(idx),2.0);
    }
    sigma=(sigma/m_NumberOfBins)-pow(mu,2.0);

    //Aj=(1/sigma)*( A1 - A2 )
    //A1=mean(Esn*sigmoid(Esn))
    //A2=mu*mean(sigmoid(Esn))
    for (int idx = 0; idx < histogram->GetSize()[0]; ++idx) {
        A1=A1+(histogram->GetFrequency(idx)*sigmoid(histogram->GetFrequency(idx),alpha, beta));
        A2=A2+sigmoid(histogram->GetFrequency(idx),alpha,beta);
    }
    A1=A1/m_NumberOfBins;
    A2=mu*(A2/m_NumberOfBins);
    Aj=(1/sigma)*(A1-A2);

    return Aj;
}

template<typename TImage>
double
SigmoidParametersOptimizationImageCalculator<TImage>
::Bj(ImageToHistogramType histogram, double alpha, double beta){
    double mu=0.0, Bj=0.0;

    //mu
    for (int idx = 0; idx < histogram->GetSize()[0]; ++idx) {
        mu=mu+histogram->GetFrequency(idx);
    }
    mu=mu/m_NumberOfBins;

    //Bj=mean(sigmoid(Esn) - mu*Aj)
    for (int idx = 0; idx < histogram->GetSize()[0]; ++idx) {
        Bj=Bj+histogram->GetFrequency(idx);
    }
    Bj=Bj/m_NumberOfBins;

    return Bj - mu*Aj(histogram, alpha, beta);
}

template<typename TImage>
double
SigmoidParametersOptimizationImageCalculator<TImage>
::sigmoid(double x, double alpha, double beta) {
    return 1/(1+std::exp((x-beta)/alpha));
}

} // end namespace itk

#endif
