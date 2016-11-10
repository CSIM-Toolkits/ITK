#ifndef __itkLogisticContrastEnhancementImageFilter_hxx
#define __itkLogisticContrastEnhancementImageFilter_hxx
#include "itkLogisticContrastEnhancementImageFilter.h"


//Threshold methods
#include <itkMaximumEntropyThresholdImageFilter.h>
#include <itkOtsuThresholdImageFilter.h>
#include <itkIsoDataThresholdImageFilter.h>
#include <itkRenyiEntropyThresholdImageFilter.h>
#include <itkMomentsThresholdImageFilter.h>
#include <itkYenThresholdImageFilter.h>


#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkStatisticsImageFilter.h>
#include <itkSigmoidImageFilter.h>

#include <math.h>

using namespace std;

namespace itk
{
template< typename TInput, typename TOutput>
LogisticContrastEnhancementImageFilter< TInput, TOutput >
::LogisticContrastEnhancementImageFilter()
{
    this->m_FlipObjectArea=false;
    this->m_Alpha=0.0;
    this->m_Beta=0.0;
    this->m_MaximumOutput=1.0;
    this->m_MinimumOutput=0.0;
    this->m_Tolerance=1;
    this->m_ManualTolerance=false;
    this->m_NumberOfBins=128;
    this->m_ThresholdMethod=1;
}

template< typename TInput, typename TOutput >
void
LogisticContrastEnhancementImageFilter< TInput, TOutput >
::GenerateData()
{
    checkTolerance(m_Tolerance);

    //Input image
    typename InputImageType::ConstPointer input = this->GetInput();
    //Output image
    typename OutputImageType::Pointer output = this->GetOutput();
    output->SetRegions(input->GetBufferedRegion());
    output->Allocate();

    //Image threshold
    typedef itk::MaximumEntropyThresholdImageFilter<InputImageType, InputImageType>  MaxEntropyThresholdType;
    typedef itk::OtsuThresholdImageFilter<InputImageType, InputImageType>  OstuThresholdType;
    typedef itk::RenyiEntropyThresholdImageFilter<InputImageType, InputImageType>  RenyiThresholdType;
    typedef itk::IsoDataThresholdImageFilter<InputImageType, InputImageType>  IsoDataThresholdType;
    typedef itk::MomentsThresholdImageFilter<InputImageType, InputImageType>  MomentsThresholdType;
    typedef itk::YenThresholdImageFilter<InputImageType, InputImageType>  YenThresholdType;

    typename MaxEntropyThresholdType::Pointer maxEnThreshold = MaxEntropyThresholdType::New();
    typename OstuThresholdType::Pointer otsuThreshold = OstuThresholdType::New();

    typename RenyiThresholdType::Pointer renyiThreshold = RenyiThresholdType::New();
    typename IsoDataThresholdType::Pointer isoDataThreshold = IsoDataThresholdType::New();
    typename MomentsThresholdType::Pointer momentsThreshold = MomentsThresholdType::New();
    typename YenThresholdType::Pointer yenThreshold = YenThresholdType::New();

    double thr=0.0;
    switch (m_ThresholdMethod) {
    case MAXENTROPY:
        maxEnThreshold->SetInput(input);
        maxEnThreshold->SetInsideValue(0);
        maxEnThreshold->SetOutsideValue(1);
        maxEnThreshold->SetNumberOfHistogramBins(m_NumberOfBins);
        maxEnThreshold->Update();
        thr=maxEnThreshold->GetThreshold();
        break;
    case OTSU:
        otsuThreshold->SetInput(input);
        otsuThreshold->SetInsideValue(0);
        otsuThreshold->SetOutsideValue(1);
        otsuThreshold->SetNumberOfHistogramBins(m_NumberOfBins);
        otsuThreshold->Update();
        thr=otsuThreshold->GetThreshold();
        break;
    case RENYI:
        renyiThreshold->SetInput(input);
        renyiThreshold->SetInsideValue(0);
        renyiThreshold->SetOutsideValue(1);
        renyiThreshold->SetNumberOfHistogramBins(m_NumberOfBins);
        renyiThreshold->Update();
        thr=renyiThreshold->GetThreshold();
        break;
    case MOMENTS:
        momentsThreshold->SetInput(input);
        momentsThreshold->SetInsideValue(0);
        momentsThreshold->SetOutsideValue(1);
        momentsThreshold->SetNumberOfHistogramBins(m_NumberOfBins);
        momentsThreshold->Update();
        thr=momentsThreshold->GetThreshold();
        break;
    case ISODATA:
        isoDataThreshold->SetInput(input);
        isoDataThreshold->SetInsideValue(0);
        isoDataThreshold->SetOutsideValue(1);
        isoDataThreshold->SetNumberOfHistogramBins(m_NumberOfBins);
        isoDataThreshold->Update();
        thr=isoDataThreshold->GetThreshold();
        break;
    case YEN:
        yenThreshold->SetInput(input);
        yenThreshold->SetInsideValue(0);
        yenThreshold->SetOutsideValue(1);
        yenThreshold->SetNumberOfHistogramBins(m_NumberOfBins);
        yenThreshold->Update();
        thr=yenThreshold->GetThreshold();
        break;
    default:
        std::cout<<"ERROR: The threshold method is not valid! Choose the options available in the ThresholdMethod enumeration."<<std::endl;
        exit(EXIT_FAILURE);
        break;
    }

    //Set Beta
    typedef itk::StatisticsImageFilter<InputImageType> StatisticsType;
    typename StatisticsType::Pointer imageStatistics = StatisticsType::New();
    imageStatistics->SetInput(input);
    imageStatistics->Update();

    double beta=0.0;
    if (m_FlipObjectArea) {
        beta = thr/2.0;
    }else{
        beta = ((imageStatistics->GetMaximum()-thr)/2.0)+thr;
    }


    //Adjust automatic tolerance
    if (m_ManualTolerance) {
        //Set Alpha
        double alpha=0.0;
        if (m_FlipObjectArea) {
            alpha=((-1)*(thr)+beta)/(log((100.0-static_cast<double>(m_Tolerance))/static_cast<double>(m_Tolerance)));
        }else{
            alpha=((-1)*(thr)+beta)/(log((100.0-static_cast<double>(m_Tolerance))/static_cast<double>(m_Tolerance)));
        }

        //Apply sigmoid on input image
        typedef itk::SigmoidImageFilter<InputImageType, OutputImageType> SigmoidFilterType;
        typename SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();
        sigmoid->SetInput(input);
        sigmoid->SetOutputMinimum(m_MinimumOutput);
        sigmoid->SetOutputMaximum(m_MaximumOutput);
        sigmoid->SetAlpha(alpha);
        sigmoid->SetBeta(beta);
        sigmoid->Update();

        itk::ImageRegionConstIterator<TInput> sigmoidIterator(sigmoid->GetOutput(), sigmoid->GetOutput()->GetBufferedRegion());
        itk::ImageRegionIterator<TOutput> outputIterator(output, output->GetBufferedRegion());

        sigmoidIterator.IsAtBegin();
        outputIterator.IsAtBegin();
        while (!sigmoidIterator.IsAtEnd()) {
            outputIterator.Set(sigmoidIterator.Get());
            ++sigmoidIterator;
            ++outputIterator;
        }

        //Output the (alpha,beta) parameters
        m_Alpha=alpha;
        m_Beta=beta;
    }else{

        //Set Alpha
        double alpha=0.0;
        if (m_FlipObjectArea) {
            alpha=((-1)*(imageStatistics->GetMaximum())+beta)/(log(0.99/(1.0-0.99)));
        }else{
            alpha=((-1)*(imageStatistics->GetMaximum())+beta)/(log((1.0-0.99)/0.99));
        }

        //Apply sigmoid on input image
        typedef itk::SigmoidImageFilter<InputImageType, OutputImageType> SigmoidFilterType;
        typename SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();
        sigmoid->SetInput(input);
        sigmoid->SetOutputMinimum(m_MinimumOutput);
        sigmoid->SetOutputMaximum(m_MaximumOutput);
        sigmoid->SetAlpha(alpha);
        sigmoid->SetBeta(beta);
        sigmoid->Update();

        itk::ImageRegionConstIterator<TInput> sigmoidIterator(sigmoid->GetOutput(), sigmoid->GetOutput()->GetBufferedRegion());
        itk::ImageRegionIterator<TOutput> outputIterator(output, output->GetBufferedRegion());

        sigmoidIterator.IsAtBegin();
        outputIterator.IsAtBegin();
        while (!sigmoidIterator.IsAtEnd()) {
            outputIterator.Set(sigmoidIterator.Get());
            ++sigmoidIterator;
            ++outputIterator;
        }

        //Output the (alpha,beta) parameters
        m_Alpha=alpha;
        m_Beta=beta;
    }
}

template<typename TInput, typename TOutput>
void
LogisticContrastEnhancementImageFilter<TInput, TOutput>
::checkTolerance(char tolerance) {
    if (tolerance > 100) {
        std::cout<<"Tolerance is out of bound. Set in percentual (0 to 100)"<<std::endl;
        exit(EXIT_FAILURE);
    }
}

template<typename TInput, typename TOutput>
double
LogisticContrastEnhancementImageFilter<TInput, TOutput>
::sigmoid(double x, double alpha, double beta) {
    return 1/(1+std::exp((x-beta)/alpha));
}

} // end namespace itk

#endif
