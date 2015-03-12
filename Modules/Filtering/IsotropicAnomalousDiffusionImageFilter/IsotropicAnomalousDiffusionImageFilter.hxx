#ifndef __itkIsotropicAnomalousDiffusionImageFilter_hxx
#define __itkIsotropicAnomalousDiffusionImageFilter_hxx
#include "IsotropicAnomalousDiffusionImageFilter.h"

#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkConstNeighborhoodIterator.h>
#include <itkLaplacianOperator.h>

namespace itk
{
template< typename TInputImage, typename TOutputImage >
IsotropicAnomalousDiffusionImageFilter< TInputImage, TOutputImage >
::IsotropicAnomalousDiffusionImageFilter()
{
    m_Q = 1.0;
    m_GeneralizedDiffusion = 1.0;
    m_Iterations = 1;
    m_TimeStep = (1.0 / std::pow(2.0,static_cast<double>(InputImageDimension)+1));
}

template< typename TInputImage, typename TOutputImage >
void
IsotropicAnomalousDiffusionImageFilter< TInputImage, TOutputImage >
::ThreadedGenerateData(const OutputImageRegionType &region, ThreadIdType )
{
    //Test stability for the discrete time iteration. This definition is passed by the TimeStep variable.
    TimeStepTestStability();

    typedef itk::ImageRegionIterator< OutputImageType>              IteratorType;
    typedef itk::ImageRegionConstIterator< InputImageType>          ConstIteratorType;
    typedef itk::ConstNeighborhoodIterator< OutputImageType >       NeighborIteratorType;
    typename InputImageType::ConstPointer input = this->GetInput();
    typename OutputImageType::Pointer output = this->GetOutput();
    typename InputImageType::Pointer aux = InputImageType::New();
    aux->SetRegions(region);
    aux->Allocate();

    //    Copy input image
    IteratorType outputIt(output, region);
    IteratorType auxIt(aux, region);
    ConstIteratorType inputIt(input, region);

    outputIt.GoToBegin();
    auxIt.GoToBegin();
    inputIt.GoToBegin();
    while (!outputIt.IsAtEnd()) {
        outputIt.Set(inputIt.Get());
        auxIt.Set(inputIt.Get());
        ++outputIt;
        ++auxIt;
        ++inputIt;
    }

    typedef itk::LaplacianOperator< OutputPixelType , OutputImageDimension > LaplacianOperatorType;
    LaplacianOperatorType laplaceOp;
    laplaceOp.SetRadius(InputImageDimension);
    laplaceOp.CreateOperator();
    NeighborIteratorType      laplaceIt(laplaceOp.GetRadius(), output,  region);

    double neighborAux = 0.0;
    for(int i=0; i<m_Iterations; i++){
        auxIt.GoToBegin();
        laplaceIt.GoToBegin();

        while( !auxIt.IsAtEnd() )
        {
            neighborAux = static_cast<InputPixelType>(0.0);

            for (unsigned int idx = 0; idx < pow(laplaceIt.GetSize()[0],OutputImageDimension); ++idx) {
                neighborAux += std::pow(laplaceIt.GetPixel(idx), 2.0 - m_Q)*laplaceOp.GetElement(idx);
            }
            auxIt.Set(neighborAux*this->GeneralizedDiffCurve()*m_TimeStep + laplaceIt.GetCenterPixel());

            ++auxIt;
            ++laplaceIt;
        }

        outputIt.GoToBegin();
        auxIt.GoToBegin();
        while (!outputIt.IsAtEnd()) {
            outputIt.Set(auxIt.Get());
            ++outputIt;
            ++auxIt;
        }
    }
}

template< typename TInputImage, typename TOutputImage >
void IsotropicAnomalousDiffusionImageFilter< TInputImage, TOutputImage >
::TimeStepTestStability()
{
    if ( m_TimeStep >  ( 1.0 / std::pow(2.0, static_cast< double >( InputImageDimension ) +1) ))
    {
        itkWarningMacro( << "Isotropic diffusion unstable time step: "
                         << m_TimeStep << std::endl
                         << "Stable time step for this image must be smaller than "
                         << 1.0 / std::pow( 2.0, static_cast< double >( InputImageDimension ) +1 ) );
    }
}

template< typename TInputImage, typename TOutputImage >
double IsotropicAnomalousDiffusionImageFilter<TInputImage, TOutputImage >
::GeneralizedDiffCurve()
{

    double d = 0.0;
    double alpha = (2.0 - m_Q)*(3.0 - m_Q);
    if (m_Q < 1.0) {
        //            d = 2*Math.pow(alpha, 2/(3-q))*Math.pow(Math.sqrt((1 - q) / Math.PI) * (gamma(1 + (1 / (1 - q))) / gamma(3 / 2 + (1 / (1 - q)))), ((2 - 2 * q) / (3 - q)));
        d = m_GeneralizedDiffusion * std::exp((-1.0) * (std::pow(m_Q - 1.0, 2.0)) / 0.08);
    } else if (m_Q >= 1.0 && m_Q < 2.0) {
        d = m_GeneralizedDiffusion * std::pow(alpha, 2.0 / (3.0 - m_Q)) * std::pow(sqrt((m_Q - 1.0) / M_PI) * (tgamma(1.0 / (m_Q - 1.0)) / tgamma((1.0 / (m_Q - 1.0)) - 1.0 / 2.0)), ((2.0 - 2.0 * m_Q) / (3.0 - m_Q)));
        //            d = percentD * Math.exp((-1) * (Math.pow(q - 1.0, 2.0)) / 0.4);
    } else if (m_Q == 1.0) {
        d = m_GeneralizedDiffusion;
    }

    return static_cast<InputPixelType>(d);
}


} // end namespace itk

#endif
