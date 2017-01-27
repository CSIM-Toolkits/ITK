#ifndef __itkAnisotropicAnomalousDiffusionImageFilter_hxx
#define __itkAnisotropicAnomalousDiffusionImageFilter_hxx
#include "itkAnisotropicAnomalousDiffusionImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkLaplacianOperator.h>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

namespace itk
{
template< typename TInputImage, typename TOutputImage >
AnisotropicAnomalousDiffusionImageFilter< TInputImage, TOutputImage >
::AnisotropicAnomalousDiffusionImageFilter()
{
    m_Q = 1.0;
    m_Conductance = 1.0;
    m_Iterations = 1;
    m_GeneralizedDiffusion=1.0;
    m_TimeStep = (1.0 / pow(2.0,static_cast<double>(InputImageDimension) + 1));
}

template<typename TInputImage, typename TOutputImage>
void
AnisotropicAnomalousDiffusionImageFilter< TInputImage, TOutputImage >
::ThreadedGenerateData(const OutputImageRegionType & region, ThreadIdType )
{

    typedef itk::ImageRegionIterator< OutputImageType>              IteratorType;
    typedef itk::ImageRegionConstIterator< InputImageType>          ConstIteratorType;
    typename InputImageType::ConstPointer input = this->GetInput();
    typename OutputImageType::Pointer output = this->GetOutput();
    typename InputImageType::Pointer aux = InputImageType::New();
    aux->SetRegions(region);
    aux->Allocate();

    //Test the algorithm stability
    TimeStepTestStability();

    //Copy input image
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
    itk::Size<OutputImageDimension> radius;
    radius.Fill(1);
    laplaceOp.CreateToRadius(radius);
    NeighborIteratorType      laplaceIt(laplaceOp.GetRadius(), output,  region);

    InputPixelType neighborAux = 0.0, center_value, next_value, edgeController, derivative, diffusionPixelCorrected;
    for(int i=0; i<m_Iterations; i++){
        auxIt.GoToBegin();
        laplaceIt.GoToBegin();

        while( !auxIt.IsAtEnd() )
        {
            if (auxIt.Get()!=static_cast<InputPixelType>(0)) {

                neighborAux = static_cast<InputPixelType>(0);
                for (unsigned int idx = 0; idx < pow(laplaceIt.GetSize()[0],static_cast< double >(InputImageDimension)); ++idx) {
                    center_value = laplaceIt.GetCenterPixel();
                    next_value =  laplaceIt.GetPixel(idx);
                    edgeController = this->EdgeWeightedController(next_value, center_value);
                    derivative = pow(next_value, 2.0 - m_Q) - pow(center_value, 2.0 - m_Q);
                    neighborAux += (idx%2==0?1.0:0.5)*edgeController*derivative;
                }
                diffusionPixelCorrected = neighborAux*static_cast<InputPixelType>(m_TimeStep) + laplaceIt.GetCenterPixel();
                if(diffusionPixelCorrected>static_cast<InputPixelType>(0)){
                    auxIt.Set(diffusionPixelCorrected);
                }else{
                    auxIt.Set(static_cast<InputPixelType>(meanNeighbors(laplaceIt)));
                }

                ++auxIt;
                ++laplaceIt;
            }else{
                ++auxIt;
                ++laplaceIt;
            }
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

template<typename TInputImage, typename TOutputImage>
double AnisotropicAnomalousDiffusionImageFilter<TInputImage, TOutputImage >
::meanNeighbors(NeighborIteratorType neighbors){
    up = neighbors.GetPixel(2);
    down = neighbors.GetPixel(8);
    left = neighbors.GetPixel(4);
    right = neighbors.GetPixel(6);

    return ((up+down+left+right)/static_cast<InputPixelType>(4));
}

template< typename TInputImage, typename TOutputImage>
double AnisotropicAnomalousDiffusionImageFilter<TInputImage, TOutputImage >
::EdgeWeightedController(InputPixelType idxValue, InputPixelType centerValue)
{
    return  GeneralizedDiffCurve()*exp((-1.0)*pow(static_cast<InputPixelType>(abs((idxValue - centerValue)))/static_cast<InputPixelType>(m_Conductance),  2.0));
}

template< typename TInputImage, typename TOutputImage >
void AnisotropicAnomalousDiffusionImageFilter< TInputImage, TOutputImage >
::TimeStepTestStability()
{
    if ( m_TimeStep >  ( 1.0 / pow(2.0, static_cast< double >( InputImageDimension ) +1) ))
    {
        itkWarningMacro( << "Anisotropic diffusion unstable time step: "
                         << m_TimeStep << endl
                         << "Stable time step for this image must be smaller than "
                         << 1.0 / pow( 2.0, static_cast< double >( InputImageDimension ) +1) );
    }
}

template< typename TInputImage, typename TOutputImage >
double AnisotropicAnomalousDiffusionImageFilter<TInputImage, TOutputImage >
::GeneralizedDiffCurve()
{

    float d = 0.0;
    float alpha = (2.0 - m_Q)*(3.0 - m_Q);
    if (m_Q < 1.0) {
        //            d = 2*Math.pow(alpha, 2/(3-q))*Math.pow(Math.sqrt((1 - q) / Math.PI) * (gamma(1 + (1 / (1 - q))) / gamma(3 / 2 + (1 / (1 - q)))), ((2 - 2 * q) / (3 - q)));
        d = m_GeneralizedDiffusion * exp((-1.0) * (pow(m_Q - 1.0, 2.0)) / 0.08);
    } else if (m_Q >= 1.0 && m_Q < 2.0) {
        d = m_GeneralizedDiffusion * pow(alpha, 2.0 / (3.0 - m_Q)) * pow(sqrt((m_Q - 1.0) / M_PI) * (tgamma(1.0 / (m_Q - 1.0)) / tgamma((1.0 / (m_Q - 1.0)) - 1.0 / 2.0)), ((2.0 - 2.0 * m_Q) / (3.0 - m_Q)));
        //            d = percentD * Math.exp((-1) * (Math.pow(q - 1.0, 2.0)) / 0.4);
    } else if (m_Q == 1.0) {
        d = m_GeneralizedDiffusion;
    }

    return d;
}
} // end namespace itk

#endif
