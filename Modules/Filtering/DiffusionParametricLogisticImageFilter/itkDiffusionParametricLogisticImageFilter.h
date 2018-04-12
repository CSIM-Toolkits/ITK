/*
   Copyright 2017 Antonio Carlos da Silva Senra Filho and Fabricio Henrique Simozo

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
 */
#ifndef __itkDiffusionParametricLogisticImageFilter_h
#define __itkDiffusionParametricLogisticImageFilter_h
#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkNumericTraits.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkHistogram.h"
#include "itkDenseFrequencyContainer2.h"
#include <math.h>

namespace itk
{

template< typename TInputImage, typename TOutputImage, typename TInputMask=Image<unsigned char, 3> >
class ITK_EXPORT DiffusionParametricLogisticImageFilter:
        public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
    /** Extract dimension from inputs images, where it is assumed there are with the same type. */
    itkStaticConstMacro(InputImageDimension, unsigned int,
                        TInputImage::ImageDimension);
    itkStaticConstMacro(OutputImageDimension, unsigned int,
                        TOutputImage::ImageDimension);

    /** Convenient typedefs for simplifying declarations. */
    typedef TInputImage     InputImageType;
    typedef TOutputImage    OutputImageType;
    typedef TInputMask      MaskImageType;

    /** Standard class typedefs. */
    typedef DiffusionParametricLogisticImageFilter          Self;
    typedef ImageToImageFilter< TInputImage, TOutputImage >       Superclass;
    typedef SmartPointer< Self >                                  Pointer;
    typedef SmartPointer< const Self >                            ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(DiffusionParametricLogisticImageFilter, ImageToImageFilter)

    typedef typename InputImageType::PixelType                  InputPixelType;
    typedef typename OutputImageType::PixelType                 OutputPixelType;
    typedef typename MaskImageType::PixelType                   MaskPixelType;

    void SetInputImage(const TInputImage* image);
    void SetDiffusionSpace(const TInputMask* mask);


    /** Debug mode is used to inform some messages in the standard output. */
    itkBooleanMacro(DebugMode)
    itkSetMacro(DebugMode, bool)

    /** Set the kappa value that regularize the noise standard deviation threshold. */
    itkSetMacro(Kappa, double)

    /** Set the maximum proportion used in the noise attenuation process. The more close to 1.0, the more intense will be the filtering. */
    itkSetMacro(MaximumProportion, double)

    itkGetMacro(DebugMode, bool)
    itkGetMacro(Kappa, double)
    itkGetMacro(MaximumProportion, double)

#ifdef ITK_USE_CONCEPT_CHECKING
    // Begin concept checking
    itkConceptMacro( InputHasNumericTraitsCheck,
                     ( Concept::HasNumericTraits< InputPixelType > ) );
    itkConceptMacro( SameDimensionCheck,
                     ( Concept::SameDimension< InputImageDimension, OutputImageDimension > ) );
#endif

protected:
    DiffusionParametricLogisticImageFilter();
    virtual ~DiffusionParametricLogisticImageFilter() {}
    virtual void GenerateOutputInformation(void) ITK_OVERRIDE;

    typename TInputImage::ConstPointer GetDWIImage();
    typename TInputMask::ConstPointer GetDiffusionSpace();
    void GenerateData();
private:
    DiffusionParametricLogisticImageFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
    void createDiffusionSpace(typename InputImageType::Pointer diffImg, typename InputImageType::ConstPointer inputImg, std::vector<unsigned int> gradientsList);
    void createDiffusionWeightedValues(typename InputImageType::Pointer diffAcquitions, typename InputImageType::Pointer diffImg, typename InputImageType::Pointer nonDiffImg, unsigned int numberOfGradients, unsigned int b0);
    void estimateImageSignalToNoiseRatio(typename InputImageType::Pointer diffImg, typename MaskImageType::Pointer mask, double& SNR);
    void attenuateDWINoise(typename OutputImageType::Pointer output, typename InputImageType::Pointer diffImg, typename MaskImageType::Pointer mask, double SNR, unsigned int numberOfGradients);
    void createFilteredDiffusionAcquisition(typename OutputImageType::Pointer output, typename InputImageType::Pointer diffImg, typename InputImageType::Pointer nonDiffImg, unsigned int b0);
    double sigmoid(double x, double alpha, double beta);
    double m_Kappa, m_MaximumProportion;
    bool m_DebugMode;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiffusionParametricLogisticImageFilter.hxx"
#endif

#endif
