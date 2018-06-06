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
#ifndef __itkDiffusionSelfInformationMappingImageFilter_hxx
#define __itkDiffusionSelfInformationMappingImageFilter_hxx
#include "itkDiffusionSelfInformationMappingImageFilter.h"

#include <itkNrrdImageIO.h>
#include <itkImageRegionIterator.h>
#include <itkDiffusionTensor3DReconstructionImageFilter.h>
#include <itkMetaDataObject.h>

#include <iostream>
using namespace std;

namespace itk
{
template< typename TInput, typename TOutput, typename TMask >
DiffusionSelfInformationMappingImageFilter< TInput, TOutput, TMask >
::DiffusionSelfInformationMappingImageFilter()
{
    m_QValue=1.0;
    m_HistogramBins=2;
    m_UseManualNumberOfBins=false;
    this->SetNumberOfRequiredInputs(1);
}

template< typename TInput, typename TOutput, typename TMask >
void
DiffusionSelfInformationMappingImageFilter<TInput, TOutput, TMask>::SetInputImage(const TInput* image)
{
    this->SetNthInput(0, const_cast<TInput*>(image));
}

template< typename TInput, typename TOutput, typename TMask>
void
DiffusionSelfInformationMappingImageFilter<TInput, TOutput, TMask>::SetDiffusionSpace(const TMask* mask)
{
    this->SetNthInput(1, const_cast<TMask*>(mask));
}

template< typename TInput, typename TOutput, typename TMask >
typename TInput::ConstPointer
DiffusionSelfInformationMappingImageFilter<TInput, TOutput, TMask>::GetDWIImage()
{
    return static_cast< const TInput * >
            ( this->ProcessObject::GetInput(0) );
}

template< typename TInput, typename TOutput, typename TMask >
typename TMask::ConstPointer
DiffusionSelfInformationMappingImageFilter<TInput, TOutput, TMask>::GetDiffusionSpace()
{
    return static_cast< const TMask * >
            ( this->ProcessObject::GetInput(1) );
}


template< typename TInput, typename TOutput, typename TMask >
void
DiffusionSelfInformationMappingImageFilter< TInput, TOutput, TMask >
::GenerateOutputInformation(void){
    OutputImageType *output = this->GetOutput();
    output->CopyInformation(this->GetDWIImage());
    output->SetRegions(this->GetDWIImage()->GetRequestedRegion());
    output->SetNumberOfComponentsPerPixel( 1 );

    output->Allocate();
    output->FillBuffer(0);
}

template< typename TInput, typename TOutput, typename TMask >
void
DiffusionSelfInformationMappingImageFilter< TInput, TOutput, TMask >
::GenerateData()
{
    unsigned int numberOfImages = 0;
    unsigned int numberOfGradientImages = 0;
    bool readb0 = false;
    double b0 = 0;
    std::vector<unsigned int> gradients;
    typename InputImageType::ConstPointer input = this->GetDWIImage();
    typename OutputImageType::Pointer output = this->GetOutput();
    typename MaskImageType::Pointer usedMaskSpace;
    usedMaskSpace = MaskImageType::New();
    usedMaskSpace->CopyInformation(output);
    usedMaskSpace->SetRegions(output->GetRequestedRegion());
    usedMaskSpace->Allocate();

    if (this->GetNumberOfIndexedInputs()==1) {
        //If no mask is provided, the entire image is processed.
        usedMaskSpace->FillBuffer(1);
    }else{
        typedef itk::ImageRegionIterator<MaskImageType>             InputRegionIteratorType;
        typedef itk::ImageRegionConstIterator<MaskImageType>        InputRegionConstIteratorType;
        InputRegionConstIteratorType   inputMaskIt(this->GetDiffusionSpace(), this->GetDiffusionSpace()->GetRequestedRegion());
        InputRegionIteratorType        usedMaskIt(usedMaskSpace, usedMaskSpace->GetRequestedRegion());

        inputMaskIt.GoToBegin();
        usedMaskIt.GoToBegin();
        while (!inputMaskIt.IsAtEnd()) {
            usedMaskIt.Set(inputMaskIt.Get());
            ++inputMaskIt;
            ++usedMaskIt;
        }
    }

    typedef itk::DiffusionTensor3DReconstructionImageFilter<OutputPixelType, OutputPixelType, double > TensorReconstructionImageFilterType;
    // -------------------------------------------------------------------------
    // Parse the Nrrd headers to get the B value and the gradient directions used
    // for diffusion weighting.
    //
    // The Nrrd headers should look like :
    // The tags specify the B value and the gradient directions. If gradient
    // directions are (0,0,0), it indicates that it is a reference image.
    //
    // DWMRI_b-value:=800
    // DWMRI_gradient_0000:= 0 0 0
    // DWMRI_gradient_0001:=-1.000000       0.000000        0.000000
    // DWMRI_gradient_0002:=-0.166000       0.986000        0.000000
    // DWMRI_gradient_0003:=0.110000        0.664000        0.740000
    // ...
    //
    itk::MetaDataDictionary imgMetaDictionary = input->GetMetaDataDictionary();
    std::vector<std::string> imgMetaKeys = imgMetaDictionary.GetKeys();
    std::vector<std::string>::const_iterator itKey = imgMetaKeys.begin();
    std::string metaString;
    typename TensorReconstructionImageFilterType::GradientDirectionType vect3d;
    typename TensorReconstructionImageFilterType::GradientDirectionContainerType::Pointer DiffusionVectors = TensorReconstructionImageFilterType::GradientDirectionContainerType::New();

    for (; itKey != imgMetaKeys.end(); ++itKey)
    {
        double x,y,z;
        itk::ExposeMetaData<std::string> (imgMetaDictionary, *itKey, metaString);
        if (itKey->find("DWMRI_gradient") != std::string::npos)
        {
            if (m_DebugMode) {
                std::cout << *itKey << " ---> " << metaString << std::endl;
            }
            sscanf(metaString.c_str(), "%lf %lf %lf\n", &x, &y, &z);
            vect3d[0] = x; vect3d[1] = y; vect3d[2] = z;
            DiffusionVectors->InsertElement( numberOfImages, vect3d );
            ++numberOfImages;
            // If the direction is 0.0, this is a reference image
            if (vect3d[0] == 0.0 &&
                    vect3d[1] == 0.0 &&
                    vect3d[2] == 0.0)
            {
                gradients.push_back(0);
                continue;
            }
            //            Add a new gradient volume
            gradients.push_back(1);
            ++numberOfGradientImages;
        }
        else if (itKey->find("DWMRI_b-value") != std::string::npos)
        {
            if (m_DebugMode) {
                std::cout << *itKey << " ---> " << metaString << std::endl;
            }
            readb0 = true;
            b0 = atof(metaString.c_str());
        }
    }
    if (m_DebugMode) {
        std::cout << "Number of gradient images: "
                  << numberOfGradientImages
                  << " and Number of reference images: "
                  << numberOfImages - numberOfGradientImages
                  << std::endl;
    }
    if(!readb0)
    {
        std::cerr << "BValue not specified in header file" << std::endl;
    }

    //First: take the mean S0 image if there are more than one non-diffusion volume in the image
    //    numberOfReferenceImages = numberOfImages - numberOfGradientImages;

    //First: take the mean S0 image if there are more than one non-diffusion volume in the image
    //Creating the new input image with only one b0 volume (the mean of the N non-diffusion images)
    typename InputImageType::Pointer diffusionAcquisitionImage = InputImageType::New();
    diffusionAcquisitionImage->CopyInformation(input);
    diffusionAcquisitionImage->SetRegions(input->GetRequestedRegion());
    diffusionAcquisitionImage->SetVectorLength( numberOfGradientImages + 1 );
    diffusionAcquisitionImage->Allocate();

    createDiffusionSpace(diffusionAcquisitionImage, input, gradients);

    //Secondly: Calculates the diffusion weights for each direction (D = (ln(Si) - ln(S0))/b )
    typename InputImageType::Pointer diffusionImage = InputImageType::New();
    diffusionImage->CopyInformation(input);
    diffusionImage->SetRegions(input->GetRequestedRegion());
    diffusionImage->SetVectorLength( numberOfGradientImages );
    diffusionImage->Allocate();

    createDiffusionWeightedValues(diffusionAcquisitionImage, diffusionImage, numberOfGradientImages, b0);

    //Setting the number of bins depending on a bin function.
    if (m_UseManualNumberOfBins) {
        if (m_DebugMode) {
            cout<<"Manual number of bins: "<<m_HistogramBins<<endl;
        }
    }else{
        m_HistogramBins = automaticHistogramBinCalculation(numberOfGradientImages);
        if (m_DebugMode) {
            cout<<"Automatic number of bins: "<<m_HistogramBins<<endl;
        }
    }

    //Finding input minimum and maximum values only on the gradient volumes and also construct the global probability distribution
    OutputPixelType minimumInputValue = NumericTraits<OutputPixelType>::max();
    OutputPixelType maximumInputValue = NumericTraits<OutputPixelType>::min();
    getSpaceMaximumMinimumDiffusion(diffusionImage, usedMaskSpace, maximumInputValue, minimumInputValue);

    if (m_DebugMode) {
        std::cout<<"Diffusion space range (max: "<<maximumInputValue<<" - min: "<<minimumInputValue<<")"<<std::endl;
    }

    typename HistogramType::Pointer prioryProbabilityDistribution = HistogramType::New();
    typename HistogramType::SizeType size(1);

    size.Fill(m_HistogramBins);

    typename HistogramType::MeasurementVectorType lowerBound;
    lowerBound.SetSize(m_HistogramBins);
    lowerBound.Fill(minimumInputValue);

    typename HistogramType::MeasurementVectorType upperBound;
    upperBound.SetSize(m_HistogramBins);
    upperBound.Fill(maximumInputValue);
    prioryProbabilityDistribution->SetMeasurementVectorSize(1);
    prioryProbabilityDistribution->Initialize(size, lowerBound, upperBound );

    typename HistogramType::IndexType index(1);
    typename HistogramType::MeasurementVectorType mv(1);

    createPriorProbabilityDistribution(diffusionImage, prioryProbabilityDistribution, index, mv);

    if (m_DebugMode) {
        for (int var = 0; var < prioryProbabilityDistribution->GetSize()[0]; ++var) {
            std::cout<<"bin["<<var<<"]: "<<prioryProbabilityDistribution->GetFrequency(var)<<std::endl;
        }
        std::cout<<"Total frequencies: "<<prioryProbabilityDistribution->GetTotalFrequency()<<std::endl;
    }

    calculatesEntropyMapping(output, diffusionImage, usedMaskSpace, prioryProbabilityDistribution, index, mv);
}

template< typename TInput, typename TOutput, typename TMask >
unsigned int
DiffusionSelfInformationMappingImageFilter< TInput, TOutput, TMask >
::automaticHistogramBinCalculation(unsigned int n)
{
    switch (m_HistogramBins) {
    case SQUAREROOT:
        return static_cast<unsigned int>(sqrt(n));
        break;
    case STURGES:
        return static_cast<unsigned int>(log2(n)+1.0);
        break;
    case RICE:
        return static_cast<unsigned int>(2.0*pow(n,(1.0/3.0)));
        break;
    default:
        break;
    }
}

template< typename TInput, typename TOutput, typename TMask >
void
DiffusionSelfInformationMappingImageFilter< TInput, TOutput, TMask >
::createDiffusionSpace(typename InputImageType::Pointer diffImg, typename InputImageType::ConstPointer inputImg, std::vector<unsigned int> gradientsList)
{
    unsigned int numberOfGradientImages=0, numberOfReferenceImages=0;
    typedef itk::ImageRegionConstIterator<InputImageType>   InputRegionConstIteratorType;
    typedef itk::ImageRegionIterator<InputImageType>        InputRegionIteratorType;
    InputRegionConstIteratorType   inputIterator(inputImg,inputImg->GetRequestedRegion());
    InputRegionIteratorType    rawDiffIt(diffImg, diffImg->GetRequestedRegion());

    //Detect the number of gradients images
    for (int i = 0; i < gradientsList.size(); ++i) {
        if (gradientsList[i]!=0) {
            numberOfGradientImages++;
        }
    }

    numberOfReferenceImages=inputImg->GetNumberOfComponentsPerPixel() - numberOfGradientImages;

    typedef itk::VariableLengthVector<OutputPixelType> DiffusionVectorType;

    rawDiffIt.GoToBegin();
    inputIterator.GoToBegin();
    while (!inputIterator.IsAtEnd()) {
        DiffusionVectorType       weights;
        weights.SetSize(numberOfGradientImages + 1);

        //Capturing the mean non-diffusion values
        OutputPixelType meanNonDiff=static_cast<OutputPixelType>(0);
        for (unsigned int n = 0; n < inputImg->GetNumberOfComponentsPerPixel(); ++n) {
            if (gradientsList[n]==0){
                meanNonDiff+=inputIterator.Get()[n];
            }
        }
        meanNonDiff/=numberOfReferenceImages;
        weights[0]=meanNonDiff;

        unsigned int index=1;
        for (unsigned int n = 0; n < inputImg->GetNumberOfComponentsPerPixel(); ++n) {
            if (gradientsList[n]!=0){
                weights[index]=inputIterator.Get()[n];
                ++index;
            }
        }
        rawDiffIt.Set(weights);

        ++rawDiffIt;
        ++inputIterator;
    }
}

template< typename TInput, typename TOutput, typename TMask >
void
DiffusionSelfInformationMappingImageFilter< TInput, TOutput, TMask >
::getSpaceMaximumMinimumDiffusion(typename InputImageType::Pointer diffImg, typename MaskImageType::Pointer mask, OutputPixelType& maximum, OutputPixelType& minimum)
{
    typedef itk::ImageRegionIterator<InputImageType>        InputRegionIteratorType;
    typedef itk::ImageRegionConstIterator<MaskImageType>        InputMaskRegionConstIteratorType;
    InputRegionIteratorType    diffIt(diffImg, diffImg->GetRequestedRegion());
    InputMaskRegionConstIteratorType    maskIt(mask, mask->GetRequestedRegion());

    diffIt.GoToBegin();
    maskIt.GoToBegin();
    while (!diffIt.IsAtEnd()) {
        if (maskIt.Get()!=static_cast<MaskPixelType>(0)) {
            for (int min = 0; min < diffImg->GetNumberOfComponentsPerPixel(); ++min) {
                if (diffIt.Get()[min]<minimum && diffIt.Get()[min]>0){
                    minimum=diffIt.Get()[min];
                }
            }

            for (int max = 0; max < diffImg->GetNumberOfComponentsPerPixel(); ++max) {
                if (diffIt.Get()[max]>maximum && !isinf(diffIt.Get()[max])){
                    maximum=diffIt.Get()[max];
                }
            }
        }
        ++diffIt;
        ++maskIt;
    }
}

template< typename TInput, typename TOutput, typename TMask >
void
DiffusionSelfInformationMappingImageFilter< TInput, TOutput, TMask >
::createDiffusionWeightedValues(typename InputImageType::Pointer diffAcquitions, typename InputImageType::Pointer diffImg, unsigned int numberOfGradientImages, unsigned int b0)
{
    typedef itk::ImageRegionIterator<InputImageType>        InputRegionIteratorType;
    InputRegionIteratorType    rawDiffIt(diffAcquitions, diffAcquitions->GetRequestedRegion());
    InputRegionIteratorType    diffIt(diffImg, diffImg->GetRequestedRegion());

    typedef itk::VariableLengthVector<OutputPixelType> DiffusionVectorType;

    diffIt.GoToBegin();
    rawDiffIt.GoToBegin();
    while (!rawDiffIt.IsAtEnd()) {
        DiffusionVectorType       diffValues;
        diffValues.SetSize(numberOfGradientImages);

        unsigned int index=0;
        for (unsigned int n = 1; n < diffAcquitions->GetNumberOfComponentsPerPixel(); ++n) {
            if (rawDiffIt.Get()[0]!=0 || rawDiffIt.Get()[n]!=0) {
                diffValues[index]=std::abs((std::log(rawDiffIt.Get()[n]) - std::log(rawDiffIt.Get()[0]))/static_cast<OutputPixelType>(b0));
            }else{
                diffValues[index]=0;
            }
            ++index;
        }
        diffIt.Set(diffValues);

        ++diffIt;
        ++rawDiffIt;
    }
}

template< typename TInput, typename TOutput, typename TMask >
void
DiffusionSelfInformationMappingImageFilter< TInput, TOutput, TMask >
::createPriorProbabilityDistribution(typename InputImageType::Pointer diffImg, typename HistogramType::Pointer prioryProbabilityDistribution, typename HistogramType::IndexType index, typename HistogramType::MeasurementVectorType mv)
{
    typedef itk::ImageRegionIterator<InputImageType>        InputRegionIteratorType;
    InputRegionIteratorType    diffIt(diffImg, diffImg->GetRequestedRegion());
    diffIt.GoToBegin();
    while (!diffIt.IsAtEnd()) {
        //Calculate the global probability distribution
        for (int j = 0; j < diffImg->GetNumberOfComponentsPerPixel(); ++j) {
            if ( diffIt.Get()[j]>static_cast<OutputPixelType>(0) && !isinf(diffIt.Get()[j]) ) {
                mv[0]=diffIt.Get()[j];
                prioryProbabilityDistribution->GetIndex(mv,index);
                prioryProbabilityDistribution->IncreaseFrequencyOfIndex(index, 1);
            }
        }

        ++diffIt;
    }
}

template< typename TInput, typename TOutput, typename TMask >
void
DiffusionSelfInformationMappingImageFilter< TInput, TOutput, TMask >
::calculatesEntropyMapping(typename OutputImageType::Pointer output, typename InputImageType::Pointer diffImg, typename MaskImageType::Pointer mask, typename HistogramType::Pointer prioryProbabilityDistribution, typename HistogramType::IndexType index, typename HistogramType::MeasurementVectorType mv)
{
    //    Iterate over the entire diffusion space to get the information contained in each voxel
    typedef itk::ImageRegionConstIterator<InputImageType>        InputRegionIteratorType;
    typedef itk::ImageRegionIterator<OutputImageType>       OutputRegionIteratorType;
    typedef itk::ImageRegionConstIterator<MaskImageType>         MaskRegionIteratorType;
    InputRegionIteratorType         diffIt(diffImg, diffImg->GetRequestedRegion());
    OutputRegionIteratorType        outputIterator(output, output->GetRequestedRegion());
    MaskRegionIteratorType          maskIt(mask, mask->GetRequestedRegion());
    diffIt.GoToBegin();
    maskIt.GoToBegin();
    outputIterator.GoToBegin();
    while (!outputIterator.IsAtEnd()) {
        if (maskIt.Get()>static_cast<MaskPixelType>(0)) {
            //        Extract the vector data from voxel.
            //        This iterative process runs through each voxel of the input image
            //        where the diffusion gradients are allocated. With each voxel values
            //        trends, it is possible to reconstruct a local histogram and then calculate
            //        the local self information using the average log(Di).

            //Calculate entropy
            OutputPixelType seflInfo = 0, p;
            for (int j = 0; j < diffImg->GetNumberOfComponentsPerPixel(); ++j) {
                mv[0]=diffIt.Get()[j];
                prioryProbabilityDistribution->GetIndex(mv,index);
                p = static_cast<double>(prioryProbabilityDistribution->GetFrequency(index))/static_cast<double>(prioryProbabilityDistribution->GetTotalFrequency());
                if (m_QValue != 1.0) {
                    if( p > NumericTraits<double>::min() )
                    {
                        seflInfo += pow<double>(p,m_QValue);
                    }
                }else{
                    if( p > NumericTraits<double>::min() )
                    {
                        seflInfo += - (p)
                                * (std::log( p ) / std::log( 2.0 ));
                    }
                }
            }

            if (m_QValue != 1.0) {
                outputIterator.Set((1.0 - seflInfo)/(m_QValue - 1.0));
            }else{
                outputIterator.Set(seflInfo);
            }
        }else{
            outputIterator.Set(0);
        }

        ++diffIt;
        ++outputIterator;
        ++maskIt;
    }
}



} // end namespace itk

#endif
