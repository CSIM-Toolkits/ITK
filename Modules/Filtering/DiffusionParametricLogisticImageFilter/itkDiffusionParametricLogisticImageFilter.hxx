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
#ifndef __itkDiffusionParametricLogisticImageFilter_hxx
#define __itkDiffusionParametricLogisticImageFilter_hxx
#include "itkDiffusionParametricLogisticImageFilter.h"

#include <itkNrrdImageIO.h>
#include <itkImageRegionIterator.h>
#include <itkDiffusionTensor3DReconstructionImageFilter.h>
#include <itkMetaDataObject.h>

#include "itkImageFileWriter.h"

#include <iostream>
using namespace std;

namespace itk
{
template< typename TInput, typename TOutput, typename TMask >
DiffusionParametricLogisticImageFilter< TInput, TOutput, TMask >
::DiffusionParametricLogisticImageFilter()
{
    m_Kappa = 3.0;
    m_MaximumProportion=0.999;
    this->SetNumberOfRequiredInputs(1);
}

template< typename TInput, typename TOutput, typename TMask >
void
DiffusionParametricLogisticImageFilter<TInput, TOutput, TMask>::SetInputImage(const TInput* image)
{
    this->SetNthInput(0, const_cast<TInput*>(image));
}

template< typename TInput, typename TOutput, typename TMask>
void
DiffusionParametricLogisticImageFilter<TInput, TOutput, TMask>::SetDiffusionSpace(const TMask* mask)
{
    this->SetNthInput(1, const_cast<TMask*>(mask));
}

template< typename TInput, typename TOutput, typename TMask >
typename TInput::ConstPointer
DiffusionParametricLogisticImageFilter<TInput, TOutput, TMask>::GetDWIImage()
{
    return static_cast< const TInput * >
            ( this->ProcessObject::GetInput(0) );
}

template< typename TInput, typename TOutput, typename TMask >
typename TMask::ConstPointer
DiffusionParametricLogisticImageFilter<TInput, TOutput, TMask>::GetDiffusionSpace()
{
    return static_cast< const TMask * >
            ( this->ProcessObject::GetInput(1) );
}


template< typename TInput, typename TOutput, typename TMask >
void
DiffusionParametricLogisticImageFilter< TInput, TOutput, TMask >
::GenerateOutputInformation(void){
    OutputImageType *output = this->GetOutput();
    output->CopyInformation(this->GetDWIImage());
    output->SetRegions(this->GetDWIImage()->GetRequestedRegion());

    output->Allocate();
}

template< typename TInput, typename TOutput, typename TMask >
void
DiffusionParametricLogisticImageFilter< TInput, TOutput, TMask >
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
        //TODO Fazer 2*Otsu na imagem B0 para estimar o CSF.
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

    typedef itk::DiffusionTensor3DReconstructionImageFilter<unsigned short, unsigned short, double > TensorReconstructionImageFilterType;
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

    //Step 1: take the mean S0 image if there are more than one non-diffusion volume in the image
    //Creating the new input image with only one b0 volume (the mean of the N non-diffusion images)
    typename InputImageType::Pointer diffusionAcquisitionImage = InputImageType::New();
    diffusionAcquisitionImage->CopyInformation(input);
    diffusionAcquisitionImage->SetRegions(input->GetRequestedRegion());
    diffusionAcquisitionImage->SetVectorLength( numberOfGradientImages + 1 );
    diffusionAcquisitionImage->Allocate();

    createDiffusionSpace(diffusionAcquisitionImage, input, gradients);

    //Step 2: Calculates the diffusion weights for each direction (D = (ln(Si) - ln(S0))/b )
    typename InputImageType::Pointer diffusionImage = InputImageType::New();
    diffusionImage->CopyInformation(input);
    diffusionImage->SetRegions(input->GetRequestedRegion());
    diffusionImage->SetVectorLength( numberOfGradientImages );
    diffusionImage->Allocate();

    typename InputImageType::Pointer filteredDiffusions = InputImageType::New();
    filteredDiffusions->CopyInformation(input);
    filteredDiffusions->SetRegions(input->GetRequestedRegion());
    filteredDiffusions->SetVectorLength( numberOfGradientImages );
    filteredDiffusions->Allocate();

    typename OutputImageType::Pointer nonDiffImage = OutputImageType::New();
    nonDiffImage->CopyInformation(input);
    nonDiffImage->SetRegions(input->GetRequestedRegion());
    nonDiffImage->SetVectorLength( 1 );
    nonDiffImage->Allocate();

    createDiffusionWeightedValues(diffusionAcquisitionImage, diffusionImage, nonDiffImage, numberOfGradientImages, b0);

    //    Step 3: Estimate the image SNR
    double SNR=0.0;
    estimateImageSignalToNoiseRatio(diffusionImage, usedMaskSpace, SNR);

    if (m_DebugMode) {
        std::cout << "Image SNR="<< SNR << std::endl;
    }

    attenuateDWINoise(filteredDiffusions, diffusionImage, usedMaskSpace, SNR, numberOfGradientImages);

    createFilteredDiffusionAcquisition(output, filteredDiffusions, nonDiffImage, b0);
}

template< typename TInput, typename TOutput, typename TMask >
void
DiffusionParametricLogisticImageFilter< TInput, TOutput, TMask >
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
    typedef itk::VariableLengthVector<unsigned short> DiffusionVectorType;

    rawDiffIt.GoToBegin();
    inputIterator.GoToBegin();
    while (!inputIterator.IsAtEnd()) {
        DiffusionVectorType       weights;
        weights.SetSize(numberOfGradientImages + 1);

        //Capturing the mean non-diffusion values
        unsigned short meanNonDiff=static_cast<unsigned short>(0);
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
DiffusionParametricLogisticImageFilter< TInput, TOutput, TMask >
::estimateImageSignalToNoiseRatio(typename InputImageType::Pointer diffImg, typename MaskImageType::Pointer mask, double& SNR)
{
    typedef itk::ImageRegionIterator<InputImageType>            InputRegionIteratorType;
    typedef itk::ImageRegionConstIterator<MaskImageType>        InputMaskRegionConstIteratorType;
    InputRegionIteratorType             diffIt(diffImg, diffImg->GetRequestedRegion());
    InputMaskRegionConstIteratorType    maskIt(mask, mask->GetRequestedRegion());

    double mean=0.0, std=0.0, N=0.0;

    diffIt.GoToBegin();
    maskIt.GoToBegin();
    while (!diffIt.IsAtEnd()) {
        if (maskIt.Get()>static_cast<MaskPixelType>(0)) {
            //Calculating the mean value of the mask region
            for (unsigned int i = 0; i < diffImg->GetNumberOfComponentsPerPixel(); ++i) {
                if (diffIt.Get()[i] > static_cast<double>(0) && !isinf(diffIt.Get()[i])) {
                    mean+=diffIt.Get()[i];
                    N++;
                }
            }
        }
        ++diffIt;
        ++maskIt;
    }
    mean/=N;

    diffIt.GoToBegin();
    maskIt.GoToBegin();
    while (!diffIt.IsAtEnd()) {
        if (maskIt.Get()>static_cast<MaskPixelType>(0)) {
            //Calculating the standard deviation value of the mask region
            for (unsigned int i = 0; i < diffImg->GetNumberOfComponentsPerPixel(); ++i) {
                if (diffIt.Get()[i] > static_cast<double>(0) && !isinf(diffIt.Get()[i])) {
                    std+=pow(diffIt.Get()[i]-mean,2.0);
                }
            }
        }
        ++diffIt;
        ++maskIt;
    }
    std=sqrt(std/(N-1));

    SNR=mean/std;
}

template< typename TInput, typename TOutput, typename TMask >
void
DiffusionParametricLogisticImageFilter< TInput, TOutput, TMask >
::createDiffusionWeightedValues(typename InputImageType::Pointer diffAcquitions, typename InputImageType::Pointer diffImg, typename OutputImageType::Pointer nonDiffImg, unsigned int numberOfGradientImages, unsigned int b0)
{
    typedef itk::ImageRegionIterator<InputImageType>        InputRegionIteratorType;
    typedef itk::ImageRegionIterator<OutputImageType>        OutputRegionIteratorType;
    InputRegionIteratorType     rawDiffIt(diffAcquitions, diffAcquitions->GetRequestedRegion());
    InputRegionIteratorType     diffIt(diffImg, diffImg->GetRequestedRegion());
    OutputRegionIteratorType    nonDiffIt(nonDiffImg, nonDiffImg->GetRequestedRegion());

    typedef itk::VariableLengthVector<double> DiffusionVectorType;
    typedef itk::VariableLengthVector<unsigned short>  NonDiffusionVectorType;

    diffIt.GoToBegin();
    nonDiffIt.GoToBegin();
    rawDiffIt.GoToBegin();
    while (!rawDiffIt.IsAtEnd()) {
        DiffusionVectorType     diffValues;
        NonDiffusionVectorType  nonDiffValues;
        diffValues.SetSize(numberOfGradientImages);
        nonDiffValues.SetSize(1);

        nonDiffValues[0]=static_cast<unsigned short>(rawDiffIt.Get()[0]);

        unsigned int index=0;
        for (unsigned int n = 1; n < diffAcquitions->GetNumberOfComponentsPerPixel(); ++n) {
            if (rawDiffIt.Get()[0]!=0 || rawDiffIt.Get()[n]!=0) {
                diffValues[index]=std::abs<double>((std::log(rawDiffIt.Get()[n]) - std::log(rawDiffIt.Get()[0]))/static_cast<double>(b0));
            }else{
                diffValues[index]=0;
            }
            ++index;
        }
        diffIt.Set(diffValues);
        nonDiffIt.Set(nonDiffValues);

        ++diffIt;
        ++rawDiffIt;
        ++nonDiffIt;
    }
}

template< typename TInput, typename TOutput, typename TMask >
void
DiffusionParametricLogisticImageFilter< TInput, TOutput, TMask >
::attenuateDWINoise(typename InputImageType::Pointer filteredDiffusions, typename InputImageType::Pointer diffImg, typename MaskImageType::Pointer mask, double SNR, unsigned int numberOfGradients)
{
    typedef itk::ImageRegionIterator<InputImageType>        InputRegionIteratorType;
    typedef itk::ImageRegionIterator<MaskImageType>         MaskRegionIteratorType;
    InputRegionIteratorType     diffIt(diffImg, diffImg->GetRequestedRegion());
    InputRegionIteratorType    outputIterator(filteredDiffusions, filteredDiffusions->GetRequestedRegion());
    MaskRegionIteratorType      maskIterator(mask, mask->GetRequestedRegion());

    typedef itk::VariableLengthVector<double> DiffusionVectorType;
    diffIt.GoToBegin();
    outputIterator.GoToBegin();
    maskIterator.GoToBegin();

    while (!outputIterator.IsAtEnd()) {
        if (maskIterator.Get()>static_cast<MaskPixelType>(0)) {
            //Find the mean value of the i-eth voxel
            double mean_i=0.0;
            unsigned int N=0.0;
            for (int j = 0; j < diffImg->GetNumberOfComponentsPerPixel(); ++j) {
                if ( diffIt.Get()[j]>static_cast<double>(0) && !isinf(diffIt.Get()[j]) ) {
                    mean_i+=diffIt.Get()[j];
                    N++;
                }
            }
            mean_i/=N;

            //Attenuate noise using the i-eth mean value and the image SNR
            DiffusionVectorType diffVoxelValues;
            diffVoxelValues.SetSize(numberOfGradients + 1);
            double sigma_i = mean_i/SNR;
             double alpha = log((1.0-m_MaximumProportion)/m_MaximumProportion)/(m_Kappa*sigma_i);
             double beta = m_Kappa*sigma_i;
            for (int j = 0; j < diffImg->GetNumberOfComponentsPerPixel(); ++j) {
                double dev=diffIt.Get()[j]-mean_i;
                double p = sigmoid(abs(dev),alpha,beta);
                diffVoxelValues[j]=diffIt.Get()[j] - p*dev;
            }
            outputIterator.Set(diffVoxelValues);
        }

        ++diffIt;
        ++outputIterator;
        ++maskIterator;
    }
}

template< typename TInput, typename TOutput, typename TMask >
void
DiffusionParametricLogisticImageFilter< TInput, TOutput, TMask >
::createFilteredDiffusionAcquisition(typename OutputImageType::Pointer output, typename InputImageType::Pointer diffImg, typename OutputImageType::Pointer nonDiffImg, unsigned int b0)
{
    typedef itk::ImageRegionIterator<InputImageType>        InputRegionIteratorType;
    typedef itk::ImageRegionIterator<OutputImageType>       OutputRegionIteratorType;
    InputRegionIteratorType     diffIt(diffImg, diffImg->GetRequestedRegion());
    OutputRegionIteratorType    nonDiffIt(nonDiffImg, nonDiffImg->GetRequestedRegion());
    OutputRegionIteratorType    outputIterator(output, output->GetRequestedRegion());

    typedef itk::VariableLengthVector<unsigned short> DiffusionVectorType;
    unsigned int numberOfImages = diffImg->GetNumberOfComponentsPerPixel() + 1;
    diffIt.GoToBegin();
    nonDiffIt.GoToBegin();
    outputIterator.GoToBegin();

    while (!outputIterator.IsAtEnd()) {
        DiffusionVectorType diffValues;
        diffValues.SetSize(numberOfImages);
        for (int i = 0; i < numberOfImages; ++i) {
            if (i==0) {
                diffValues[i] = static_cast<unsigned short>(nonDiffIt.Get()[0]);
            }else{
                diffValues[i] = static_cast<unsigned short>(static_cast<double>(nonDiffIt.Get()[0])*exp(-static_cast<double>(b0)*diffIt.Get()[i-1]));
            }
        }
        outputIterator.Set(diffValues);

        ++outputIterator;
        ++diffIt;
        ++nonDiffIt;
    }
}


template< typename TInput, typename TOutput, typename TMask >
double
DiffusionParametricLogisticImageFilter< TInput, TOutput, TMask >
::sigmoid(double x, double alpha, double beta) {
    return 1/(1+std::exp(-alpha*(x-beta)));
}



} // end namespace itk

#endif
