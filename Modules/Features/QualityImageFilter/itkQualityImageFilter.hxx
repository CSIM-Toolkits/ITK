#ifndef __itkQualityImageFilter_hxx
#define __itkQualityImageFilter_hxx
#include "itkQualityImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkStatisticsImageFilter.h>
#include <math.h>

namespace itk
{
template< typename TImage >
QualityImageFilter< TImage >
::QualityImageFilter()
{
    this->SetNumberOfRequiredInputs(2);
}

template< typename TImage >
void
QualityImageFilter< TImage >
::SetReferenceImage(const TImage *reference){
    this->SetNthInput(0, const_cast<TImage*>(reference));
}

template< typename TImage >
void
QualityImageFilter< TImage >
::SetCompareImage(const TImage *compared){
    this->SetNthInput(1, const_cast<TImage*>(compared));
}


template< typename TImage >
double
QualityImageFilter< TImage >
::SNR(){
    typedef itk::ImageRegionConstIterator< TImage>          ConstInputIteratorType;
    typename TImage::ConstPointer referenceImage = static_cast<const TImage*>(this->ProcessObject::GetInput(0));
    typename TImage::ConstPointer compareImage = static_cast<const TImage*>(this->ProcessObject::GetInput(1));

    ConstInputIteratorType referenceImageIt(referenceImage, referenceImage->GetRequestedRegion());
    ConstInputIteratorType compareImageIt(compareImage, compareImage->GetRequestedRegion());

    double answer_SNR = 0.0d;

    int N = 1;
    for (int dimN = 0; dimN < TImageDimension; ++dimN) {
        N *= referenceImageIt.GetRegion().GetSize()[dimN];
    }

    double refPix = 0.0d;
    double inPix = 0.0d;
    referenceImageIt.GoToBegin();
    compareImageIt.GoToBegin();

    while (!referenceImageIt.IsAtEnd()) {
        inPix = static_cast<double>(referenceImageIt.Get());
        refPix = static_cast<double>(compareImageIt.Get());

        if(refPix!=inPix && refPix!=0)
            answer_SNR += (10.0 * log10(pow(refPix, 2)/(pow(refPix - inPix, 2))))/static_cast<double>(N);

        ++referenceImageIt;
        ++compareImageIt;
    }

    return answer_SNR;
}

template< typename TImage >
double
QualityImageFilter< TImage >
::RMSE(){
    typedef itk::ImageRegionConstIterator< TImage>          ConstInputIteratorType;
    typename TImage::ConstPointer referenceImage = static_cast<const TImage*>(this->ProcessObject::GetInput(0));
    typename TImage::ConstPointer compareImage = static_cast<const TImage*>(this->ProcessObject::GetInput(1));

    ConstInputIteratorType referenceImageIt(referenceImage, referenceImage->GetRequestedRegion());
    ConstInputIteratorType compareImageIt(compareImage, compareImage->GetRequestedRegion());

    double answer_RMSE = 0.0d;

    int N = 1;
    for (int dimN = 0; dimN < TImageDimension; ++dimN) {
        N *= referenceImageIt.GetRegion().GetSize()[dimN];
    }

    double refPix = 0.0d;
    double inPix = 0.0d;
    referenceImageIt.GoToBegin();
    compareImageIt.GoToBegin();

    while (!referenceImageIt.IsAtEnd()) {
        inPix = static_cast<double>(referenceImageIt.Get());
        refPix = static_cast<double>(compareImageIt.Get());

        if(refPix!=inPix && refPix!=0)
            answer_RMSE += pow(refPix - inPix, 2);

        ++referenceImageIt;
        ++compareImageIt;
    }

    return sqrt(answer_RMSE/static_cast<double>(N));
}

template< typename TImage >
double
QualityImageFilter< TImage >
::MAE(){
    typedef itk::ImageRegionConstIterator< TImage>          ConstInputIteratorType;
    typename TImage::ConstPointer referenceImage = static_cast<const TImage*>(this->ProcessObject::GetInput(0));
    typename TImage::ConstPointer compareImage = static_cast<const TImage*>(this->ProcessObject::GetInput(1));

    ConstInputIteratorType referenceImageIt(referenceImage, referenceImage->GetRequestedRegion());
    ConstInputIteratorType compareImageIt(compareImage, compareImage->GetRequestedRegion());

    double answer_MAE = 0.0d;

    int N = 1;
    for (int dimN = 0; dimN < TImageDimension; ++dimN) {
        N *= referenceImageIt.GetRegion().GetSize()[dimN];
    }

    double refPix = 0.0d;
    double inPix = 0.0d;
    referenceImageIt.GoToBegin();
    compareImageIt.GoToBegin();

    while (!referenceImageIt.IsAtEnd()) {
        inPix = static_cast<double>(referenceImageIt.Get());
        refPix = static_cast<double>(compareImageIt.Get());

        if(refPix!=inPix && refPix!=0)
            answer_MAE += static_cast<double>(abs(refPix - inPix))/static_cast<double>(N);

        ++referenceImageIt;
        ++compareImageIt;
    }

    return answer_MAE;
}

template< typename TImage >
double
QualityImageFilter< TImage >
::SSIM()
{
    typedef itk::ImageRegionConstIterator< TImage>          ConstInputIteratorType;
    typename TImage::ConstPointer img1 = static_cast<const TImage*>(this->ProcessObject::GetInput(0));
    typename TImage::ConstPointer img2 = static_cast<const TImage*>(this->ProcessObject::GetInput(1));

    ConstInputIteratorType img1It(img1, img1->GetRequestedRegion());
    ConstInputIteratorType img2It(img2, img2->GetRequestedRegion());

    typedef itk::StatisticsImageFilter<TImage> StatisticsImageFilterType;
    typename  StatisticsImageFilterType::Pointer statInput1ImageFilter = StatisticsImageFilterType::New ();
    typename  StatisticsImageFilterType::Pointer statInput2ImageFilter = StatisticsImageFilterType::New ();
    statInput1ImageFilter->SetInput(img1);
    statInput2ImageFilter->SetInput(img2);
    statInput1ImageFilter->Update();
    statInput2ImageFilter->Update();

    double answer_SSIM = 0.0d;

    int N = 1;
    for (int dimN = 0; dimN < TImageDimension; ++dimN) {
        N *= img1It.GetRegion().GetSize()[dimN];
    }

    double mu1, mu2, mu1Pow, mu2Pow;
    double sigma1, sigma2, sigma1Pow, sigma2Pow;
    double K1 = 0.01;
    double K2 = 0.03;
    double C1 = pow(K1 * 255, 2);
    double C2 = pow(K2 * 255, 2);

    mu1 = static_cast<double>(statInput1ImageFilter->GetMean());
    mu2 = static_cast<double>(statInput2ImageFilter->GetMean());

    mu1Pow = pow(mu1, 2);
    mu2Pow = pow(mu2, 2);

    sigma1 = static_cast<double>(statInput1ImageFilter->GetSigma());
    sigma2 = static_cast<double>(statInput2ImageFilter->GetSigma());

    sigma1Pow = pow(sigma1, 2);
    sigma2Pow = pow(sigma2, 2);

    answer_SSIM = ((2.0 * mu1 * mu2 + C1) * (2.0 * stdCorrelation() + C2)) / ((mu1Pow + mu2Pow + C1) * (sigma1Pow + sigma2Pow + C2));

    return answer_SSIM;
}

template< typename TImage >
void
QualityImageFilter< TImage >
::GenerateData()
{

}

template<typename TImage>
double
QualityImageFilter<TImage>
::stdCorrelation() {
    typename TImage::ConstPointer img1 = static_cast<const TImage*>(this->ProcessObject::GetInput(0));
    typename TImage::ConstPointer img2 = static_cast<const TImage*>(this->ProcessObject::GetInput(1));
    typedef itk::ImageRegionConstIterator< TImage>          ConstInputIteratorType;
    ConstInputIteratorType img1It(img1, img1->GetRequestedRegion());
    ConstInputIteratorType img2It(img2, img2->GetRequestedRegion());
    double mu1, mu2;
    double sigma = 0;
    int M = 1;
    for (int dimN = 0; dimN < TImageDimension; ++dimN) {
        M *= img1It.GetRegion().GetSize()[dimN];
    }
    typedef itk::StatisticsImageFilter<TImage> StatisticsImageFilterType;
    typename  StatisticsImageFilterType::Pointer stat1 = StatisticsImageFilterType::New ();
    typename  StatisticsImageFilterType::Pointer stat2 = StatisticsImageFilterType::New ();
    stat1->SetInput(img1);
    stat2->SetInput(img2);
    stat1->Update();
    stat2->Update();

    mu1 = static_cast<double>(stat1->GetMean());
    mu2 = static_cast<double>(stat2->GetMean());

    img1It.GoToBegin();
    img2It.GoToBegin();

    while (!img1It.IsAtEnd()) {
        sigma += (static_cast<double>(img1It.Get()) - mu1)*(static_cast<double>(img2It.Get()) - mu2);
        ++img1It;
        ++img2It;
    }

    return (sigma / M);
}

} // end namespace itk

#endif
