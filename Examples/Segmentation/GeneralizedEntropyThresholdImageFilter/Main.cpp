#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkGeneralizedEntropyThresholdImageFilter.h"

using namespace std;

int main(int argc, char* argv[])
{
    if ( argc < 4 )
    {
        std::cerr << "Missing parameters. " << std::endl;
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0]
                << " inputImageFileName outputImageFileName qValue "
                << std::endl;
        return -1;
    }

    const unsigned int Dimension = 2;

    typedef float                       PixelType;
    typedef itk::Image<PixelType, Dimension>    InputImageType;
    typedef itk::Image<PixelType, Dimension>    OutputImageType;

    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typedef itk::ImageFileWriter<OutputImageType> WriterType;

    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(argv[1]);
    try
    {
        reader->Update();
    }
    catch ( itk::ExceptionObject &err)
    {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return -1;
    }

    typedef itk::GeneralizedEntropyThresholdImageFilter<InputImageType, OutputImageType> MaxEntropyType;
    typename MaxEntropyType::Pointer maxEn = MaxEntropyType::New();
    maxEn->SetInput(reader->GetOutput());

    maxEn->SetQ(atof(argv[3]));
    maxEn->SetOutsideValue(1);
    maxEn->SetInsideValue(0);
    maxEn->SetNumberOfHistogramBins(256);
    maxEn->Update();

    cout<<"Threshold: "<<maxEn->GetThreshold()<<endl;
    cout<<"Q: "<<maxEn->GetQ()<<endl;

    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(argv[2]);
    writer->SetInput(maxEn->GetOutput());
    try
    {
        writer->Update();
    }
    catch ( itk::ExceptionObject &err)
    {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return -1;
    }

    return EXIT_SUCCESS;
}
