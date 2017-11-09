#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkLogisticContrastEnhancementImageFilter.h"

int main(int argc, char* argv[])
{
    if ( argc < 4 )
    {
        std::cerr << "Missing parameters. " << std::endl;
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0]
                << " inputImageFileName outputImageFileName [flipObject:(1 is true or 0 is false)] thrMethod"
                << std::endl;
        return -1;
    }

    const unsigned int Dimension = 3;

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

    typedef itk::LogisticContrastEnhancementImageFilter<InputImageType, OutputImageType>  FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(reader->GetOutput());
    filter->SetNumberOfBins(128);
    if (atof(argv[3])==1) {
        filter->FlipObjectAreaOn();
    }
    filter->SetThresholdMethod(FilterType::ThresholdMethod(atoi(argv[4])));
    filter->Update();

    std::cout<<"Alpha: "<<filter->GetAlpha()<<std::endl;
    std::cout<<"Beta: "<<filter->GetBeta()<<std::endl;

    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(argv[2]);
    writer->SetInput(filter->GetOutput());
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
