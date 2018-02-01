#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkGlobalLogisticContrastEnhancementImageFilter.h"

int main(int argc, char* argv[])
{
    if ( argc < 3 )
    {
        std::cerr << "Missing parameters. " << std::endl;
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0]
                << " inputImageFileName outputImageFileName flipContrastAdjustment(0/1)"
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

    typedef itk::GlobalLogisticContrastEnhancementImageFilter<InputImageType, OutputImageType>  FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(reader->GetOutput());
    if (atof(argv[3])==1) {
        filter->FlipObjectAreaOn();
    }
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
