#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkDiffusionSelfInformationMappingImageFilter.h"

int main(int argc, char* argv[])
{
    if ( argc < 4 )
    {
        std::cerr << "Missing parameters. " << std::endl;
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0]
                << " inputImageFileName diffusionMaskFileName outputImageFileName qValue HistogramsBinsType"
                << std::endl;
        return -1;
    }

    const unsigned int InputDimension = 3;
    const unsigned int OutputDimension = 3;

    typedef float                       PixelType;
    typedef itk::VectorImage<PixelType, InputDimension>    InputImageType;
    typedef itk::Image<unsigned char, InputDimension>    InputMaskType;
    typedef itk::Image<PixelType, OutputDimension>    OutputImageType;
    double q = 0.0;
    if (argc>4) {
        q = atof(argv[4]);
    }else{
        q = atof(argv[3]);
    }
    unsigned char HistBins = atof(argv[5]);

    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typedef itk::ImageFileReader<InputMaskType> ReaderMaskType;
    typedef itk::ImageFileWriter<OutputImageType> WriterType;
    typedef itk::ImageFileWriter<InputImageType> VectorWriterType;

    ReaderType::Pointer reader = ReaderType::New();
    ReaderMaskType::Pointer readerMask = ReaderMaskType::New();
    reader->SetFileName(argv[1]);

    if (argc>4) {
        readerMask->SetFileName(argv[2]);
        readerMask->Update();
    }

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

    cout<<"reader output region: "<<reader->GetOutput()->GetRequestedRegion()<<" Number of Components: "<<reader->GetOutput()->GetNumberOfComponentsPerPixel()<<endl;

    typedef itk::DiffusionSelfInformationMappingImageFilter<InputImageType>  FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(reader->GetOutput());
    if (argc>4) {
        filter->SetDiffusionSpace(readerMask->GetOutput());
    }

    filter->SetQValue(q);
    filter->DebugModeOn();
    filter->Update();

    cout<<"Filter output region: "<<filter->GetOutput()->GetRequestedRegion()<<"Filter Number of Components: "<<filter->GetOutput()->GetNumberOfComponentsPerPixel()<<endl;
    typename WriterType::Pointer writer = WriterType::New();
    if (argc>4) {
        writer->SetFileName(argv[3]);
    }else{
        writer->SetFileName(argv[2]);
    }

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
