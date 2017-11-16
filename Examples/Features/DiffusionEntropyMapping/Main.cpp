#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkDiffusionEntropyMappingImageFilter.h"

int main(int argc, char* argv[])
{
    if ( argc < 4 )
    {
        std::cerr << "Missing parameters. " << std::endl;
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0]
                << " inputImageFileName outputImageFileName qValue[Default=1.0] HistogramsBinsType[Default=2]"
                << std::endl;
        return -1;
    }

    const unsigned int InputDimension = 3;
    const unsigned int OutputDimension = 3;

    typedef double                       PixelType;
    typedef itk::VectorImage<PixelType, InputDimension>    InputImageType;
    typedef itk::Image<PixelType, OutputDimension>    OutputImageType;
    double q = (argc>=4?atof(argv[3]):1.0);
    unsigned char HistBins = (argc>=5?atof(argv[4]):2);

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

    cout<<"reader output region: "<<reader->GetOutput()->GetRequestedRegion()<<" Number of Components: "<<reader->GetOutput()->GetNumberOfComponentsPerPixel()<<endl;

    typedef itk::DiffusionEntropyMappingImageFilter<InputImageType>  FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(reader->GetOutput());

    filter->SetHistogramBins(HistBins);
    filter->SetQValue(q);
    filter->DebugModeOn();
    filter->Update();

    cout<<"filter output region: "<<filter->GetOutput()->GetRequestedRegion()<<"Filter Number of Components: "<<filter->GetOutput()->GetNumberOfComponentsPerPixel()<<endl;
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
