#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkBrainLogisticSegmentationImageFilter.h"

int main(int argc, char* argv[])
{
    if ( argc < 3 )
    {
        std::cerr << "Missing parameters. " << std::endl;
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0]
                << " inputImageFileName outputImageFileName nTissues"
                << std::endl;
        return -1;
    }

    const unsigned int Dimension = 3;

    typedef float                       PixelType;
    typedef itk::Image<PixelType, Dimension>    InputImageType;
    typedef itk::VectorImage<PixelType, Dimension>    OutputImageType;

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

    typedef itk::BrainLogisticSegmentationImageFilter<InputImageType>  FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(reader->GetOutput());
    filter->SetNumberOfTissues(atoi(argv[3]));
    filter->DebugModeOn();
    filter->Update();

    cout<<"Histogram Peaks: [ ";
    for (int i = 0; i < filter->GetHistogramPeaks().size(); ++i) {
        cout<<filter->GetHistogramPeaks()[i]<<" ";
    }
    cout<<" ]"<<endl;

    cout<<"Histogram Valleys: [ ";
    for (int i = 0; i < filter->GetHistogramValleys().size(); ++i) {
        cout<<filter->GetHistogramValleys()[i]<<" ";
    }
    cout<<" ]"<<endl;

    cout<<"Betas: [ ";
    for (int i = 0; i < filter->GetBetas().size(); ++i) {
        cout<<filter->GetBetas()[i]<<" ";
    }
    cout<<" ]"<<endl;

    cout<<"Alphas: [ ";
    for (int i = 0; i < filter->GetAlphas().size(); ++i) {
        cout<<filter->GetAlphas()[i]<<" ";
    }
    cout<<" ]"<<endl;

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
