#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkIsotropicAnomalousDiffusionImageFilter.h"

int main(int argc, char* argv[])
{
 if ( argc < 5 )
        {
          std::cerr << "Missing parameters. " << std::endl;
          std::cerr << "Usage: " << std::endl;
          std::cerr << argv[0]
                    << " inputImageFile outputImageFile GeneralizedDiffusionCoefficient Qvalue NumberOfIteration TimeStep"
                    << std::endl;
          return -1;
        }

    const unsigned int Dimension = 3;

    typedef float                       PixelType;
    typedef itk::Image<PixelType, Dimension>    ImageType;

    typedef itk::ImageFileReader<ImageType> ReaderType;
    typedef itk::ImageFileWriter<ImageType> WriterType;

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

    typedef itk::IsotropicAnomalousDiffusionImageFilter<ImageType, ImageType> FilterType;
    FilterType::Pointer filter = FilterType::New();

    filter->SetInput(reader->GetOutput());
    filter->SetGeneralizedDiffusion(std::atof(argv[3]));
    filter->SetQ(std::atof(argv[4]));
    filter->SetIterations(std::atoi(argv[5]));
    filter->SetTimeStep(std::atof(argv[6]));
    filter->Update();

      WriterType::Pointer writer = WriterType::New();
      writer->SetFileName(argv[2]);
      writer->SetInput( filter->GetOutput() );


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

