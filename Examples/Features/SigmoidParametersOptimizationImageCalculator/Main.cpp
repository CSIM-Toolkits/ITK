#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkSigmoidParametersOptimizationImageCalculator.h"

int main(int argc, char* argv[])
{
 if ( argc < 2 )
        {
          std::cerr << "Missing parameters. " << std::endl;
          std::cerr << "Usage: " << std::endl;
          std::cerr << argv[0]
                    << " inputImageFile"
                    << std::endl;
          return -1;
        }

    const unsigned int Dimension = 3;

    typedef float                       PixelType;
    typedef itk::Image<PixelType, Dimension>    ImageType;

    typedef itk::ImageFileReader<ImageType> ReaderType;

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

    typedef itk::SigmoidParametersOptimizationImageCalculator<ImageType> FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(reader->GetOutput());
    filter->SetMaximumAlpha(1.0);
    filter->SetMinimumAlpha(0.0);
    filter->SetMaximumBeta(1.0);
    filter->SetMinimumBeta(0.0);

    filter->SetNumberOfBins(256);
    filter->Update();

    std::cout<<"Alpha: "<<filter->GetOptimumAlpha()<<std::endl;
    std::cout<<"Beta: "<<filter->GetOptimumBeta()<<std::endl;

  return EXIT_SUCCESS;
}
