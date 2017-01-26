#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkDiffusionEdgeOptimizationImageCalculator.h"

int main(int argc, char* argv[])
{
    if ( argc < 3 )
    {
        std::cerr << "Missing parameters. " << std::endl;
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0]
                << " inputImageFileName optimizationFunction (1,2,3)"
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

    typedef itk::DiffusionEdgeOptimizationImageCalculator<InputImageType>  CalculatorType;
    CalculatorType::Pointer calculator = CalculatorType::New();
    calculator->SetImage(reader->GetOutput());
    int optFunction = atoi(argv[2]);
    switch (optFunction) {
    case 1:
        calculator->SetOptimizationMethod(CalculatorType::CANNY);
        break;
    case 2:
        calculator->SetOptimizationMethod(CalculatorType::MAD);
        break;
    case 3:
        calculator->SetOptimizationMethod(CalculatorType::MORPHOLOGICAL);
        break;
    default:
        break;
    }
    calculator->Compute();

    std::cout<<"Kappa: "<<calculator->GetKappa()<<std::endl;

    return EXIT_SUCCESS;
}
