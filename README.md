# Comprehensive Image Processing & Filtering Toolkit

A complete image processing suite implementing various filtering techniques, morphological operations, and image enhancement methods using Python (Tkinter) and MATLAB.

## Academic Context

This project was developed as part of a comprehensive Computer Vision course practical work (TP). The implementation includes detailed theoretical foundations and experimental results documented in the accompanying academic report.


## Features Overview

### Python Application - Morphological Operations (Interactive GUI)
- **Morphological Operations**: Erosion, Dilation, Opening, Closing
- **Advanced Operations**: Top Hat (Bright/Dark), Gradient, Median Filter
- **Contour Detection**: Internal, External, and Morphological Contours
- **Interactive Controls**: Real-time threshold adjustment
- **Image Management**: Load, save, and reset functionality
- **Dual View**: Side-by-side comparison of original and processed images

### MATLAB Implementation - Complete Filtering Suite
A comprehensive collection of image processing filters and enhancement techniques:

#### **Spatial Domain Filters**
- **Low-Pass Filters**: 
  - `FiltragePassePas_Callback` - Basic low-pass filtering
  - `Moyenneur3_3_Callback` - 3x3 averaging filter
  - `Moyenneur5_5_Callback` - 5x5 averaging filter
  - `Gaussien3_3_Callback` - 3x3 Gaussian filter
  - `Gaussien5_5_Callback` - 5x5 Gaussian filter
  - `Conique_Callback` - Conical filter
  - `Pyramidal_Callback` - Pyramidal filter
  - `TriangleDePascal_Callback` - Pascal's triangle filter

- **High-Pass Filters**:
  - `FiltragePasseHaut` - Basic high-pass filtering
  - `Laplacien` - Laplacian edge detection
  - `Gradient` - Gradient-based edge detection

#### **Frequency Domain Filters**
- **FFT Operations**: `FFT` - Fast Fourier Transform processing
- **Frequency Low-Pass**:
  - `FiltrePasseBasFrequentiel` - Frequency domain low-pass
  - `FiltrePasseBasIdeal_Callback` - Ideal low-pass filter
  - `FiltrePasseBasDeButterwoth_Callback` - Butterworth low-pass filter
- **Frequency High-Pass**:
  - `FiltrePasseHautFrequentiel` - Frequency domain high-pass
  - `FiltrePasseHautIdeal_Callback` - Ideal high-pass filter
  - `FiltrePasseHautDeButterworth_Callback` - Butterworth high-pass filter
- **Band Filters**:
  - `FiltrePasseBandeIdeal_Callback` - Ideal band-pass filter
  - `FiltrePasseBandedeButterworth_Callback` - Butterworth band-pass filter
  - `FiltreARejectiondeBande_Callback` - Band-reject filter
- **Specialized Filters**:
  - `FiltreSpectraleLocale_Callback` - Local spectral filtering
  - `FiltreHomomorphique_Callback` - Homomorphic filtering

#### **Edge Detection Algorithms**
- **Classical Methods**:
  - `Prewitt` - Prewitt edge detector
  - `Roberts_Callback` - Roberts cross-gradient
  - `Sobel_Callback` - Sobel edge detection
  - `Kirsh` - Kirsh edge detection
- **Advanced Methods**:
  - `Marr_Hildreth_Callback` - Marr-Hildreth edge detection
  - `Canny` - Canny edge detection

#### **Image Enhancement & Processing**
- **Brightness & Contrast**:
  - `Luminosite_Callback` - Brightness adjustment
  - `Contraste_Callback` - Contrast enhancement
  - `Negatif_Callback` - Image negative
- **Color Processing**:
  - `NiveauDeGris` - Grayscale conversion
  - `Histogramme_Callback` - Histogram analysis
  - `Binarisation_Callback` - Image binarization
- **Filtering Types**:
  - `Lineaire_Callback` - Linear filtering
  - `NonLineaire_Callback` - Non-linear filtering

#### **Noise Reduction**
- **Noise Types & Removal**:
  - `PoivreEtSel` - Salt and pepper noise handling
  - `Gaussien` - Gaussian noise processing
  - `Median_Callback` - Median filtering for noise reduction

## Project Structure
```
comprehensive-image-processing/
├── README.md                   # This documentation
├── requirements.txt            # Python dependencies
├── .gitignore                 # Git ignore rules
├── src/
│   └── python/
│       └── Bajjou_Abdelilah_Morphological_Operators.py
├── matlab/
│   ├── image_filters.fig      # MATLAB GUI interface
│   ├── image_filters.m        # Complete filtering implementation
│   └── functions/             # Individual filter functions
├── screenshots/
│   ├── morphological/         # Python app screenshots
│   └── matlab_interface/      # MATLAB GUI screenshots (when available)
└── docs/
    ├── filter_documentation.md
    └── usage_examples.md
```

## Installation & Setup

### Python Environment
git clone https://github.com/abdelilahbajjou/comprehensive-image-processing.git
cd comprehensive-image-processing

# Install Python dependencies
pip install -r requirements.txt

# Run the morphological operations tool
python src/python/Bajjou_Abdelilah_Morphological_Operators.py
```

### MATLAB Environment
```matlab
% Open MATLAB and navigate to the matlab/ directory
cd matlab/

% Run the main GUI
run('image_filters.fig')
% OR
image_filters
```

## Usage Examples

### Python Morphological Operations
1. **Load Image**: Import your target image
2. **Select Operation**: Choose from 11 morphological operations
3. **Adjust Parameters**: Use threshold slider for fine-tuning
4. **Chain Operations**: Apply multiple operations sequentially
5. **Save Results**: Export processed images

### MATLAB Comprehensive Filtering
The MATLAB interface provides access to 40+ different filtering and processing techniques:

- **Basic Filtering**: Start with low-pass or high-pass filters
- **Noise Reduction**: Apply median, Gaussian, or averaging filters
- **Edge Detection**: Use Sobel, Canny, or Marr-Hildreth detectors
- **Frequency Analysis**: Perform FFT and frequency domain filtering
- **Enhancement**: Adjust brightness, contrast, and apply histogram operations

## Documentation

- **[Complete Academic Report](report/academic_report.pdf)** - Detailed theoretical analysis and experimental results

## Technical Implementation

### Morphological Operations (Python)
- **Kernel**: 5x5 structuring element
- **Backend**: OpenCV for optimized performance
- **Interface**: Tkinter for cross-platform compatibility

### Filtering Suite (MATLAB)
- **Spatial Domain**: Direct convolution and correlation methods
- **Frequency Domain**: FFT-based filtering with various transfer functions
- **Edge Detection**: Gradient-based and zero-crossing methods
- **Enhancement**: Histogram manipulation and point operations

## Filter Categories

| Category | Count | Examples |
|----------|--------|----------|
| **Spatial Low-Pass** | 8 | Gaussian, Averaging, Conical |
| **Spatial High-Pass** | 3 | Laplacian, Gradient, High-pass |
| **Frequency Domain** | 9 | Ideal, Butterworth, Band-pass |
| **Edge Detection** | 6 | Sobel, Canny, Prewitt, Roberts |
| **Enhancement** | 6 | Contrast, Brightness, Histogram |
| **Noise Handling** | 4 | Median, Gaussian, Salt & Pepper |
| **Morphological** | 11 | Erosion, Dilation, Opening, Closing |

## Applications

This toolkit is suitable for:
- **Computer Vision Research**: Comprehensive preprocessing pipeline
- **Medical Image Analysis**: Noise reduction and enhancement
- **Digital Photography**: Artistic filters and enhancement
- **Industrial Inspection**: Edge detection and morphological analysis
- **Educational Purposes**: Learning image processing concepts

## Future Enhancements

- [ ] Port MATLAB filters to Python for unified interface
- [ ] Add real-time video processing capabilities
- [ ] Implement machine learning-based filters
- [ ] Create batch processing functionality
- [ ] Add performance benchmarking tools

## Requirements

### Python Dependencies
- `opencv-python >= 4.8.0`
- `numpy >= 1.24.0`
- `Pillow >= 10.0.0`
- `tkinter` (usually included with Python)

### MATLAB Requirements
- MATLAB R2013a or later
- Image Processing Toolbox
- Signal Processing Toolbox (for frequency domain operations)

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-filter`)
3. Commit your changes (`git commit -m 'Add amazing new filter'`)
4. Push to the branch (`git push origin feature/amazing-filter`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

**Abdelilah Bajjou**
- GitHub: [@abdelilahbajjou](https://github.com/abdelilahbajjou)
- Email: abdelilahbajjou@gmail.com
- Field: Computer Vision & Image Processing
- Institution: La Faculté des Sciences Dhar El Mahraz –FSDM-
- Course: Computer Vision - Practical Work (TP)

## Acknowledgments

- Computer Vision course instructor and materials
- Academic supervisor for guidance and feedback
- OpenCV community for excellent image processing libraries
- MATLAB Image Processing Toolbox documentation
- Peer reviewers and collaborators

## Academic Contributions

This project demonstrates:
- **Theoretical Understanding**: Comprehensive grasp of image processing fundamentals
- **Practical Implementation**: Working code in multiple programming environments
- **Experimental Analysis**: Detailed results and performance evaluation
- **Documentation**: Professional academic reporting standards
- **Reproducibility**: Complete source code and documentation for replication

## Project Statistics

- **Total Functions**: 40+ filtering and processing functions
- **Implementation Languages**: Python + MATLAB
- **Filter Categories**: 7 major categories
- **GUI Interfaces**: 2 (Python Tkinter + MATLAB GUI)
- **Image Formats Supported**: PNG, JPEG, BMP, TIFF

---

*A comprehensive toolkit for image processing enthusiasts, researchers, and students*
