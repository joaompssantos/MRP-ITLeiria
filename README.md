# Minimum Rate Predictors family of encoders

Working implementations of various codecs based on Minimum Rate Predictors (MRP) resulting from my research activities at Instituto de Telecomunicações - Leiria Branch.

Codecs implemented in this code:
- MRP Video [1]
- 4D-MRP [2]
- DT-4D-MRP [2]
- M-MRP [3]
- H-MRP [4]

## Requirements

In the source code of each prediction mode a CMakeLists.txt file is provided. The only requirements are a working CMake and C/C++ environment.

## How to use?

- Each implemented codec has a different git branch.
- After selecting the desired branch create a sub folder to run CMake, example for linux:

~~~
mkdir build
cd build
cmake ../
~~~

- Each binary can present a help text when run without options, such as the following example of H-MRP:

~~~
IT - Leiria: Hierarchical Minimum Rate Predictors
encmrp/decmrp version 1.1.0 (November 2021)
usage: encmrp [options] infile outfile
options:
    -J 2 * num  Views dimensions (in pixels) [H W]*
    -K 2 * num  Dimensions of the array of views [H W]*
    -L 2 * num  Prediction order (Intra, Inter) [-1 1]
    -c str      Hierarchical encoding configuration file*
    -R num      Number of Random Access regions in {0, 2, 4, 5, 9} [0]
    -S num      Number of reference SAIs [4]
    -T num      Maximum reference SAI distance [8]
    -v          Use disparity compensation
    -B 2 * num  Set the block size and the maximum disparity for the disparity compensation [8 8]
    -s          Use the current level causal SAI as reference
    -b num      Bit depth [8]
    -E num      Endianness: little-endian = 0, big-endian = 1. Default: little-endian
    -D num      Distance between views [1]
    -M num      Number of predictors [-1]
    -P num      Precision of prediction coefficients (fractional bits) [6]
    -V num      Number of probability models [16]
    -A num      Accuracy of probability models [3]
    -I num      Maximum number of iterations [100]
    -m          Use MMSE predictors
    -h          Use Huffman coding
    -f          Fixed block-size for adaptive prediction
    -u          Deactivate the histogram packing procedures
    -d          Create extra debug output (coefficients, partitions, etc.)
    -r str      Light field file format [SAI]*. Supported formats:
                    MIA; --> Not yet implemented
                    PVS;
                    SAI.
infile:         Input file (must be in a raw YUV format)
outfile:        Output file

Note: * stands for a mandatory option.
~~~

- The MRP family encoders compress each image component separately.
- Each component should be provided as a raw file, different codecs have different requirements:
    - MRP Video: single image or video sequence (the resolution height and width must be multiple of 8);
    - M-MRP: LF re-arranged as an array of SAIs (the resolution height and width must be multiple of 8);
    - 4D-MRP, DT-4D-MRP, and H-MRP: LF re-arranged either as a pseudo-video sequence of SAIs or as an array of SAIs.
- The configurations used for each codec can be found on their associated publications.
- An example of a hierarchical layers configuration file is provided in `hi-mrp-example.cfg`.
- The following show some example on how to run H-MRP:

~~~
# Encode each component separately, Y
./encmrp.exe -E 0 -K 13 13 -J 434 625 -c h-mrp-sai-13-raster-scalable-hl-4.cfg -S 2 -L 12 13 -D 1 -b 11 Bikes_LE_GRAY_Y.yuv Bikes_LE_GRAY_Y.mrp

# Encode each component separately, U
./encmrp.exe -E 0 -K 13 13 -J 434 625 -c h-mrp-sai-13-raster-scalable-hl-4.cfg -S 2 -L 12 13 -D 1 -b 11 Bikes_LE_GRAY_U.yuv Bikes_LE_GRAY_U.mrp

# Encode each component separately, V
./encmrp.exe -E 0 -K 13 13 -J 434 625 -c h-mrp-sai-13-raster-scalable-hl-4.cfg -S 2 -L 12 13 -D 1 -b 11 Bikes_LE_GRAY_V.yuv Bikes_LE_GRAY_V.mrp
~~~

## References

[1] J. M. Santos, A. F. R. Guarda, L. A. da Silva Cruz, N. M. M. Rodrigues, S. M. M. Faria, **[Compression of medical images using MRP with bi-directional prediction and histogram packing](https://ieeexplore.ieee.org/document/7906386)**, in: Picture Coding Symposium (PCS), Nuremberg, Germany, 2016, pp. 1–5.

[2] J. M. Santos, L. A. Thomaz, P. A. A. Assunção, L. A. da Silva Cruz, L. Távora, S. M. M. Faria, Lossless Coding of Light Fields based on 4D Minimum Rate Predictors, [Online]. Available: https://arxiv.org/abs/2104.06252 (2021). arXiv:2104.06252 (In review process).

[3] J. M. Santos, P. A. A. Assuncao, L. A. da S. Cruz, L. M. N. Tavora, R. Fonseca-Pinto, S. M. M. Faria,  **[Lossless compression of Light Fields using multi-reference Minimum Rate Predictors](https://ieeexplore.ieee.org/document/8712634)**, in: Data Compression Conference (DCC), Snowbird, UT, USA, 2019, pp. 408–417.

[4] Submitted for publishing.
