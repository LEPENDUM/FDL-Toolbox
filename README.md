
# FDL toolbox

## Description

This is the FDL toolbox providing matlab functions for Light Fields processing in the Fourier Disparity Layer representation (see the [project page](https://v-sense.scss.tcd.ie/research/a-fourier-disparity-layer-representation-for-light-fields/)).

The main features are listed here (detailed documentations of the referenced functions are available in the matlab (.m) files):

### 1. Fast Fourier transform (FFT) and inverse transform tools (see ./utils/FFT_tools)

- **fftImgs** and **ifftImgs** (both cpu and gpu):
	- Functions for 2D FFT and inverse FFT for a stack of images.
	- Options for signal extension (padding and windowing of padded borders).
	- Supports hexagonal sampling of the input images (e.g. RAW Lytro data).
	
- **refocusFFT** (both cpu and gpu):
	- Fourier domain implementation of the "Shift and sum" light field refocusing algorithm.
	
### 2. FDL Calibration algorithms (see ./Calibration):

- **CalibrateFDL_UVD_cpu** and **CalibrateFDL_UVD_gpu**:
	- Gradient descent based calibration to determine view positions on the camera plane (U, V) and disparity values (D) of the layers (see algorithm description in the TIP 2019 paper [1]).
	- Support for independent U,V,D parameters per color channel (e.g. for chromatic aberration estimation).
	- Options for estimating only D knowing U,V, or only U,V knowing D.
	- Options for improved robustness (e.g. see 'UseSignFlip', 'numIterLowFreqs' options in the function documentation).
	
- **CalibrateFDL_NoR1_cpu** and **CalibrateFDL_NoR1_gpu**:
	- Calibration with a relaxed version of the FDL model with no Rank-1 constraint on the parameter matrix containing the shift of each layer to reconstruct each view (see details in [1]).
	- FDL constructed from this relaxed model are only suitable to reconstruct images at the same view positions as the input views (e.g. for denoising or spatial super-resolution).
	- It is recommended to initialise the matrices of parameters using the calibration results of the rank 1 constrained calibration (CalibrateFDL_UVD_xxx functions).


### 3. FDL construction algorithms (see ./FDL_Construction/):

- **ComputeFDL_cpu** and **ComputeFDL_gpu**: 
	- Simple FDL construction from the TIP 2019 paper [1].
	- Use l2 and 2nd order view regularisation.
	
- **ComputeFDL_SparseReg_cpu** (cpu only): 
	- Other version of the FDL contruction. 
	- Use l1 and l2 regularisation.
	
- **ComputeFDL_SuperRes_gpu** (gpu only):
	-  Super-resolution algorithm from the ICCP 2020 paper [2].
	- Supports hexagonal input sampling (e.g. RAW Lytro data).
	- Supports color regularisation to reduce color noise/artifacts (see details in [2]).
	- Supports spatial and angular preconditioning (see details in [2] for spatial pre-conditioning).
	- This method can also be used for FDL construction without super-resolution (super-resolution factor = 1).
	
- **FDL_Complete_SuperRes** (gpu only):
	- Completion algorithm from the ICCP 2020 paper [2]:
	- Supports all the tools from the super-resolution algorithm (e.g. completion jointly with super-resolution and color regularisation for demosaicing).


### 4. Interactive Light Field Rendering application from a FDL model (see ./Render):

-  **RenderModel** (both cpu and gpu)
	- Class for  fast light field rendering from a FDL model.
	- Control of viewpoint, focus, aperture shape, aperture size.
	- Possibility to save the FDL model (using mat file system).

- **RenderAppMain** (both cpu and gpu)
	- GUI application for FDL rendering based on the RenderModel class.
	- User interation for controlling viewpoint, focus, aperture shape and aperture size.
	- Automatic refocusing on click using a disparity map (if not available, a fast disparity estimation can be made).
	- Visualisation of the Fourier magnitude spectrum of the rendered result.
	- Possibility to save the FDL model (using mat file system).

### 6. FDL Tree Codec (see ./FDLTree_Codec):

- Encoder/Decoder for FDL models using a tree structure as described in the paper: [M. Le Pendu, C. Ozcinar and A. Smolic, "Hierarchical Fourier Disparity Layer Transmission for Light Field Streaming," ICIP 2020.](https://v-sense.scss.tcd.ie/wp-content/uploads/2020/05/FDL_Stream_final.pdf)
	- See documentation in ./FDLTree_Codec.README.md.

### 5. Example application scripts (see ./Demo):

- **Demo_FDL_Simple** (both cpu and gpu)
	- Demo script showing an example usage of the FDL processing chain (method without super-resolution from [1]):
	- Loads views, apply Fourier Transform, FDL calibration, FDL construction and launch Rendering GUI.

- **Demo_FDL_SR** (gpu only)
	- Similar to Demo_FDL_Simple but using the FDL construction method with super-resolution from [2].

- **Demo_FDL_View_Extract** (gpu only)
	- Script for Lytro camera view extraction using FDL completion, super-resolution and  color demosaicing in [2]. The [V-Sense Light Field Toolbox](https://github.com/V-Sense/LFToolbox-CLIM_VSENSE) is required to previously convert RAW sensor data into incomplete views along their masks indicating the missing pixels.

- **Demo_FDLTree_Codec** (both cpu and gpu)
	- Demo script of the FDL Tree codec. The demo generates a FDL with the Demo_FDL_Simple script, encodes/decodes it in the FDL Tree representation, and displays decoded results for different levels of the tree.


## Requirements

The toolbox has been developped using Matlab 2017a and Windows 10.

For the GPU versions of the tools, Matlab's Parallel Processing Toolbox and a CUDA GPU are required.

## References

When using the FDL toolbox in your reasearch, please cite the relevant papers:

[\[1\] M. Le Pendu, C. Guillemot and A. Smolic, "A Fourier Disparity Layer Representation for Light Fields," in IEEE Transactions on Image Processing, vol. 28, no. 11, pp. 5740-5753, Nov. 2019.](https://v-sense.scss.tcd.ie/wp-content/uploads/2019/05/FDL_preprint.pdf)

[\[2\] M. Le Pendu, A. Smolic "High Resolution Light Field Recovery with Fourier Disparity Layer Completion, Demosaicing, and Super-Resolution", International Conference on Computational Photography (ICCP) 2020.](https://v-sense.scss.tcd.ie/wp-content/uploads/2017/10/FDLSR_ICCP_preprint.pdf)
