
# SliceMatchingTransport

This repository is the companion repository for our work **Measure transport via stochastic slicing and matching** by Shiying Li and Caroline Moosmueller. 

## Introduction
This project includes sample Matlab codes for morphing a source image to a target image using stochastic slice-and-matching transport.

- **ImageMorph_SliceMatching.m**: Tests the iterative scheme in our study for morphing a source image to a target image. It produces plots of relative error from the morphed images at different iterations compared to the target image.

- **Supporting Functions**: 
  - **slicetransport_matrix.m**: Performs one-step *matrix-slice transport* from a source image to a target image using two orthogonal angles (resulting in an orthogonal matrix).

  - **slicetransport_theta.m**: Executes one-step *single-slice transport*, involving only one angle.
  - **plot_error.m**: Plot the relative errors measured in sliced-Wasserstein distances
  - **SW.m**: Computes the sliced-Wasserstein between two images. 
