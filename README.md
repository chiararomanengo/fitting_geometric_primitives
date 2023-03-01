# fitting_geometric_primitives

This repository contains a Matlab library for fitting simple geometric primitives on point clouds. As simple primitives we mean the classical surface primitives derived from constructive solid geometry, i.e., planes, spheres, cylinders, cones and tori. 

The ```code``` directory contains our method based on the Hough Transform (see [1]). This method has been compared with state of the art algorithms, and the results can be found in [2].

Instead, the ```test``` directory contains 100 point clouds that can be used to perform tests and evaluate performance. These have been extracted from [SHREC22](https://github.com/chiararomanengo/SHREC2022) dataset.


## Basic Usage
To use this method, you can simply run the ```main.m``` file in Matlab. The input is a .txt file containing the point cloud. The point cloud is then processed and the method find the best fitting primitive type and its geometric descriptors.  

If the point cloud is noisy, you can run the ```mainNoise.m```. This is configured to deal with perturbed point clouds, i.e., thresholds are properly tuned w.r.t. the ```main.m``` file.

# References

[1] A. Raffo, C. Romanengo, B. Falcidieno, S. Biasotti, "Fitting and recognition of geometric primitives in segmented 3D point clouds using a localized voting procedure", Computer Aided Geometric Design, Vol. 97, 2022, pp. 102123, https://doi.org/10.1016/j.cagd.2022.102123.

[2] C. Romanengo, A. Raffo, S. Biasotti, B. Falcidieno, V. Fotis, I. Romanelis, E. Psatha, K. Moustakas, I. Sipiran, Q.-T. Nguyen, C.-B. Chu,K.-N. Nguyen-Ngoc, D.-K. Vo, T.-A. To, N.-T. Nguyen, N.-Q. Le-Pham, H.-D. Nguyen, M.-T. Tran, Y. Qie, N. Anwer, "SHREC 2022: Fitting and recognition of simple geometric primitives on point clouds", Volume 107, 2022, Pages 32-49, https://doi.org/10.1016/j.cag.2022.07.004.
