# Rotation_angle_estimation_harmonic  
We propose aggregating a candidate peak and its harmonics as a more robust metric in determining the rotation-dependent peak. This work will be presented in the **Session : Multimedia Forensics, 6/8/2023 8:15:00 (EEST), Greece**, IEEE ICASSP 2023. You are very welcomed to discuss with us face to face. 

![image](https://github.com/zengh5/Rotation_angle_estimation_harmonic/blob/main/Figs/Aggregating_harmonics.png)

Compared with the state of the arts, the accuracy of image rotation angle estimation can be improved consistently.

![image](https://github.com/zengh5/Rotation_angle_estimation_harmonic/blob/main/Figs/Overall_comparison.png)

W. Wei, S. Wang, X. Zhang and Z. Tang, “Estimation of image rotation angle using interpolation-related spectral signatures with application to blind detection of image forgery,” IEEE Trans. Inf. Forensics Secur., 5(3): 507–517, 2010.  
D. Vázquez-Padín, C. Mosquera, and F. Pérez-González, “Two-dimensional statistical test for the presence of almost cyclostationarity on images,” in Int. Conf. Image Processing, 2010, pp. 1745–1748.  
C. Chen, J. Ni and Z. Shen, “Effective estimation of image rotation angle using spectral method,” IEEE Signal Processing Letters, 21(7): 890–894, 2014. (https://github.com/ChenglongChen/image-rotation-angle-estimation)  
S. Dai, Y. Zhang, W. Song, et al. “Rotation angle estimation of JPEG compressed image by cyclic spectrum analysis,” Electronics, 2019, 8(12): 1431.

In the supplementary file 'supp.pdf', we provide more detailed results:

- Rotation angle estimation accuracy comparison on uncompresed, JPEG compressed, and scaled-then-rotated images;
- Estimation accuracy for images with different sizes; 
- Ablation study of the parameter N_har. 

## Usage
Please run Onesample_harmonics.m to see the difference between the proposed method and previous ones.  
