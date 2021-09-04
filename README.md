# Panel_ChangeinMean
This simple package provides the Matlab code to implement the V_{N,T}(1), V_{N,T}(2), H_{N,T}(1) and H_{N,T}(2) tests for detecting mean change points in high dimensional panel dataset, with both large time-series dimension T and cross-section dimension N and no specific relationship between N and T is required.
Reference:
Horváth, L., Liu, Z., Rice, G. and Zhao, Y., (2021). Detecting common breaks in the means of high dimensional cross-dependent panels. The Econometrics Journal.
The package contains two main files, VNT_tests.m and HNT_tests.m

## VNT_tests.m
The file provides codes for computing the test statistic V_{N,T} and the two types of bootstrapping critical values Algorithms (1) and (2) for assessing the null hypothesis of no change in mean.

## HNT_tests.m
The file provides codes for computing the test statistic H_{N,T} and the two types of bootstrapping critical values Algorithms (1) and (2) for assessing the null hypothesis of no change in mean.

## Note
Regarding the tests J_{N,T} and C_{N,T} discussed in the paper, I do not publicly provide their source code as they refer to different references, but they are available upon request.
