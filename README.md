dfmm
====

Directional Fast Multipole Method (dFMM)

The dFMM is the implementation of a numerical method which reduces the
quadratic O(N^2) cost the matrix-vector product for a given class of
matrices. By using it, matrices, whose generating kernel function is
asymptotically smooth (respectively, oscillatory), can be applied to vectors
with linear O(N) (respectively, almost linear O(N logN)) cost.

The term directional stems from the fact that oscillatory kernel (eg., the
Helmholtz kernel) functions become smooth if multiplied by a plane wave. This
smoothness is the prerequisite of our method: we use a Chebyshev interpolation
scheme to interpolate the kernel function what paves the road for the
construction of a fast method. Thus, we do not modify the kernel and that is
why our approach is to a large extend kernel-independent.

The C++ implementation comprises both, the bbFMM (black-box FMM) for
asymptotically smooth kernels and the dFMM (directional FMM) for oscillatory
kernels. Moreover, eight different variants for the M2L operator are available
(SA, SArcmp, NA, NAsym, NAblk, IA, IAsym, IAblk). For a study of their
individual weaknesses and strenghts we refer to a forthcoming paper titled
"Optimized M2L Kernels for the Chebyshev Interpolation based Fast Multipole
Method". Moreover, there are three different low-rank approximation schemes to
choose from (truncated SVD, fully and partially pivoted ACA with an optional
subsequent truncated SVD).

The dFMM implementation is under the "BSD 2-Clause" license. Feel free to use
and distribute the code. If you have remarks on bugs, improvements, etc. I
would be very happy if you let me know. Moreover, if you would like to commit
changes let me know too, I will give you the necessary rights.
