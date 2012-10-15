/*! 

@mainpage

<a name="dfmm"></a>
<h2>What is dFMM?</h2> 

dFMM stands for [d]irectional [F]ast [M]ultipole [M]ethod and is a C++
implementation of the the method presented in ''Fast directional multilevel
summation for oscillatory kernels based on chebyshev interpolation''
in <em>Journal of Computational Physics, 231(4):1175 – 1196, 2012</em> by
M. Messner, M. Schanz and E. Darve. The library is heavily based on template
programming techniques in order to achieve extensibility. The main goal is to
offer a flexible set of independent modules that allow a straight forward use
in existing libraries.

For completeness we insert the sections of the abstract of the above mentioned
paper:<br />

<em>Many applications lead to large systems of linear equations with dense
matrices.  Direct matrix-vector products become prohibitive, since the
computational cost increases quadratically with the size of the problem. ...By
exploiting specific kernel properties fast algorithms can be
constructed. ...We present a directional multilevel algorithm based on the
%Chebyshev interpolation for translation-invariant oscillatory kernels of the
type \f$K(x,y) = G(x-y) \, \mathsf{e}^{\imath k |x-y|}\f$, with \f$G(x-y)\f$
being any smooth kernel. We consider the modification of the kernel \f$K^u =
K(x,y) \, e^{- \imath k u \cdot (x-y)}\f$, where \f$u\f$ is a unit vector. We
proved that the kernel \f$K^u\f$ can be interpolated efficiently when
\f$x-y\f$ lies in a cone of direction \f$u\f$.
</em>

Indeed the mathematical scheme of the presented FMM is a generalization of
''The black-box fast multipole method'' in <em>Journal of Computational
Physics, 228(23):8712–8725, 2009</em> by W. Fong and E. Darve to oscillatory
kernel functions. The main references for the directional ideas are
''Multilevel computations of integral transforms and particle interactions
with oscillatory kernels'' in <em>Computer Physics Communications,
65(1-3):24–38, 1991</em> by A. Brandt and ''Fast directional multilevel
algorithms for oscillatory kernels'' in <em>SIAM Journal on Scientific
Computing, 29(4):1710–1737, 2007</em> by B. Engquist and L. Ying.


<a name="offers"></a>
<h2>What offers dFMM?</h2> 

dFMM offers a generic FMM implementation which computes the matrix-vector
product in \f$\mathcal{O}(N)\f$ (class of matrices based non-oscillatory
kernel functions), respectively, in \f$\mathcal{O}(N\log N)\f$ (class of
matrices based on oscillatory kernel functions) operations. Moreover, it
offers
<ol>
	<li> Support for two and three space dimensions (given at compile-time as
	DIM (template parameter)).

	<li> Different accuracies (given at compile time as interpolation order
	\f$\ell\f$ (ORDER - template parameter) and at run-time by the accuracy of
	the low-rank approximation of the M2L kernels a \f$\varepsilon\f$ (eps))
		
	<li> <strong>TO CONTINUE ...</strong>

</ol>


<a name="modules"></a>
<h2>How does dFMM work?</h2>

This is the main page for class and function documentation in dFMM. Many of
the classes in dFMM can be grouped into modules (see
the <a href="modules.html">Modules page</a> or the corresponding entry in the
menu at the top of this page).

<strong>TO CONTINUE ...</strong>

<a name="opensource"></a>
<h2>Is dFMM free software?</h2>

Yes, dFMM is free software. It is provided under the <a
href="http://opensource.org/licenses/BSD-2-Clause">BSD 2-Clause
License</a>. If you would like to help adding new features or removing bugs in
the code feel free to contact <a href="mailto:burgerbua@gmail.com">us</a>. We
appreciate any suggestion and recommendation on dFMM and every bug news since
there will certainly be loads of them.

*/
