/*
 * CXSparse: a Concise Sparse matrix package.
 * Copyright (C) 2006-2011, Timothy A. Davis.
 * Copyright (C) 2011-2012, Richard W. Lincoln.
 * http://www.cise.ufl.edu/research/sparse/CXSparse
 *
 * -------------------------------------------------------------------------
 *
 * CXSparseJ is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * CXSparseJ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this Module; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
 *
 */

package edu.emory.mathcs.csparsej.tdcomplex.test ;

import java.io.InputStream;
import java.util.Random;

import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_add.cs_add;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_chol.cs_chol;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex.cs_cmult;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex.cs_cone;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_ipvec.cs_ipvec;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_lsolve.cs_lsolve;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_ltsolve.cs_ltsolve;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_multiply.cs_multiply;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_permute.cs_permute;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_pinv.cs_pinv;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_pvec.cs_pvec;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_schol.cs_schol;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_transpose.cs_transpose;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_updown.cs_updown;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_util.cs_spalloc;

import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcs;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcsa;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcsn;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcss;

/**
 * Read a matrix, solve a linear system, update/downdate.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class DZcs_test3 extends DZcs_test {

	/**
	 * Cholesky update/downdate
	 *
	 * @param prob
	 *            problem
	 * @return true if successful, false on error
	 */
	protected static boolean test3(DZproblem prob)
	{
		DZcs A, C, W = null, WW, WT, E = null, W2 ;
		int n, k, Li[], Lp[], Wi[], Wp[], p1, p2, p[] = null ;
		boolean ok ;
		DZcsa b, x, resid, y = null, Lx, Wx ;
		double[] s ;
		double t, t1 ;
		DZcss S = null ;
		DZcsn N = null ;
		if (prob == null || prob.sym == 0 || prob.A.n == 0) return (false) ;
		A = prob.A ; C = prob.C ; b = prob.b ; x = prob.x ; resid = prob.resid ;
		n = A.n ;
		if (prob.sym == 0 || n == 0) return (true) ;
		rhs (x, b, n) ;						/* compute right-hand side */
		System.out.print("\nchol then update/downdate ") ;
		print_order (1) ;
		y = new DZcsa(n) ;
		t = tic() ;
		S = cs_schol (1, C) ;					/* symbolic Chol, amd(A+A') */
		System.out.printf("\nsymbolic chol time %8.2f ms\n", toc (t)) ;
		t = tic() ;
		N = cs_chol (C, S) ;					/* numeric Cholesky */
		System.out.printf("numeric  chol time %8.2f ms\n", toc (t)) ;
		if (S == null || N == null) return (false) ;
		t = tic() ;
		cs_ipvec (S.pinv, b, y, n) ;				/* y = P*b */
		cs_lsolve (N.L, y) ;					/* y = L\y */
		cs_ltsolve (N.L, y) ;					/* y = L'\y */
		cs_pvec (S.pinv, y, x, n) ;				/* x = P'*y */
		System.out.printf("solve    chol time %8.2f ms\n", toc (t)) ;
		System.out.printf("original: ") ;
		print_resid (true, C, x, b, resid, prob) ;		/* print residual */
		k = n / 2 ;						/* construct W  */
		W = cs_spalloc (n, 1, n, true, false) ;
		Lp = N.L.p ; Li = N.L.i ; Lx = new DZcsa(N.L.x) ;
		Wp = W.p ; Wi = W.i ; Wx = new DZcsa (W.x) ;
		Wp [0] = 0 ;
		p1 = Lp [k] ;
		Wp [1] = Lp [k+1] - p1 ;
		s = Lx.get(p1) ;
		Random r = new Random(1) ;
		for ( ; p1 < Lp [k+1] ; p1++)
		{
			p2 = p1 - Lp [k] ;
			Wi [p2] = Li [p1] ;
			Wx.set(p2, cs_cmult(s, r.nextDouble())) ;
		}
		t = tic() ;
		ok = cs_updown (N.L, +1, W, S.parent) ;			/* update: L*L'+W*W' */
		t1 = toc (t) ;
		System.out.printf("update:   time: %8.2f ms\n", t1) ;
		if (!ok) return (false) ;
		t = tic() ;
		cs_ipvec (S.pinv, b, y, n) ;				/* y = P*b */
		cs_lsolve (N.L, y) ;					/* y = L\y */
		cs_ltsolve (N.L, y) ;					/* y = L'\y */
		cs_pvec (S.pinv, y, x, n) ;				/* x = P'*y */
		t = toc (t) ;
		p = cs_pinv (S.pinv, n) ;
		W2 = cs_permute (W, p, null, true) ;			/* E = C + (P'W)*(P'W)' */
		WT = cs_transpose (W2, true) ;
		WW = cs_multiply (W2, WT) ;
		WT = null ;
		W2 = null ;
		E = cs_add (C, WW, cs_cone(), cs_cone()) ;
		WW = null ;
		if (E == null || p == null) return (false) ;
		System.out.printf("update:   time: %8.2f ms(incl solve) ", t1 + t) ;
		print_resid (true, E, x, b, resid, prob) ;		/* print residual */
		N = null ;						/* clear N */
		t = tic() ;
		N = cs_chol (E, S) ;					/* numeric Cholesky */
		if (N == null) return (false) ;
		cs_ipvec (S.pinv, b, y, n) ;				/* y = P*b */
		cs_lsolve (N.L, y) ;					/* y = L\y */
		cs_ltsolve (N.L, y) ;					/* y = L'\y */
		cs_pvec (S.pinv, y, x, n) ;				/* x = P'*y */
		t = toc (t) ;
		System.out.printf("rechol:   time: %8.2f ms(incl solve) ", t) ;
		print_resid (true, E, x, b, resid, prob) ;		/* print residual */
		t = tic() ;
		ok = cs_updown (N.L, -1, W, S.parent) ;			/* downdate: L*L'-W*W' */
		t1 = toc (t) ;
		if (!ok) return (false) ;
		System.out.printf("downdate: time: %8.2f\n", t1) ;
		t = tic() ;
		cs_ipvec (S.pinv, b, y, n) ;				/* y = P*b */
		cs_lsolve (N.L, y) ;					/* y = L\y */
		cs_ltsolve (N.L, y) ;					/* y = L'\y */
		cs_pvec (S.pinv, y, x, n) ;				/* x = P'*y */
		t = toc (t) ;
		System.out.printf("downdate: time: %8.2f ms(incl solve) ", t1 + t) ;
		print_resid (true, C, x, b, resid, prob) ;		/* print residual */
		return (true) ;
	}

	public void test_c4()
	{
		InputStream in = get_stream (C4) ;
		DZproblem prob = get_problem (in, 0) ;

		assert_problem(prob, 4, 4, 10, -1, 16, 7.37e+01) ;

		test3(prob) ;

		double x_norm = 8.3862 ;
		assertEquals(x_norm, prob.norms.get(0), DELTA) ;
		x_norm = 9.0064 ;
		assertEquals(x_norm, prob.norms.get(1), 0.3) ;  // TODO: tighten accuracy
		assertEquals(x_norm, prob.norms.get(2), 0.3) ;
		x_norm = 8.3862 ;
		assertEquals(x_norm, prob.norms.get(3), DELTA) ;
	}

	public void test_mhd1280b()
	{
		InputStream in = get_stream (MHD1280B) ;
		DZproblem prob = get_problem (in, 0) ;

		assert_problem(prob, 1280, 1280, 12029, -1, 22778, 79.9740) ;

		test3(prob) ;

		double x_norm = 76270143066.4197 ;
		assertEquals(x_norm, prob.norms.get(0), DELTA) ;
		assertEquals(x_norm, prob.norms.get(1), DELTA) ;
		assertEquals(x_norm, prob.norms.get(2), DELTA) ;
		assertEquals(x_norm, prob.norms.get(3), DELTA) ;
	}

}
