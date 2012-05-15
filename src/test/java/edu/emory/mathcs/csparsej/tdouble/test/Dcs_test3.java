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

package edu.emory.mathcs.csparsej.tdouble.test ;

import java.io.InputStream;
import java.util.Random;

import static edu.emory.mathcs.csparsej.tdouble.Dcs_add.cs_add;
import static edu.emory.mathcs.csparsej.tdouble.Dcs_chol.cs_chol;
import static edu.emory.mathcs.csparsej.tdouble.Dcs_ipvec.cs_ipvec;
import static edu.emory.mathcs.csparsej.tdouble.Dcs_lsolve.cs_lsolve;
import static edu.emory.mathcs.csparsej.tdouble.Dcs_ltsolve.cs_ltsolve;
import static edu.emory.mathcs.csparsej.tdouble.Dcs_multiply.cs_multiply;
import static edu.emory.mathcs.csparsej.tdouble.Dcs_permute.cs_permute;
import static edu.emory.mathcs.csparsej.tdouble.Dcs_pinv.cs_pinv;
import static edu.emory.mathcs.csparsej.tdouble.Dcs_pvec.cs_pvec;
import static edu.emory.mathcs.csparsej.tdouble.Dcs_schol.cs_schol;
import static edu.emory.mathcs.csparsej.tdouble.Dcs_transpose.cs_transpose;
import static edu.emory.mathcs.csparsej.tdouble.Dcs_updown.cs_updown;
import static edu.emory.mathcs.csparsej.tdouble.Dcs_util.cs_spalloc;

import edu.emory.mathcs.csparsej.tdouble.Dcs_common.Dcs;
import edu.emory.mathcs.csparsej.tdouble.Dcs_common.Dcsn;
import edu.emory.mathcs.csparsej.tdouble.Dcs_common.Dcss;

/**
 * Read a matrix, solve a linear system, update/downdate.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class Dcs_test3 extends Dcs_test {

	/**
	 * Cholesky update/downdate
	 *
	 * @param prob
	 *            problem
	 * @return true if successful, false on error
	 */
	protected static boolean test3(Dproblem prob)
	{
		Dcs A, C, W = null, WW, WT, E = null, W2 ;
		int n, k, Li[], Lp[], Wi[], Wp[], p1, p2, p[] = null ;
		boolean ok ;
		double[] b, x, resid, y = null, Lx, Wx ;
		double s ;
		double t, t1 ;
		Dcss S = null ;
		Dcsn N = null ;
		if (prob == null || prob.sym == 0 || prob.A.n == 0) return (false) ;
		A = prob.A ; C = prob.C ; b = prob.b ; x = prob.x ; resid = prob.resid ;
		n = A.n ;
		if (prob.sym == 0 || n == 0) return (true) ;
		rhs (x, b, n) ;						/* compute right-hand side */
		System.out.print("\nchol then update/downdate ") ;
		print_order (1) ;
		y = new double[n] ;
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
		Lp = N.L.p ; Li = N.L.i ; Lx = N.L.x ;
		Wp = W.p ; Wi = W.i ; Wx = W.x ;
		Wp [0] = 0 ;
		p1 = Lp [k] ;
		Wp [1] = Lp [k+1] - p1 ;
		s = Lx[p1] ;
		Random r = new Random(1) ;
		for ( ; p1 < Lp [k+1] ; p1++)
		{
			p2 = p1 - Lp [k] ;
			Wi [p2] = Li [p1] ;
			Wx[p2] = s * r.nextDouble() ;
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
		E = cs_add (C, WW, 1, 1) ;
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

	public void test_bcsstk01()
	{
		InputStream in = get_stream (BCSSTK01) ;
		Dproblem prob = get_problem (in, 0) ;

		assert_problem(prob, 48, 48, 224, -1, 400, 3.5709480746974373e+09) ;

		test3(prob) ;

		double x_norm = 0.0005 ;
		assertEquals(x_norm, prob.norms.get(0), 1e-04) ;
		assertEquals(x_norm, prob.norms.get(1), 1e-04) ;
		assertEquals(x_norm, prob.norms.get(2), 1e-04) ;
		assertEquals(x_norm, prob.norms.get(3), 1e-04) ;
	}

	public void test_bcsstk16()
	{
		InputStream in = get_stream (BCSSTK16) ;
		Dproblem prob = get_problem (in, 0) ;

		assert_problem(prob, 4884, 4884, 147631, -1, 290378, 7.008379365769155e+09) ;

		test3(prob) ;

		double x_norm = 1.9998 ;
		assertEquals(x_norm, prob.norms.get(0), DELTA) ;
		assertEquals(x_norm, prob.norms.get(1), DELTA) ;
		assertEquals(x_norm, prob.norms.get(2), DELTA) ;
		assertEquals(x_norm, prob.norms.get(3), DELTA) ;
	}

}
