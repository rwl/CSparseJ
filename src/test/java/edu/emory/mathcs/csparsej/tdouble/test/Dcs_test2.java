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

import static edu.emory.mathcs.csparsej.tdouble.Dcs_cholsol.cs_cholsol ;
import static edu.emory.mathcs.csparsej.tdouble.Dcs_dmperm.cs_dmperm ;
import static edu.emory.mathcs.csparsej.tdouble.Dcs_lusol.cs_lusol ;
import static edu.emory.mathcs.csparsej.tdouble.Dcs_qrsol.cs_qrsol ;

import edu.emory.mathcs.csparsej.tdouble.Dcs_common.Dcs ;
import edu.emory.mathcs.csparsej.tdouble.Dcs_common.Dcsd ;

/**
 * Read a matrix from a file and solve a linear system.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class Dcs_test2 extends Dcs_test {

	/**
	 * Solves a linear system using Cholesky, LU, and QR, with various
	 * orderings.
	 *
	 * @param prob problem
	 * @return true if successful, false on error
	 */
	protected static boolean test2(Dproblem prob)
	{
		Dcs A, C ;
		double[] b, x, resid ;
		double t, tol ;
		int k, m, n, order, nb, ns, r[], s[], rr[], sprank ;
		boolean ok ;
		Dcsd D ;
		if (prob == null) return (false) ;
		A = prob.A ; C = prob.C ; b = prob.b ; x = prob.x ; resid = prob.resid ;
		m = A.m ; n = A.n ;
		tol = prob.sym != 0 ? 0.001 : 1 ;	/* partial pivoting tolerance */
		D = cs_dmperm (C, 1) ;			/* randomized dmperm analysis */
		if (D == null) return (false) ;
		prob.nb = nb = D.nb ; r = D.r ; s = D.s ; rr = D.rr ;
		prob.sprank = sprank = rr [3] ;
		for (ns = 0, k = 0 ; k < nb ; k++)
		{
			if ((r [k+1] == r [k] + 1) && (s [k+1] == s [k] + 1))
			{
				ns++ ;
			}
		}
		prob.ns = ns ;
		System.out.printf("blocks: %d singletons: %d structural rank: %d\n", nb, ns, sprank) ;
		D = null ;
		for (order = 0 ; order <= 3 ; order += 3)	/* natural and amd(A'*A) */
		{
			if (order == 0 && m > 1000) continue ;
			System.out.print("QR   ") ;
			print_order (order) ;
			rhs (x, b, m) ;				/* compute right-hand side */
			t = tic() ;
			ok = cs_qrsol (order, C, x) ;		/* min norm(Ax-b) with QR */
			System.out.printf("time: %8.2f ms ", toc (t)) ;
			print_resid (ok, C, x, b, resid, prob) ;	/* print residual */
		}
		if (m != n || sprank < n) return (true) ;	/* return if rect. or singular*/
		for (order = 0 ; order <= 3 ; order++)		/* try all orderings */
		{
			if (order == 0 && m > 1000) continue ;
			System.out.print("LU   ") ;
			print_order (order) ;
			rhs (x, b, m) ;				/* compute right-hand side */
			t = tic() ;
			ok = cs_lusol (order, C, x, tol) ;	/* solve Ax=b with LU */
			System.out.printf("time: %8.2f ms ", toc (t)) ;
			print_resid (ok, C, x, b, resid, prob) ;  /* print residual */
		}
		if (prob.sym == 0) return (true) ;
		for (order = 0 ; order <= 1 ; order++)		/* natural and amd(A+A') */
		{
			if (order == 0 && m > 1000) continue ;
			System.out.print("Chol ") ;
			print_order (order) ;
			rhs (x, b, m) ;				/* compute right-hand side */
			t = tic() ;
			ok = cs_cholsol (order, C, x) ;		/* solve Ax=b with Cholesky */
			System.out.printf("time: %8.2f ms ", toc (t)) ;
			print_resid (ok, C, x, b, resid, prob) ;  /* print residual */
		}
		return (true) ;
	}

	public void test_ash219()
	{
		InputStream in = get_stream (ASH219) ;
		Dproblem prob = get_problem (in, DROP_TOL) ;

		test2(prob) ;

		assert_problem(prob, 219, 85, 438, 0, 0, 9.0) ;
		assert_structure(prob, 1, 0, 85) ;

		assertEquals(1.0052, prob.norms.get(0), DELTA) ;
		assertEquals(1.0052, prob.norms.get(1), DELTA) ;
	}

	public void test_bcsstk01()
	{
		InputStream in = get_stream (BCSSTK01) ;
		Dproblem prob = get_problem (in, DROP_TOL) ;

		test2(prob) ;
		assert_problem(prob, 48, 48, 224, -1, 400, 3.57094807469e+09) ;
		assert_structure(prob, 1, 0, 48) ;

		double x_norm = 0.0005 ;
		assertEquals(x_norm, prob.norms.get(0), DELTA) ;
		assertEquals(x_norm, prob.norms.get(1), DELTA) ;
		assertEquals(x_norm, prob.norms.get(2), DELTA) ;
		assertEquals(x_norm, prob.norms.get(3), DELTA) ;
		assertEquals(x_norm, prob.norms.get(4), DELTA) ;
		assertEquals(x_norm, prob.norms.get(5), DELTA) ;
		assertEquals(x_norm, prob.norms.get(6), DELTA) ;
		assertEquals(x_norm, prob.norms.get(7), DELTA) ;
	}

	public void test_bcsstk16()
	{
		InputStream in = get_stream (BCSSTK16) ;
		Dproblem prob = get_problem (in, DROP_TOL) ;

		test2(prob) ;

		assert_problem(prob, 4884, 4884, 147631, -1, 290378, 7.008379365769155e+09) ;
		assert_structure(prob, 75, 74, 4884) ;

		double x_norm = 1.9998 ;
		assertEquals(x_norm, prob.norms.get(0), DELTA) ;
		assertEquals(x_norm, prob.norms.get(1), DELTA) ;
		assertEquals(x_norm, prob.norms.get(2), DELTA) ;
		assertEquals(x_norm, prob.norms.get(3), DELTA) ;
		assertEquals(x_norm, prob.norms.get(4), DELTA) ;
	}

	public void test_fs_183_1()
	{
		InputStream in = get_stream (FS_183_1) ;
		Dproblem prob = get_problem (in, DROP_TOL) ;

		test2(prob) ;

		assert_problem(prob, 183, 183, 988, 0, 0, 1.7031774210073e+09) ;
		assert_dropped(prob, 71, 10) ;
		assert_structure(prob, 38, 37, 183) ;

		double x_norm = 212022.2099 ;
		assertEquals(x_norm, prob.norms.get(0), DELTA) ;
		assertEquals(x_norm, prob.norms.get(1), DELTA) ;
		assertEquals(x_norm, prob.norms.get(2), DELTA) ;
		assertEquals(x_norm, prob.norms.get(3), DELTA) ;
		assertEquals(x_norm, prob.norms.get(4), DELTA) ;
		assertEquals(x_norm, prob.norms.get(5), DELTA) ;
	}

	public void test_ibm32a()
	{
		InputStream in = get_stream (IBM32A) ;
		Dproblem prob = get_problem (in, DROP_TOL) ;

		test2(prob) ;

		assert_problem(prob, 32, 31, 123, 0, 0, 7.0) ;
		assert_structure(prob, 1, 0, 31) ;

		double x_norm = 5.5800 ;
		assertEquals(x_norm, prob.norms.get(0), DELTA) ;
		assertEquals(x_norm, prob.norms.get(1), DELTA) ;
	}

	public void test_ibm32b()
	{
		InputStream in = get_stream (IBM32B) ;
		Dproblem prob = get_problem (in, DROP_TOL) ;

		test2(prob) ;

		assert_problem(prob, 31, 32, 123, 0, 0, 8.0) ;
		assert_structure(prob, 1, 0, 31) ;

		double x_norm = 5.3348 ;
		assertEquals(x_norm, prob.norms.get(0), DELTA) ;
		assertEquals(x_norm, prob.norms.get(1), DELTA) ;
	}

	public void test_lp_afiro()
	{
		InputStream in = get_stream (LP_AFIRO) ;
		Dproblem prob = get_problem (in, DROP_TOL) ;

		test2(prob) ;

		assert_problem(prob, 27, 51, 102, 0, 0, 3.43) ;
		assert_structure(prob, 1, 0, 27) ;

		assertEquals(2.4534, prob.norms.get(0), DELTA) ;
		assertEquals(2.4534, prob.norms.get(1), DELTA) ;
	}

	public void test_mbeacxc()
	{
		InputStream in = get_stream (MBEACXC) ;
		Dproblem prob = get_problem (in, DROP_TOL) ;

		test2(prob) ;

		assert_problem(prob, 492, 490, 49920, 0, 0, 9.29e-01) ;
		assert_structure(prob, 10, 8, 448) ;

		assertEquals(Double.NaN, prob.norms.get(0), DELTA) ;
		assertEquals(Double.NaN, prob.norms.get(1), DELTA) ;
	}

	public void test_t1()
	{
		InputStream in = get_stream (T1) ;
		Dproblem prob = get_problem (in, DROP_TOL) ;

		test2(prob) ;

		assert_problem(prob, 4, 4, 10, 0, 0, 1.11e+01) ;
		assert_structure(prob, 1, 0, 4) ;

		double x_norm = 2.4550 ;
		assertEquals(x_norm, prob.norms.get(0), DELTA) ;
		assertEquals(x_norm, prob.norms.get(1), DELTA) ;
		assertEquals(x_norm, prob.norms.get(2), DELTA) ;
		assertEquals(x_norm, prob.norms.get(3), DELTA) ;
		assertEquals(x_norm, prob.norms.get(4), DELTA) ;
		assertEquals(x_norm, prob.norms.get(5), DELTA) ;
	}

	public void test_west0067()
	{
		InputStream in = get_stream (WEST0067) ;
		Dproblem prob = get_problem (in, DROP_TOL) ;

		test2(prob) ;

		assert_problem(prob, 67, 67, 294, 0, 0, 6.14) ;
		assert_structure(prob, 2, 1, 67) ;

		double x_norm = 21.9478 ;
		assertEquals(x_norm, prob.norms.get(0), DELTA) ;
		assertEquals(x_norm, prob.norms.get(1), DELTA) ;
		assertEquals(x_norm, prob.norms.get(2), DELTA) ;
		assertEquals(x_norm, prob.norms.get(3), DELTA) ;
		assertEquals(x_norm, prob.norms.get(4), DELTA) ;
		assertEquals(x_norm, prob.norms.get(5), DELTA) ;
	}

}
