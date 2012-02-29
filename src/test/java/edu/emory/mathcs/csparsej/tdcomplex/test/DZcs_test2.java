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

import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_cholsol.cs_cholsol ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_dmperm.cs_dmperm ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_lusol.cs_lusol ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_qrsol.cs_qrsol ;

import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcs ;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcsa ;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcsd ;

/**
 * Read a matrix from a file and solve a linear system.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class DZcs_test2 extends DZcs_test {

	/**
	 * Solves a linear system using Cholesky, LU, and QR, with various
	 * orderings.
	 *
	 * @param prob
	 *            problem
	 * @return true if successful, false on error
	 */
	protected static boolean test2(DZproblem prob)
	{
		DZcs A, C ;
		DZcsa b, x, resid ;
		double t, tol ;
		int k, m, n, order, nb, ns, r[], s[], rr[], sprank ;
		boolean ok ;
		DZcsd D ;
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

	public void test_c_ibm32a()
	{
		InputStream in = get_stream (C_IBM32A) ;
		DZproblem prob = get_problem (in, DROP_TOL) ;

		test2(prob) ;

		assert_problem(prob, 32, 31, 123, 0, 0, 9.90) ;
		assert_structure(prob, 1, 0, 31) ;

		assertEquals(1.92e-03, prob.residuals.get(0), DELTA) ;
		assertEquals(1.92e-03, prob.residuals.get(1), DELTA) ;
	}

	public void test_c_ibm32b()
	{
		InputStream in = get_stream (C_IBM32B) ;
		DZproblem prob = get_problem (in, DROP_TOL) ;

		test2(prob) ;

		assert_problem(prob, 31, 32, 123, 0, 0, 1.13e+01) ;
		assert_structure(prob, 1, 0, 31) ;

		assertEquals(1.22e-16, prob.residuals.get(0), DELTA) ;
		assertEquals(9.45e-17, prob.residuals.get(1), DELTA) ;
	}

	public void test_c_mbeacxc()
	{
		InputStream in = get_stream (C_MBEACXC) ;
		DZproblem prob = get_problem (in, DROP_TOL) ;

		test2(prob) ;

		assert_problem(prob, 492, 490, 49920, 0, 0, 9.29e-01) ;
		assert_structure(prob, 10, 8, 448) ;

		assertNull(prob.residuals.get(0)) ;
		assertNull(prob.residuals.get(1)) ;
	}

	public void test_c_west0067()
	{
		InputStream in = get_stream (C_WEST0067) ;
		DZproblem prob = get_problem (in, DROP_TOL) ;

		test2(prob) ;

		assert_problem(prob, 67, 67, 294, 0, 0, 6.17) ;
		assert_structure(prob, 2, 1, 67) ;

		assertEquals(8.16e-17, prob.residuals.get(0), DELTA) ;
		assertEquals(5.34e-17, prob.residuals.get(1), DELTA) ;
		assertEquals(4.54e-17, prob.residuals.get(2), DELTA) ;
		assertEquals(5.90e-17, prob.residuals.get(3), DELTA) ;
		assertEquals(5.28e-17, prob.residuals.get(4), DELTA) ;
		assertEquals(4.76e-17, prob.residuals.get(5), DELTA) ;
	}

	public void test_c4()
	{
		InputStream in = get_stream (C4) ;
		DZproblem prob = get_problem (in, DROP_TOL) ;

		test2(prob) ;

		assert_problem(prob, 4, 4, 10, -1, 16, 7.37e+01) ;
		assert_structure(prob, 1, 0, 4) ;

		assertEquals(5.85e-17, prob.residuals.get(0), DELTA) ;
		assertEquals(5.85e-17, prob.residuals.get(1), DELTA) ;

		assertEquals(2.29e-17, prob.residuals.get(2), DELTA) ;
		assertEquals(2.29e-17, prob.residuals.get(3), DELTA) ;
		assertEquals(2.29e-17, prob.residuals.get(4), DELTA) ;
		assertEquals(2.29e-17, prob.residuals.get(5), DELTA) ;

		assertEquals(6.88e-17, prob.residuals.get(6), DELTA) ;
		assertEquals(6.88e-17, prob.residuals.get(7), DELTA) ;
	}

	public void test_czero()
	{
		InputStream in = get_stream (CZERO) ;
		DZproblem prob = get_problem (in, DROP_TOL) ;

		test2(prob) ;

		assert_problem(prob, 1, 1, 0, 1, 0, 0.0) ;
		assert_dropped(prob, 1, 0) ;
		assert_structure(prob, 2, 0, 0) ;

		assertNull(prob.residuals.get(0)) ;
		assertNull(prob.residuals.get(1)) ;
	}

	public void test_mhd1280b()
	{
		InputStream in = get_stream (MHD1280B) ;
		DZproblem prob = get_problem (in, DROP_TOL) ;

		test2(prob) ;

		assert_problem(prob, 1280, 1280, 11963, -1, 22646, 8.00e+01) ;
		assert_dropped(prob, 0, 66) ;
		assert_structure(prob, 20, 14, 1280) ;

		assertEquals(6.15e-25, prob.residuals.get(0), DELTA) ;
		assertEquals(2.33e-25, prob.residuals.get(1), DELTA) ;
		assertEquals(3.96e-25, prob.residuals.get(2), DELTA) ;
		assertEquals(3.96e-25, prob.residuals.get(3), DELTA) ;
		assertEquals(1.58e-25, prob.residuals.get(4), DELTA) ;
	}

	public void test_neumann()
	{
		InputStream in = get_stream (NEUMANN) ;
		DZproblem prob = get_problem (in, DROP_TOL) ;

		test2(prob) ;

		assert_problem(prob, 1600, 1600, 7840, 0, 0, 1.41e+01) ;
		assert_structure(prob, 1, 0, 1600) ;

		assertEquals(1.04e-15, prob.residuals.get(0), DELTA) ;
		assertEquals(3.55e-16, prob.residuals.get(1), DELTA) ;
		assertEquals(4.03e-16, prob.residuals.get(2), DELTA) ;
		assertEquals(4.03e-16, prob.residuals.get(3), DELTA) ;
	}

	public void test_qc324()
	{
		InputStream in = get_stream (QC324) ;
		DZproblem prob = get_problem (in, DROP_TOL) ;

		test2(prob) ;

		assert_problem(prob, 324, 324, 26730, 0, 0, 1.71) ;
		assert_structure(prob, 1, 0, 324) ;

		assertEquals(9.42e-17, prob.residuals.get(0), DELTA) ;
		assertEquals(8.94e-17, prob.residuals.get(1), DELTA) ;
		assertEquals(6.01e-17, prob.residuals.get(2), DELTA) ;
		assertEquals(4.05e-17, prob.residuals.get(3), DELTA) ;
		assertEquals(4.71e-17, prob.residuals.get(4), DELTA) ;
		assertEquals(4.71e-17, prob.residuals.get(5), DELTA) ;
	}

	public void test_t2()
	{
		InputStream in = get_stream (T2) ;
		DZproblem prob = get_problem (in, DROP_TOL) ;

		test2(prob) ;

		assert_problem(prob, 4, 4, 10, 0, 0, 1.06e+02) ;
		assert_structure(prob, 1, 0, 4) ;

		assertEquals(2.06e-17, prob.residuals.get(0), DELTA) ;
		assertEquals(6.36e-18, prob.residuals.get(1), DELTA) ;
		assertEquals(4.88e-18, prob.residuals.get(2), DELTA) ;
		assertEquals(2.11e-18, prob.residuals.get(3), DELTA) ;
		assertEquals(6.36e-18, prob.residuals.get(4), DELTA) ;
		assertEquals(2.11e-18, prob.residuals.get(5), DELTA) ;
	}

	public void test_t3()
	{
		InputStream in = get_stream (T3) ;
		DZproblem prob = get_problem (in, DROP_TOL) ;

		test2(prob) ;

		assert_problem(prob, 3, 4, 12, 0, 0, 3.06) ;
		assert_structure(prob, 1, 0, 3) ;

		assertEquals(1.05e-16, prob.residuals.get(0), DELTA) ;
		assertEquals(1.05e-16, prob.residuals.get(1), DELTA) ;
	}

	public void test_t4()
	{
		InputStream in = get_stream (T4) ;
		DZproblem prob = get_problem (in, DROP_TOL) ;

		test2(prob) ;

		assert_problem(prob, 2, 2, 3, 1, 4, 2.83) ;
		assert_structure(prob, 1, 0, 2) ;

		assertEquals(5.65e-17, prob.residuals.get(0), DELTA) ;
		assertEquals(5.65e-17, prob.residuals.get(1), DELTA) ;
		assertEquals(0.0, prob.residuals.get(2), DELTA) ;
		assertEquals(0.0, prob.residuals.get(3), DELTA) ;
		assertEquals(0.0, prob.residuals.get(4), DELTA) ;
		assertEquals(0.0, prob.residuals.get(5), DELTA) ;
		assertNull(prob.residuals.get(6)) ;
		assertNull(prob.residuals.get(7)) ;
	}

	public void test_young1c()
	{
		InputStream in = get_stream (YOUNG1C) ;
		DZproblem prob = get_problem (in, DROP_TOL) ;

		test2(prob) ;

		assert_problem(prob, 841, 841, 4089, 0, 0, 7.30e+02) ;
		assert_structure(prob, 1, 0, 841) ;

		assertEquals(1.81e-16, prob.residuals.get(0), DELTA) ;
		assertEquals(1.57e-16, prob.residuals.get(1), DELTA) ;
		assertEquals(1.39e-16, prob.residuals.get(2), DELTA) ;
		assertEquals(2.95e-16, prob.residuals.get(3), DELTA) ;
		assertEquals(3.37e-16, prob.residuals.get(4), DELTA) ;
		assertEquals(3.37e-16, prob.residuals.get(5), DELTA) ;
	}

}
