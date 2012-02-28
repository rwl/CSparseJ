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

}
