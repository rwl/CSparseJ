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

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random ;

import junit.framework.TestCase;

import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_add.cs_add ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_chol.cs_chol ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_compress.cs_compress ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_droptol.cs_droptol ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_dropzeros.cs_dropzeros ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_dupl.cs_dupl ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_fkeep.cs_fkeep ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_gaxpy.cs_gaxpy ;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_ifkeep;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_ipvec.cs_ipvec ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_load.cs_load ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_lsolve.cs_lsolve ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_ltsolve.cs_ltsolve ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_multiply.cs_multiply ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_norm.cs_norm ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_permute.cs_permute ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_pinv.cs_pinv ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_pvec.cs_pvec ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_schol.cs_schol ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_transpose.cs_transpose ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_updown.cs_updown ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_util.cs_spalloc ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex.cs_cone ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex.cs_cabs ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex.cs_cneg ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex.cs_cmult ;

import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcs ;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcsa ;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcsn ;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcss ;

/**
 * Support routines for Zcs_demo*.java
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
abstract public class DZcs_test extends TestCase {

	protected static final double DELTA = 1e-5;
	protected static final double DROP_TOL = 1e-14;

	protected static final String C_IBM32A = "c_ibm32a";
	protected static final String C_IBM32B = "c_ibm32b";
	protected static final String C_MBEACXC = "c_mbeacxc";
	protected static final String C_WEST0067 = "c_west0067";
	protected static final String C4 = "c4";
	protected static final String CZERO = "czero";
	protected static final String MHD1280B = "mhd1280b";
	protected static final String NEUMANN = "neumann";
	protected static final String QC324 = "qc324";
	protected static final String T2 = "t2";
	protected static final String T3 = "t3";
	protected static final String YOUNG1C = "young1c";

	protected static final String DIR = "matrix";

	protected static InputStream get_stream(String name) {
		try
		{
			return DZcs_test1.class.getResource(DIR + "/" + name).openStream() ;
		}
		catch (IOException e)
		{
			return (null) ;
		}
	}

	protected static void assert_dimensions(DZcs A, int m, int n, int nzmax, int nnz, double norm1) {
		assert_dimensions(A, m, n, nzmax, nnz);
		assertEquals(norm1, cs_norm (A), DELTA);
	}

	protected static void assert_dimensions(DZcs A, int m, int n, int nzmax, int nnz) {
		assertEquals(m, A.m);
		assertEquals(n, A.n);
		assertEquals(nzmax, A.nzmax);

		int nz = (A.nz < 0) ? A.p [A.n] : A.nz ;
		assertEquals(nnz, nz);
	}

	/**
	 *
	 * A structure for a demo problem.
	 *
	 */
	protected static class DZproblem {

		public DZcs A ;
		public DZcs C ;
		public int sym ;
		public DZcsa x ;
		public DZcsa b ;
		public DZcsa resid ;

		public List<Double> residuals = new ArrayList<Double>();

		public int nb ;
		public int ns ;
		public int sprank ;

		public DZproblem () {}

	};

	/**
	 * 1 if A is square & upper tri., -1 if square & lower tri., 0 otherwise
	 */
	protected static int is_sym(DZcs A)
	{
		int j, p, n = A.n, m = A.m, Ap[] = A.p, Ai[] = A.i ;
		boolean is_upper, is_lower ;
		if (m != n) return (0) ;
		is_upper = true ;
		is_lower = true ;
		for (j = 0 ; j < n ; j++)
		{
			for (p = Ap [j] ; p < Ap [j+1] ; p++)
			{
				if (Ai [p] > j) is_upper = false ;
				if (Ai [p] < j) is_lower = false ;
			}
		}
		return (is_upper ? 1 : (is_lower ? -1 : 0)) ;
	}

	/**
	 * true for off-diagonal entries
	 */
	protected static class Dropdiag implements DZcs_ifkeep {

		public boolean fkeep(int i, int j, double [] aij, Object other)
		{
			return (i != j) ;
		}

	}

	/**
	 * C = A + triu(A,1)'
	 */
	protected static DZcs make_sym(DZcs A)
	{
		DZcs AT, C ;
		AT = cs_transpose (A, true) ; 			/* AT = A' */
		cs_fkeep (AT, new Dropdiag(), null) ;		/* drop diagonal entries from AT */
		C = cs_add (A, AT, cs_cone(), cs_cone()) ;	/* C = A+AT */
		AT = null ;
		return (C) ;
	}

	/**
	 * create a right-hand side
	 */
	protected static void rhs(DZcsa x, DZcsa b, int m)
	{
		int i;
		for (i = 0; i < m; i++) b.set(i, new double[] {1 + ((double) i) / m, 0.0}) ;
		for (i = 0; i < m; i++) x.set(i, b.get(i)) ;
	}

	/**
	 * infinity-norm of x
	 */
	protected static double norm(DZcsa x, int n)
	{
		int i ;
		double normx = 0 ;
		for (i = 0 ; i < n ; i++) normx = Math.max(normx, cs_cabs (x.get(i))) ;
		return (normx) ;
	}

	/**
	 * compute residual, norm(A*x-b,inf) / (norm(A,1)*norm(x,inf) + norm(b,inf))
	 */
	protected static void print_resid(boolean ok, DZcs A, DZcsa x, DZcsa b, DZcsa resid, DZproblem prob)
	{
		int i, m, n ;
		if (!ok)
		{
		    System.out.print("    (failed)\n") ;
		    return ;
		}
		m = A.m ; n = A.n ;
		for (i = 0 ; i < m ; i++)
			resid.set(i, cs_cneg(b.get(i))) ;	/* resid = -b */
		cs_gaxpy (A, x, resid) ;			/* resid = resid + A*x  */

		double r = norm (resid, m) / ((n == 0) ? 1 : (cs_norm (A) * norm (x, n) + norm (b, m))) ;
		prob.residuals.add(r) ;

		System.out.printf("resid: %8.2e\n", r) ;
	}

	protected static double tic()
	{
		return System.nanoTime();
	}

	protected static double toc(double t)
	{
		double s = tic() ;
		return (Math.max(0, s - t)) / 1000000.0 ;
	}

	protected static void print_order(int order)
	{
		switch (order)
		{
		case 0:
		    System.out.print ("natural    ") ;
		    break ;
		case 1:
		    System.out.print ("amd(A+A')  ") ;
		    break ;
		case 2:
		    System.out.print ("amd(S'*S)  ") ;
		    break ;
		case 3:
		    System.out.print ("amd(A'*A)  ") ;
		    break ;
		}
	}

	/**
	 * Reads a problem from a file.
	 *
	 * @param fileName
	 *            file name
	 * @param tol
	 *            drop tolerance
	 * @return problem
	 */
	protected static DZproblem get_problem(InputStream in, double tol)
	{
		DZcs T, A, C ;
		int sym, m, n, mn, nz1, nz2 ;
		DZproblem prob ;
		prob = new DZproblem() ;
		T = cs_load (in) ;				/* load triplet matrix T from a file */
		prob.A = A = cs_compress (T) ;			/* A = compressed-column form of T */
		T = null ;					/* clear T */
		if (!cs_dupl (A)) return (null) ;		/* sum up duplicates */
		prob.sym = sym = is_sym (A) ;			/* determine if A is symmetric */
		m = A.m ; n = A.n ;
		mn = Math.max (m, n) ;
		nz1 = A.p [n] ;
		cs_dropzeros (A) ;				/* drop zero entries */
		nz2 = A.p [n] ;
		if (tol > 0) cs_droptol (A, tol) ;		/* drop tiny entries (just to test) */
		prob.C = C = sym != 0 ? make_sym(A) : A ;	/* C = A + triu(A,1)', or C=A */
		if (C == null) return (null) ;
		System.out.printf("\n--- Matrix: %d-by-%d, nnz: %d (sym: %d: nnz %d), norm: %8.2e\n",
			m, n, A.p [n], sym, sym != 0 ? C.p [n] : 0, cs_norm (C)) ;
		if (nz1 != nz2) System.out.printf("zero entries dropped: %d\n", nz1 - nz2) ;
		if (nz2 != A.p [n]) System.out.printf("tiny entries dropped: %d\n", nz2 - A.p [n]) ;
		prob.b = new DZcsa (mn) ;
		prob.x = new DZcsa (mn) ;
		prob.resid = new DZcsa (mn) ;
		return prob ;
	}

	protected static void assert_problem(DZproblem prob, int m, int n, int nnz, int sym, int sym_nnz,
			double norm) {
		assertEquals(m, prob.A.m) ;
		assertEquals(n, prob.A.n) ;
		assertEquals(nnz, prob.A.p [n]) ;
		assertEquals(sym, prob.sym) ;
		assertEquals(sym_nnz, sym != 0 ? prob.C.p [n] : 0) ;
		assertEquals(norm, cs_norm (prob.C), 1e-2) ;
	}

	protected static void assert_structure(DZproblem prob, int blocks, int singletons, int rank) {
		assertEquals(blocks, prob.nb) ;
		assertEquals(singletons, prob.ns) ;
		assertEquals(rank, prob.sprank) ;
	}

	/**
	 * Cholesky update/downdate
	 *
	 * @param Prob
	 *            problem
	 * @return true if successful, false on error
	 */
	protected static boolean demo3(DZproblem Prob)
	{
		DZcs A, C, W = null, WW, WT, E = null, W2 ;
		int n, k, Li[], Lp[], Wi[], Wp[], p1, p2, p[] = null ;
		boolean ok ;
		DZcsa b, x, resid, y = null, Lx, Wx ;
		double[] s ;
		double t, t1 ;
		DZcss S = null ;
		DZcsn N = null ;
		if (Prob == null || Prob.sym == 0 || Prob.A.n == 0) return (false) ;
		A = Prob.A ; C = Prob.C ; b = Prob.b ; x = Prob.x ; resid = Prob.resid ;
		n = A.n ;
		if (Prob.sym == 0 || n == 0) return (true) ;
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
		print_resid (true, C, x, b, resid, Prob) ;			/* print residual */
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
		print_resid (true, E, x, b, resid, Prob) ;		/* print residual */
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
		print_resid (true, E, x, b, resid, Prob) ;		/* print residual */
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
		print_resid (true, C, x, b, resid, Prob) ;		/* print residual */
		return (true) ;
	}

}
