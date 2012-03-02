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

import junit.framework.TestCase;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_add.cs_add ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_compress.cs_compress ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_droptol.cs_droptol ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_dropzeros.cs_dropzeros ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_dupl.cs_dupl ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_fkeep.cs_fkeep ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_gaxpy.cs_gaxpy ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_load.cs_load ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_norm.cs_norm ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_transpose.cs_transpose ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex.cs_cone ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex.cs_cabs ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex.cs_cneg ;

import edu.emory.mathcs.csparsej.tdcomplex.DZcs_ifkeep;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcs ;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcsa ;

/**
 * Support routines for DZcs_test*.java
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
abstract public class DZcs_test extends TestCase {

	protected static final double DELTA = 1e-4;
	protected static final double DROP_TOL = 1e-14;

	protected static final String DIR = "matrix";

	protected static final String CZERO = "czero";
	protected static final String C4 = "c4";
	protected static final String T2 = "t2";
	protected static final String T3 = "t3";
	protected static final String T4 = "t4";

	/**
	 * US economy, 1972.  Dan Szyld, while at NYU.
	 * Numerically and structurally singular
	 */
	protected static final String C_MBEACXC = "c_mbeacxc";

	/** Cavett problem with 5 components (chemical eng., Westerberg) */
	protected static final String C_WEST0067 = "c_west0067";

	/** Alfven spectra in magnetohydrodynamics */
	protected static final String MHD1280B = "mhd1280b";

	/** Model of H+ in an electromagnetic field */
	protected static final String QC324 = "qc324";

	/** The Harwell/Boeing matrix ibm32, but with the last column removed. */
	protected static final String C_IBM32A = "c_ibm32a";

	/** The transpose of ibm32a */
	protected static final String C_IBM32B = "c_ibm32b";

	/** Aeronautical problem */
	protected static final String YOUNG1C = "young1c";


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
		assert_dimensions (A, m, n, nzmax, nnz);
		assertEquals (norm1, cs_norm (A), DELTA);
	}

	protected static void assert_dimensions(DZcs A, int m, int n, int nzmax, int nnz) {
		assertEquals (m, A.m);
		assertEquals (n, A.n);
		assertEquals (nzmax, A.nzmax);

		int nz = (A.nz < 0) ? A.p [A.n] : A.nz ;
		assertEquals (nnz, nz);
	}

	protected static void assert_problem(DZproblem prob, int m, int n, int nnz, int sym, int sym_nnz,
			double norm) {
		assertEquals (m, prob.A.m) ;
		assertEquals (n, prob.A.n) ;
		assertEquals (nnz, prob.A.p [n]) ;
		assertEquals (sym, prob.sym) ;
		assertEquals (sym_nnz, sym != 0 ? prob.C.p [n] : 0) ;
		assertEquals (norm, cs_norm (prob.C), 1e-2) ;
	}

	protected static void assert_structure(DZproblem prob, int blocks, int singletons, int rank) {
		assertEquals (blocks, prob.nb) ;
		assertEquals (singletons, prob.ns) ;
		assertEquals (rank, prob.sprank) ;
	}

	protected static void assert_dropped(DZproblem prob, int dropped_zeros, int dropped_tiny) {
		assertEquals (dropped_zeros, prob.dropped_zeros) ;
		assertEquals (dropped_tiny, prob.dropped_tiny) ;
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

		public List<Double> norms = new ArrayList<Double>() ;

		public int nb ;
		public int ns ;
		public int sprank ;

		public int dropped_zeros ;
		public int dropped_tiny ;

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
		System.out.printf ("resid: %8.2e", r) ;

		double nrm = norm (x, n) ;
		System.out.printf (" (norm: %8.4f, %8.4f)\n", nrm, norm (b, m)) ;
		prob.norms.add (nrm) ;
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
		if (tol > 0) cs_dropzeros (A) ;			/* drop zero entries */
		nz2 = A.p [n] ;
		if (tol > 0) cs_droptol (A, tol) ;		/* drop tiny entries (just to test) */
		prob.C = C = sym != 0 ? make_sym(A) : A ;	/* C = A + triu(A,1)', or C=A */
		if (C == null) return (null) ;
		System.out.printf("\n--- Matrix: %d-by-%d, nnz: %d (sym: %d: nnz %d), norm: %8.2e\n",
			m, n, A.p [n], sym, sym != 0 ? C.p [n] : 0, cs_norm (C)) ;
		prob.dropped_zeros = nz1 - nz2 ;
		if (nz1 != nz2) System.out.printf("zero entries dropped: %d\n", nz1 - nz2) ;
		prob.dropped_tiny = nz2 - A.p [n] ;
		if (nz2 != A.p [n]) System.out.printf("tiny entries dropped: %d\n", nz2 - A.p [n]) ;
		prob.b = new DZcsa (mn) ;
		prob.x = new DZcsa (mn) ;
		prob.resid = new DZcsa (mn) ;
		return prob ;
	}

}
