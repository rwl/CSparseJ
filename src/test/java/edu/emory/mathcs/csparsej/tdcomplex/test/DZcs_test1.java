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

package edu.emory.mathcs.csparsej.tdcomplex.test;

import java.io.InputStream;

import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_add.cs_add ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_compress.cs_compress ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_entry.cs_entry ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_load.cs_load ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_multiply.cs_multiply ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_norm.cs_norm ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_transpose.cs_transpose ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_util.cs_spalloc ;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex.cs_cone ;
//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_print.cs_print ;

import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcs ;

/**
 * Read a matrix from a file and perform basic matrix operations.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class DZcs_test1 extends DZcs_test {

	protected static DZcs demo1_load(InputStream in) {
		DZcs T = cs_load (in) ;			/* load triplet matrix T from file */
		//System.out.print("T:\n") ;
		//cs_print (T, false) ;			/* print T */
		return T ;
	}

	protected static DZcs demo1_compress(DZcs T) {
		DZcs A = cs_compress (T) ;		/* A = compressed-column form of T */
		//System.out.print("A:\n") ;
		//cs_print (A, false) ;			/* print A */
		return A ;
	}

	protected static DZcs demo1_transpose(DZcs A) {
		DZcs AT = cs_transpose (A, true) ;	/* AT = A' */
		//System.out.print("AT:\n") ;
		//cs_print (AT, false) ;			/* print AT */
		return AT ;
	}

	protected static DZcs demo1_multiply_add(DZcs A, DZcs AT) {
		DZcs T, Eye, C, D ;
		int m = A != null ? A.m : 0 ;		/* m = # of rows of A */
		T = cs_spalloc (m, m, m, true, true);	/* create triplet identity matrix */
		for (int i = 0; i < m; i++) cs_entry (T, i, i, cs_cone()) ;
		Eye = cs_compress (T) ;			/* Eye = speye (m) */
		T = null ;
		C = cs_multiply (A, AT) ;		/* C = A*A' */
		D = cs_add(C, Eye, cs_cone(),
			new double[] {cs_norm (C), 0.0}) ;  /* D = C + Eye*norm (C,1) */
		//System.out.print("D:\n") ;
		//cs_print(D, false) ;			/* print D */
		return D ;
	}

	public void test_demo1_c_ibm32a() {
		DZcs T, A, AT, D ;

		InputStream in = getStream(C_IBM32A) ;

		T = demo1_load (in) ;
		assertDimensions(T, 32, 31, 128, 123);

		A = demo1_compress(T) ;
		assertDimensions(A, 32, 31, 123, 123, 9.89949) ;

		AT = demo1_transpose(A) ;
		assertDimensions(AT, 31, 32, 123, 123, 11.3137) ;

		D = demo1_multiply_add(A, AT) ;
		assertDimensions(D, 32, 32, 386, 386, 140) ;
	}

}
