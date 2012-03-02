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

		InputStream in = get_stream(C_IBM32A) ;

		T = demo1_load (in) ;
		assert_dimensions(T, 32, 31, 128, 123);

		A = demo1_compress(T) ;
		assert_dimensions(A, 32, 31, 123, 123, 9.89949) ;

		AT = demo1_transpose(A) ;
		assert_dimensions(AT, 31, 32, 123, 123, 11.3137) ;

		D = demo1_multiply_add(A, AT) ;
		assert_dimensions(D, 32, 32, 386, 386, 140) ;
	}

	public void test_demo1_c_ibm32b() {
		DZcs T, A, AT, D ;

		InputStream in = get_stream(C_IBM32B) ;

		T = demo1_load (in) ;
		assert_dimensions(T, 31, 32, 128, 123);

		A = demo1_compress(T) ;
		assert_dimensions(A, 31, 32, 123, 123, 11.3137) ;

		AT = demo1_transpose(A) ;
		assert_dimensions(AT, 32, 31, 123, 123, 9.89949) ;

		D = demo1_multiply_add(A, AT) ;
		assert_dimensions(D, 31, 31, 373, 373, 128) ;
	}

	public void test_demo1_c_mbeacxc() {
		DZcs T, A, AT, D ;

		InputStream in = get_stream(C_MBEACXC) ;

		T = demo1_load (in) ;
		assert_dimensions(T, 492, 490, 65536, 49920);

		A = demo1_compress(T) ;
		assert_dimensions(A, 492, 490, 49920, 49920, 0.92863) ;

		AT = demo1_transpose(A) ;
		assert_dimensions(AT, 490, 492, 49920, 49920, 16.5516) ;

		D = demo1_multiply_add(A, AT) ;
		assert_dimensions(D, 492, 492, 157350, 157350, 19.6068) ;
	}

	public void test_demo1_c_west0067() {
		DZcs T, A, AT, D ;

		InputStream in = get_stream(C_WEST0067) ;

		T = demo1_load (in) ;
		assert_dimensions(T, 67, 67, 512, 299);

		A = demo1_compress(T) ;
		assert_dimensions(A, 67, 67, 299, 299, 6.16948) ;

		AT = demo1_transpose(A) ;
		assert_dimensions(AT, 67, 67, 299, 299, 6.62541) ;

		D = demo1_multiply_add(A, AT) ;
		assert_dimensions(D, 67, 67, 1041, 1041, 61.4518) ;
	}

	public void test_demo1_c4() {
		DZcs T, A, AT, D ;

		InputStream in = get_stream(C4) ;

		T = demo1_load (in) ;
		assert_dimensions(T, 4, 4, 16, 10);

		A = demo1_compress(T) ;
		assert_dimensions(A, 4, 4, 10, 10, 58.8904) ;

		AT = demo1_transpose(A) ;
		assert_dimensions(AT, 4, 4, 10, 10, 66.8771) ;

		D = demo1_multiply_add(A, AT) ;
		assert_dimensions(D, 4, 4, 16, 16, 5029.4944) ;
	}

	public void test_demo1_czero() {
		DZcs T, A, AT, D ;

		InputStream in = get_stream(CZERO) ;

		T = demo1_load (in) ;
		assert_dimensions(T, 1, 1, 1, 1);

		A = demo1_compress(T) ;
		assert_dimensions(A, 1, 1, 1, 1, 0.0) ;

		AT = demo1_transpose(A) ;
		assert_dimensions(AT, 1, 1, 1, 1, 0.0) ;

		D = demo1_multiply_add(A, AT) ;
		assert_dimensions(D, 1, 1, 1, 1, 0.0) ;
	}

	public void test_demo1_mhd1280b() {
		DZcs T, A, AT, D ;

		InputStream in = get_stream(MHD1280B) ;

		T = demo1_load (in) ;
		assert_dimensions(T, 1280, 1280, 16384, 12029);

		A = demo1_compress(T) ;
		assert_dimensions(A, 1280, 1280, 12029, 12029, 79.974) ;

		AT = demo1_transpose(A) ;
		assert_dimensions(AT, 1280, 1280, 12029, 12029, 64.1995) ;

		D = demo1_multiply_add(A, AT) ;
		assert_dimensions(D, 1280, 1280, 27642, 27642, 8519.6629) ;
	}

	public void test_demo1_qc324() {
		DZcs T, A, AT, D ;

		InputStream in = get_stream(QC324) ;

		T = demo1_load (in) ;
		assert_dimensions(T, 324, 324, 32768, 26730);

		A = demo1_compress(T) ;
		assert_dimensions(A, 324, 324, 26730, 26730, 1.70664) ;

		AT = demo1_transpose(A) ;
		assert_dimensions(AT, 324, 324, 26730, 26730, 1.70664) ;

		D = demo1_multiply_add(A, AT) ;
		assert_dimensions(D, 324, 324, 65934, 65934, 5.42006) ;
	}

	public void test_demo1_t2() {
		DZcs T, A, AT, D ;

		InputStream in = get_stream(T2) ;

		T = demo1_load (in) ;
		assert_dimensions(T, 4, 4, 16, 10);

		A = demo1_compress(T) ;
		assert_dimensions(A, 4, 4, 10, 10, 106.07516) ;

		AT = demo1_transpose(A) ;
		assert_dimensions(AT, 4, 4, 10, 10, 144.29640) ;

		D = demo1_multiply_add(A, AT) ;
		assert_dimensions(D, 4, 4, 16, 16, 25308.3283) ;
	}

	public void test_demo1_t3() {
		DZcs T, A, AT, D ;

		InputStream in = get_stream(T3) ;

		T = demo1_load (in) ;
		assert_dimensions(T, 3, 4, 16, 12);

		A = demo1_compress(T) ;
		assert_dimensions(A, 3, 4, 12, 12, 3.05601) ;

		AT = demo1_transpose(A) ;
		assert_dimensions(AT, 4, 3, 12, 12, 3.93809) ;

		D = demo1_multiply_add(A, AT) ;
		assert_dimensions(D, 3, 3, 9, 9, 21.3318) ;
	}

	public void test_demo1_t4() {
		DZcs T, A, AT, D ;

		InputStream in = get_stream(T4) ;

		T = demo1_load (in) ;
		assert_dimensions(T, 2, 2, 4, 3);

		A = demo1_compress(T) ;
		assert_dimensions(A, 2, 2, 3, 3, 2.82843) ;

		AT = demo1_transpose(A) ;
		assert_dimensions(AT, 2, 2, 3, 3, 2.82843) ;

		D = demo1_multiply_add(A, AT) ;
		assert_dimensions(D, 2, 2, 4, 4, 12.0) ;
	}

	public void test_demo1_young1c() {
		DZcs T, A, AT, D ;

		InputStream in = get_stream(YOUNG1C) ;

		T = demo1_load (in) ;
		assert_dimensions(T, 841, 841, 4096, 4089);

		A = demo1_compress(T) ;
		assert_dimensions(A, 841, 841, 4089, 4089, 730.46) ;

		AT = demo1_transpose(A) ;
		assert_dimensions(AT, 841, 841, 4089, 4089, 730.46) ;

		D = demo1_multiply_add(A, AT) ;
		assert_dimensions(D, 841, 841, 10357, 10357, 1.0671436232e+06) ;
	}

}
