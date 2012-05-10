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

package edu.emory.mathcs.csparsej.tdouble.test;

import java.io.InputStream;

import static edu.emory.mathcs.csparsej.tdouble.Dcs_add.cs_add ;
import static edu.emory.mathcs.csparsej.tdouble.Dcs_compress.cs_compress ;
import static edu.emory.mathcs.csparsej.tdouble.Dcs_entry.cs_entry ;
import static edu.emory.mathcs.csparsej.tdouble.Dcs_load.cs_load ;
import static edu.emory.mathcs.csparsej.tdouble.Dcs_multiply.cs_multiply ;
import static edu.emory.mathcs.csparsej.tdouble.Dcs_norm.cs_norm ;
import static edu.emory.mathcs.csparsej.tdouble.Dcs_transpose.cs_transpose ;
import static edu.emory.mathcs.csparsej.tdouble.Dcs_util.cs_spalloc ;
//import static edu.emory.mathcs.csparsej.tdouble.Dcs_print.cs_print ;

import edu.emory.mathcs.csparsej.tdouble.Dcs_common.Dcs ;

/**
 * Read a matrix from a file and perform basic matrix operations.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class Dcs_test1 extends Dcs_test {

	protected static Dcs demo1_load(InputStream in) {
		Dcs T = cs_load (in) ;			/* load triplet matrix T from file */
		//System.out.print("T:\n") ;
		//cs_print (T, false) ;			/* print T */
		return T ;
	}

	protected static Dcs demo1_compress(Dcs T) {
		Dcs A = cs_compress (T) ;		/* A = compressed-column form of T */
		//System.out.print("A:\n") ;
		//cs_print (A, false) ;			/* print A */
		return A ;
	}

	protected static Dcs demo1_transpose(Dcs A) {
		Dcs AT = cs_transpose (A, true) ;	/* AT = A' */
		//System.out.print("AT:\n") ;
		//cs_print (AT, false) ;			/* print AT */
		return AT ;
	}

	protected static Dcs demo1_multiply_add(Dcs A, Dcs AT) {
		Dcs T, Eye, C, D ;
		int m = A != null ? A.m : 0 ;		/* m = # of rows of A */
		T = cs_spalloc (m, m, m, true, true);	/* create triplet identity matrix */
		for (int i = 0; i < m; i++) cs_entry (T, i, i, 1) ;
		Eye = cs_compress (T) ;			/* Eye = speye (m) */
		T = null ;
		C = cs_multiply (A, AT) ;		/* C = A*A' */
		D = cs_add(C, Eye, 1, cs_norm (C)) ;	/* D = C + Eye*norm (C,1) */
		//System.out.print("D:\n") ;
		//cs_print(D, false) ;			/* print D */
		return D ;
	}

	public void test_demo1_ash219() {
		Dcs T, A, AT, D ;

		InputStream in = get_stream(ASH219) ;

		T = demo1_load (in) ;
		assert_dimensions(T, 219, 85, 512, 438);

		A = demo1_compress(T) ;
		assert_dimensions(A, 219, 85, 438, 438, 9) ;

		AT = demo1_transpose(A) ;
		assert_dimensions(AT, 85, 219, 438, 438, 2) ;

		D = demo1_multiply_add(A, AT) ;
		assert_dimensions(D, 219, 219, 2205, 2205, 32) ;
	}

	public void test_demo1_bcsstk01() {
		Dcs T, A, AT, D ;

		InputStream in = get_stream(BCSSTK01) ;

		T = demo1_load (in) ;
		assert_dimensions(T, 48, 48, 256, 224);

		A = demo1_compress(T) ;
		assert_dimensions(A, 48, 48, 224, 224, 3.00944e+09, 1e4) ;

		AT = demo1_transpose(A) ;
		assert_dimensions(AT, 48, 48, 224, 224, 3.57095e+09, 1e4) ;

		D = demo1_multiply_add(A, AT) ;
		assert_dimensions(D, 48, 48, 764, 764, 1.73403e+19, 1e14) ;
	}

	public void test_demo1_bcsstk16() {
		Dcs T, A, AT, D ;

		InputStream in = get_stream(BCSSTK16) ;

		T = demo1_load (in) ;
		assert_dimensions(T, 4884, 4884, 262144, 147631);

		A = demo1_compress(T) ;
		assert_dimensions(A, 4884, 4884, 147631, 147631, 4.91422e+09, 1e4) ;

		AT = demo1_transpose(A) ;
		assert_dimensions(AT, 4884, 4884, 147631, 147631, 5.47522e+09, 1e4) ;

		D = demo1_multiply_add(A, AT) ;
		assert_dimensions(D, 4884, 4884, 544856, 544856, 4.13336e+19, 1e14) ;
	}

	public void test_demo1_fs_183_1() {
		Dcs T, A, AT, D ;

		InputStream in = get_stream(FS_183_1) ;

		T = demo1_load (in) ;
		assert_dimensions(T, 183, 183, 2048, 1069);

		A = demo1_compress(T) ;
		assert_dimensions(A, 183, 183, 1069, 1069, 1.70318e+09, 1e4) ;

		AT = demo1_transpose(A) ;
		assert_dimensions(AT, 183, 183, 1069, 1069, 8.22724e+08, 1e3) ;

		D = demo1_multiply_add(A, AT) ;
		assert_dimensions(D, 183, 183, 19665, 19665, 2.80249e+18, 1e13) ;
	}

	public void test_demo1_ibm32a() {
		Dcs T, A, AT, D ;

		InputStream in = get_stream(IBM32A) ;

		T = demo1_load (in) ;
		assert_dimensions(T, 32, 31, 128, 123);

		A = demo1_compress(T) ;
		assert_dimensions(A, 32, 31, 123, 123, 7) ;

		AT = demo1_transpose(A) ;
		assert_dimensions(AT, 31, 32, 123, 123, 8) ;

		D = demo1_multiply_add(A, AT) ;
		assert_dimensions(D, 32, 32, 386, 386, 70) ;
	}

	public void test_demo1_ibm32b() {
		Dcs T, A, AT, D ;

		InputStream in = get_stream(IBM32B) ;

		T = demo1_load (in) ;
		assert_dimensions(T, 31, 32, 128, 123);

		A = demo1_compress(T) ;
		assert_dimensions(A, 31, 32, 123, 123, 8) ;

		AT = demo1_transpose(A) ;
		assert_dimensions(AT, 32, 31, 123, 123, 7) ;

		D = demo1_multiply_add(A, AT) ;
		assert_dimensions(D, 31, 31, 373, 373, 64) ;
	}

	public void test_demo1_lp_afiro() {
		Dcs T, A, AT, D ;

		InputStream in = get_stream(LP_AFIRO) ;

		T = demo1_load (in) ;
		assert_dimensions(T, 27, 51, 128, 102);

		A = demo1_compress(T) ;
		assert_dimensions(A, 27, 51, 102, 102, 3.429) ;

		AT = demo1_transpose(A) ;
		assert_dimensions(AT, 51, 27, 102, 102, 20.525) ;

		D = demo1_multiply_add(A, AT) ;
		assert_dimensions(D, 27, 27, 153, 153, 128.963) ;
	}

	public void test_demo1_mbeacxc() {
		Dcs T, A, AT, D ;

		InputStream in = get_stream(MBEACXC) ;

		T = demo1_load (in) ;
		assert_dimensions(T, 492, 490, 65536, 49920);

		A = demo1_compress(T) ;
		assert_dimensions(A, 492, 490, 49920, 49920, 0.928629) ;

		AT = demo1_transpose(A) ;
		assert_dimensions(AT, 490, 492, 49920, 49920, 16.5516) ;

		D = demo1_multiply_add(A, AT) ;
		assert_dimensions(D, 492, 492, 157350, 157350, 19.6068) ;
	}

	public void test_demo1_t1() {
		Dcs T, A, AT, D ;

		InputStream in = get_stream(T1) ;

		T = demo1_load (in) ;
		assert_dimensions(T, 4, 4, 16, 10);

		A = demo1_compress(T) ;
		assert_dimensions(A, 4, 4, 10, 10, 11.1) ;

		AT = demo1_transpose(A) ;
		assert_dimensions(AT, 4, 4, 10, 10, 7.7) ;

		D = demo1_multiply_add(A, AT) ;
		assert_dimensions(D, 4, 4, 16, 16, 139.58) ;
	}

	public void test_demo1_west0067() {
		Dcs T, A, AT, D ;

		InputStream in = get_stream(WEST0067) ;

		T = demo1_load (in) ;
		assert_dimensions(T, 67, 67, 512, 299);

		A = demo1_compress(T) ;
		assert_dimensions(A, 67, 67, 299, 299, 6.14337) ;

		AT = demo1_transpose(A) ;
		assert_dimensions(AT, 67, 67, 299, 299, 6.59006) ;

		D = demo1_multiply_add(A, AT) ;
		assert_dimensions(D, 67, 67, 1041, 1041, 61.0906) ;
	}

}
