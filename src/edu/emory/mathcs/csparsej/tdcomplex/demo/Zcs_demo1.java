/* ***** BEGIN LICENSE BLOCK *****
 *
 * CSparse: a Concise Sparse matrix package.
 * Copyright (c) 2006, Timothy A. Davis.
 * http://www.cise.ufl.edu/research/sparse/CSparse
 *
 * -------------------------------------------------------------------------
 *
 * CSparseJ is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * CSparseJ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this Module; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 *
 * ***** END LICENSE BLOCK ***** */

package edu.emory.mathcs.csparsej.tdcomplex.demo;

import org.apache.commons.math.complex.Complex;

import edu.emory.mathcs.csparsej.tdcomplex.Zcs_add;
import edu.emory.mathcs.csparsej.tdcomplex.Zcs_compress;
import edu.emory.mathcs.csparsej.tdcomplex.Zcs_entry;
import edu.emory.mathcs.csparsej.tdcomplex.Zcs_load;
import edu.emory.mathcs.csparsej.tdcomplex.Zcs_multiply;
import edu.emory.mathcs.csparsej.tdcomplex.Zcs_norm;
import edu.emory.mathcs.csparsej.tdcomplex.Zcs_print;
import edu.emory.mathcs.csparsej.tdcomplex.Zcs_transpose;
import edu.emory.mathcs.csparsej.tdcomplex.Zcs_util;
import edu.emory.mathcs.csparsej.tdcomplex.Zcs_common.Zcs;

/**
 * Read a matrix from a file and perform basic matrix operations.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class Zcs_demo1 {
    public static void main(String[] args) {
        Zcs T = null, A, Eye, AT, C, D;
        int i, m;
        if (args.length == 0) {
            throw new IllegalArgumentException("Usage: java edu.emory.mathcs.csparsej.tdcomplex.demo.Zcs_demo1 fileName");
        }
        T = Zcs_load.cs_load(args[0]); /* load triplet matrix T from file */
        System.out.print("T:\n");
        Zcs_print.cs_print(T, false); /* print T */
        A = Zcs_compress.cs_compress(T); /* A = compressed-column form of T */
        System.out.print("A:\n");
        Zcs_print.cs_print(A, false); /* print A */
        T = null; /* clear T */
        AT = Zcs_transpose.cs_transpose(A, true); /* AT = A' */
        System.out.print("AT:\n");
        Zcs_print.cs_print(AT, false); /* print AT */
        m = A != null ? A.m : 0; /* m = # of rows of A */
        T = Zcs_util.cs_spalloc(m, m, m, true, true); /* create triplet identity matrix */
        for (i = 0; i < m; i++)
            Zcs_entry.cs_entry(T, i, i, Complex.ONE);
        Eye = Zcs_compress.cs_compress(T); /* Eye = speye (m) */
        T = null;
        C = Zcs_multiply.cs_multiply(A, AT); /* C = A*A' */
        D = Zcs_add.cs_add(C, Eye, Complex.ONE, new Complex(Zcs_norm.cs_norm(C), 0.0)); /* D = C + Eye*norm (C,1) */
        System.out.print("D:\n");
        Zcs_print.cs_print(D, false); /* print D */
    }
}
