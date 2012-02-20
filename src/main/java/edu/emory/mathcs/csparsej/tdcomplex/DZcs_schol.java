/*
 * CSparse: a Concise Sparse matrix package.
 * Copyright (C) 2006-2011, Timothy A. Davis.
 * Copyright (C) 2011-2012, Richard W. Lincoln.
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
 *
 */

package edu.emory.mathcs.csparsej.tdcomplex;

import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcs;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcss;

/**
 * Symbolic Cholesky ordering and analysis.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class DZcs_schol {
    /**
     * Ordering and symbolic analysis for a Cholesky factorization.
     *
     * @param order
     *            ordering option (0 or 1)
     * @param A
     *            column-compressed matrix
     * @return symbolic analysis for Cholesky, null on error
     */
    public static DZcss cs_schol(int order, DZcs A) {
        int n, c[], post[], P[];
        DZcs C;
        DZcss S;
        if (!DZcs_util.CS_CSC(A))
            return (null); /* check inputs */
        n = A.n;
        S = new DZcss(); /* allocate result S */
        P = DZcs_amd.cs_amd(order, A); /* P = amd(A+A'), or natural */
        S.pinv = DZcs_pinv.cs_pinv(P, n); /* find inverse permutation */
        if (order != 0 && S.pinv == null)
            return null;
        C = DZcs_symperm.cs_symperm(A, S.pinv, false); /* C = spones(triu(A(P,P))) */
        S.parent = DZcs_etree.cs_etree(C, false); /* find etree of C */
        post = DZcs_post.cs_post(S.parent, n); /* postorder the etree */
        c = DZcs_counts.cs_counts(C, S.parent, post, false); /* find column counts of chol(C) */
        S.cp = new int[n + 1]; /* allocate result S.cp */
        S.unz = S.lnz = DZcs_cumsum.cs_cumsum(S.cp, c, n); /* find column pointers for L */
        return ((S.lnz >= 0) ? S : null);
    }
}
