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

package edu.emory.mathcs.csparsej.tdcomplex;

import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcsa;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcs;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcsn;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcss;

/**
 * Solves Ax=b where A is symmetric positive definite.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class DZcs_cholsol {

    /**
     * Solves Ax=b where A is symmetric positive definite; b is overwritten with
     * solution.
     *
     * @param order
     *            ordering method to use (0 or 1)
     * @param A
     *            column-compressed matrix, symmetric positive definite, only
     *            upper triangular part is used
     * @param b
     *            right hand side, b is overwritten with solution
     * @return true if successful, false on error
     */
    public static boolean cs_cholsol(int order, DZcs A, DZcsa b) {
        DZcsa x;
        DZcss S;
        DZcsn N;
        int n;
        boolean ok;
        if (!DZcs_util.CS_CSC(A) || b == null)
            return (false); /* check inputs */
        n = A.n;
        S = DZcs_schol.cs_schol(order, A); /* ordering and symbolic analysis */
        N = DZcs_chol.cs_chol(A, S); /* numeric Cholesky factorization */
        x = new DZcsa(n); /* get workspace */
        ok = (S != null && N != null);
        if (ok) {
            DZcs_ipvec.cs_ipvec(S.pinv, b, x, n); /* x = P*b */
            DZcs_lsolve.cs_lsolve(N.L, x); /* x = L\x */
            DZcs_ltsolve.cs_ltsolve(N.L, x); /* x = L'\x */
            DZcs_pvec.cs_pvec(S.pinv, x, b, n); /* b = P'*x */
        }
        return (ok);
    }

}
