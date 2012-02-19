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
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex;

/**
 * Sparse lower or upper triangular solve. x=G\b where G, x, and b are sparse,
 * and G upper/lower triangular.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class DZcs_spsolve {
    /**
     * solve Gx=b(:,k), where G is either upper (lo=false) or lower (lo=true)
     * triangular.
     *
     * @param G
     *            lower or upper triangular matrix in column-compressed form
     * @param B
     *            right hand side, b=B(:,k)
     * @param k
     *            use kth column of B as right hand side
     * @param xi
     *            size 2*n, nonzero pattern of x in xi[top..n-1]
     * @param x
     *            size n, x in x[xi[top..n-1]]
     * @param pinv
     *            mapping of rows to columns of G, ignored if null
     * @param lo
     *            true if lower triangular, false if upper
     * @return top, -1 in error
     */
    public static int cs_spsolve(DZcs G, DZcs B, int k, int[] xi, DZcsa x, int[] pinv, boolean lo) {
        int j, J, p, q, px, top, n, Gp[], Gi[], Bp[], Bi[];
        DZcsa Gx = new DZcsa(), Bx = new DZcsa();
        if (!DZcs_util.CS_CSC(G) || !DZcs_util.CS_CSC(B) || xi == null || x == null)
            return (-1);
        Gp = G.p;
        Gi = G.i;
        Gx.x = G.x;
        n = G.n;
        Bp = B.p;
        Bi = B.i;
        Bx.x = B.x;
        top = DZcs_reach.cs_reach(G, B, k, xi, pinv); /* xi[top..n-1]=Reach(B(:,k)) */
        for (p = top; p < n; p++)
            x.set(xi[p], DZcs_complex.cs_czero()); /* clear x */
        for (p = Bp[k]; p < Bp[k + 1]; p++)
            x.set(Bi[p], Bx.get(p)); /* scatter B */
        for (px = top; px < n; px++) {
            j = xi[px]; /* x(j) is nonzero */
            J = pinv != null ? (pinv[j]) : j; /* j maps to col J of G */
            if (J < 0)
                continue; /* column J is empty */
            x.set(j, DZcs_complex.cs_cdiv(x.get(j), Gx.get(lo ? (Gp[J]) : (Gp[J + 1] - 1))));/* x(j) /= G(j,j) */
            p = lo ? (Gp[J] + 1) : (Gp[J]); /* lo: L(j,j) 1st entry */
            q = lo ? (Gp[J + 1]) : (Gp[J + 1] - 1); /* up: U(j,j) last entry */
            for (; p < q; p++) {
                x.set(Gi[p], DZcs_complex.cs_cminus(x.get(Gi[p]), DZcs_complex.cs_cmult(Gx.get(p), x.get(j)))); /* x(i) -= G(i,j) * x(j) */
            }
        }
        return (top); /* return top of stack */
    }
}
