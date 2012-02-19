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
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex;

/**
 * Sparse Cholesky.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class DZcs_chol {
    /**
     * Numeric Cholesky factorization LL=PAP'.
     *
     * @param A
     *            column-compressed matrix, only upper triangular part is used
     * @param S
     *            symbolic Cholesky analysis, pinv is optional
     * @return numeric Cholesky factorization, null on error
     */
    public static DZcsn cs_chol(DZcs A, DZcss S) {
        double[] d, lki;
        DZcsa Lx = new DZcsa(), x, Cx = new DZcsa();
        int top, i, p, k, n, Li[], Lp[], cp[], pinv[], s[], c[], parent[], Cp[], Ci[];
        DZcs L, C;
        DZcsn N;
        if (!DZcs_util.CS_CSC(A) || S == null || S.cp == null || S.parent == null)
            return (null);
        n = A.n;
        N = new DZcsn(); /* allocate result */
        c = new int[2 * n]; /* get int workspace */
        x = new DZcsa(n); /* get complex workspace */
        cp = S.cp;
        pinv = S.pinv;
        parent = S.parent;
        C = pinv != null ? DZcs_symperm.cs_symperm(A, pinv, true) : A;
        s = c;
        int s_offset = n;
        Cp = C.p;
        Ci = C.i;
        Cx.x = C.x;
        N.L = L = DZcs_util.cs_spalloc(n, n, cp[n], true, false); /* allocate result */
        Lp = L.p;
        Li = L.i;
        Lx.x = L.x;
        for (k = 0; k < n; k++)
            Lp[k] = c[k] = cp[k];
        for (k = 0; k < n; k++) /* compute L(k,:) for L*L' = C */
        {
            /* --- Nonzero pattern of L(k,:) ------------------------------------ */
            top = DZcs_ereach.cs_ereach(C, k, parent, s, s_offset, c); /* find pattern of L(k,:) */
            x.set(k, DZcs_complex.cs_czero()); /* x (0:k) is now zero */
            for (p = Cp[k]; p < Cp[k + 1]; p++) /* x = full(triu(C(:,k))) */
            {
                if (Ci[p] <= k)
                    x.set(Ci[p], Cx.get(p));
            }
            d = x.get(k); /* d = C(k,k) */
            x.set(k, DZcs_complex.cs_czero()); /* clear x for k+1st iteration */
            /* --- Triangular solve --------------------------------------------- */
            for (; top < n; top++) /* solve L(0:k-1,0:k-1) * x = C(:,k) */
            {
                i = s[s_offset + top]; /* s [top..n-1] is pattern of L(k,:) */
                lki = DZcs_complex.cs_cdiv(x.get(i), Lx.get(Lp[i])); /* L(k,i) = x (i) / L(i,i) */
                x.set(i, DZcs_complex.cs_czero()); /* clear x for k+1st iteration */
                for (p = Lp[i] + 1; p < c[i]; p++) {
                    x.set(Li[p], DZcs_complex.cs_cminus(x.get(Li[p]), DZcs_complex.cs_cmult(lki, Lx.get(p))));
                }
                d = DZcs_complex.cs_cminus(d, DZcs_complex.cs_csquare(lki)); /* d = d - L(k,i)*L(k,i) */
                p = c[i]++;
                Li[p] = k; /* store L(k,i) in column i */
                Lx.set(p, lki);
            }
            /* --- Compute L(k,k) ----------------------------------------------- */
            if (d[0] <= 0 || d[1] != 0)
                return null; /* not pos def */
            p = c[k]++;
            Li[p] = k; /* store L(k,k) = sqrt (d) in column k */
            Lx.set(p, DZcs_complex.cs_csqrt(d));
        }
        Lp[n] = cp[n]; /* finalize L */
        return N;
    }
}
