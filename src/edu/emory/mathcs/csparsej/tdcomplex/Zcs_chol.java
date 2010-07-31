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

import edu.emory.mathcs.csparsej.tdcomplex.Zcs_common.Zcs;
import edu.emory.mathcs.csparsej.tdcomplex.Zcs_common.Zcsn;
import edu.emory.mathcs.csparsej.tdcomplex.Zcs_common.Zcss;
import org.apache.commons.math.complex.Complex;

/**
 * Sparse Cholesky.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class Zcs_chol {
    /**
     * Numeric Cholesky factorization LL=PAP'.
     *
     * @param A
     *            column-compressed matrix, only upper triangular part is used
     * @param S
     *            symbolic Cholesky analysis, pinv is optional
     * @return numeric Cholesky factorization, null on error
     */
    public static Zcsn cs_chol(Zcs A, Zcss S) {
        Complex d, lki, Lx[], x[], Cx[];
        int top, i, p, k, n, Li[], Lp[], cp[], pinv[], s[], c[], parent[], Cp[], Ci[];
        Zcs L, C;
        Zcsn N;
        if (!Zcs_util.CS_CSC(A) || S == null || S.cp == null || S.parent == null)
            return (null);
        n = A.n;
        N = new Zcsn(); /* allocate result */
        c = new int[2 * n]; /* get int workspace */
        x = new Complex[n]; /* get complex workspace */
        cp = S.cp;
        pinv = S.pinv;
        parent = S.parent;
        C = pinv != null ? Zcs_symperm.cs_symperm(A, pinv, true) : A;
        s = c;
        int s_offset = n;
        Cp = C.p;
        Ci = C.i;
        Cx = C.x;
        N.L = L = Zcs_util.cs_spalloc(n, n, cp[n], true, false); /* allocate result */
        Lp = L.p;
        Li = L.i;
        Lx = L.x;
        for (k = 0; k < n; k++)
            Lp[k] = c[k] = cp[k];
        for (k = 0; k < n; k++) /* compute L(k,:) for L*L' = C */
        {
            /* --- Nonzero pattern of L(k,:) ------------------------------------ */
            top = Zcs_ereach.cs_ereach(C, k, parent, s, s_offset, c); /* find pattern of L(k,:) */
            x[k] = Complex.ZERO; /* x (0:k) is now zero */
            for (p = Cp[k]; p < Cp[k + 1]; p++) /* x = full(triu(C(:,k))) */
            {
                if (Ci[p] <= k)
                    x[Ci[p]] = Cx[p];
            }
            d = x[k]; /* d = C(k,k) */
            x[k] = Complex.ZERO; /* clear x for k+1st iteration */
            /* --- Triangular solve --------------------------------------------- */
            for (; top < n; top++) /* solve L(0:k-1,0:k-1) * x = C(:,k) */
            {
                i = s[s_offset + top]; /* s [top..n-1] is pattern of L(k,:) */
                lki = x[i].divide(Lx[Lp[i]]); /* L(k,i) = x (i) / L(i,i) */
                x[i] = Complex.ZERO; /* clear x for k+1st iteration */
                for (p = Lp[i] + 1; p < c[i]; p++) {
                    x[Li[p]] = x[Li[p]].subtract(Lx[p].multiply(lki));
                }
                d = d.subtract(lki.multiply(lki)); /* d = d - L(k,i)*L(k,i) */
                p = c[i]++;
                Li[p] = k; /* store L(k,i) in column i */
                Lx[p] = lki;
            }
            /* --- Compute L(k,k) ----------------------------------------------- */
            if (d.getReal() <= 0 || d.getImaginary() != 0)
                return null; /* not pos def */
            p = c[k]++;
            Li[p] = k; /* store L(k,k) = sqrt (d) in column k */
            Lx[p] = d.sqrt();
        }
        Lp[n] = cp[n]; /* finalize L */
        return N;
    }
}
