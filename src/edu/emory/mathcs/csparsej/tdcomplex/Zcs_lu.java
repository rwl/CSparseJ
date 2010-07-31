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

import org.apache.commons.math.complex.Complex;

import edu.emory.mathcs.csparsej.tdcomplex.Zcs_common.Zcs;
import edu.emory.mathcs.csparsej.tdcomplex.Zcs_common.Zcsn;
import edu.emory.mathcs.csparsej.tdcomplex.Zcs_common.Zcss;

/**
 * Sparse LU factorization.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class Zcs_lu {

    /**
     * Sparse LU factorization of a square matrix, PAQ = LU.
     *
     * @param A
     *            column-compressed matrix
     * @param S
     *            symbolic LU analysis
     * @param tol
     *            partial pivoting threshold (1 for partial pivoting)
     * @return numeric LU factorization, null on error
     */
    public static Zcsn cs_lu(Zcs A, Zcss S, double tol) {
        Zcs L, U;
        Zcsn N;
        Complex pivot, Lx[], Ux[], x[];
        double a, t;
        int Lp[], Li[], Up[], Ui[], pinv[], xi[], q[], n, ipiv, k, top, p, i, col, lnz, unz;
        if (!Zcs_util.CS_CSC(A) || S == null)
            return (null); /* check inputs */
        n = A.n;
        q = S.q;
        lnz = S.lnz;
        unz = S.unz;
        x = new Complex[n]; /* get double workspace */
        xi = new int[2 * n]; /* get int workspace */
        N = new Zcsn(); /* allocate result */
        N.L = L = Zcs_util.cs_spalloc(n, n, lnz, true, false); /* allocate result L */
        N.U = U = Zcs_util.cs_spalloc(n, n, unz, true, false); /* allocate result U */
        N.pinv = pinv = new int[n]; /* allocate result pinv */
        Lp = L.p;
        Up = U.p;
        for (i = 0; i < n; i++)
            x[i] = Complex.ZERO; /* clear workspace */
        for (i = 0; i < n; i++)
            pinv[i] = -1; /* no rows pivotal yet */
        for (k = 0; k <= n; k++)
            Lp[k] = 0; /* no cols of L yet */
        lnz = unz = 0;
        for (k = 0; k < n; k++) /* compute L(:,k) and U(:,k) */
        {
            /* --- Triangular solve --------------------------------------------- */
            Lp[k] = lnz; /* L(:,k) starts here */
            Up[k] = unz; /* U(:,k) starts here */
            if (lnz + n > L.nzmax) {
                Zcs_util.cs_sprealloc(L, 2 * L.nzmax + n);
            }
            if (unz + n > U.nzmax) {
                Zcs_util.cs_sprealloc(U, 2 * U.nzmax + n);
            }
            Li = L.i;
            Lx = L.x;
            Ui = U.i;
            Ux = U.x;
            col = q != null ? (q[k]) : k;
            top = Zcs_spsolve.cs_spsolve(L, A, col, xi, x, pinv, true); /* x = L\A(:,col) */
            /* --- Find pivot --------------------------------------------------- */
            ipiv = -1;
            a = -1;
            for (p = top; p < n; p++) {
                i = xi[p]; /* x(i) is nonzero */
                if (pinv[i] < 0) /* row i is not yet pivotal */
                {
                    if ((t = x[i].abs()) > a) {
                        a = t; /* largest pivot candidate so far */
                        ipiv = i;
                    }
                } else /* x(i) is the entry U(pinv[i],k) */
                {
                    Ui[unz] = pinv[i];
                    Ux[unz++] = x[i];
                }
            }
            if (ipiv == -1 || a <= 0)
                return (null);
            if (pinv[col] < 0 && x[col].abs() >= a * tol)
                ipiv = col;
            /* --- Divide by pivot ---------------------------------------------- */
            pivot = x[ipiv]; /* the chosen pivot */
            Ui[unz] = k; /* last entry in U(:,k) is U(k,k) */
            Ux[unz++] = pivot;
            pinv[ipiv] = k; /* ipiv is the kth pivot row */
            Li[lnz] = ipiv; /* first entry in L(:,k) is L(k,k) = 1 */
            Lx[lnz++] = Complex.ONE;
            for (p = top; p < n; p++) /* L(k+1:n,k) = x / pivot */
            {
                i = xi[p];
                if (pinv[i] < 0) /* x(i) is an entry in L(:,k) */
                {
                    Li[lnz] = i; /* save unpermuted row in L */
                    Lx[lnz++] = x[i].divide(pivot); /* scale pivot column */
                }
                x[i] = Complex.ZERO; /* x [0..n-1] = 0 for next k */
            }
        }
        /* --- Finalize L and U ------------------------------------------------- */
        Lp[n] = lnz;
        Up[n] = unz;
        Li = L.i; /* fix row indices of L for final pinv */
        for (p = 0; p < lnz; p++)
            Li[p] = pinv[Li[p]];
        Zcs_util.cs_sprealloc(L, 0); /* remove extra space from L and U */
        Zcs_util.cs_sprealloc(U, 0);
        return N;
    }

}
