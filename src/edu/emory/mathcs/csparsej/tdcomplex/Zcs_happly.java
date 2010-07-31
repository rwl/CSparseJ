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

/**
 * Apply Householder reflection.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class Zcs_happly {

    /**
     * Applies a Householder reflection to a dense vector, x = (I -
     * beta*v*v')*x.
     *
     * @param V
     *            column-compressed matrix of Householder vectors
     * @param i
     *            v = V(:,i), the ith column of V
     * @param beta
     *            scalar beta
     * @param x
     *            vector x of size m
     * @return true if successful, false on error
     */
    public static boolean cs_happly(Zcs V, int i, Complex beta, Complex[] x) {
        int p, Vp[], Vi[];
        Complex Vx[], tau = Complex.ZERO;
        if (!Zcs_util.CS_CSC(V) || x == null)
            return (false); /* check inputs */
        Vp = V.p;
        Vi = V.i;
        Vx = V.x;
        for (p = Vp[i]; p < Vp[i + 1]; p++) /* tau = v'*x */
        {
            tau = tau.add(Vx[p].multiply(x[Vi[p]]));
        }
        tau = tau.multiply(beta); /* tau = beta*(v'*x) */
        for (p = Vp[i]; p < Vp[i + 1]; p++) /* x = x - v*tau */
        {
            x[Vi[p]] = x[Vi[p]].subtract(Vx[p].multiply(tau));
        }
        return (true);
    }

}
