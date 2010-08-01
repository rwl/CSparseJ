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

import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcs;

/**
 * Sparse matrix 1-norm.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class DZcs_norm {

    /**
     * Computes the 1-norm of a sparse matrix = max (sum (abs (A))), largest
     * column sum.
     *
     * @param A
     *            column-compressed matrix
     * @return the 1-norm if successful, -1 on error
     */
    public static double cs_norm(DZcs A) {
        int p, j, n, Ap[];
        Complex Ax[];
        double norm = 0, s;
        if (!DZcs_util.CS_CSC(A) || A.x == null)
            return (-1); /* check inputs */
        n = A.n;
        Ap = A.p;
        Ax = A.x;
        for (j = 0; j < n; j++) {
            for (s = 0, p = Ap[j]; p < Ap[j + 1]; p++)
                s += Ax[p].abs();
            norm = Math.max(norm, s);
        }
        return (norm);
    }
}
