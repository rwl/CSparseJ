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
 * Convert a triplet form to compressed-column form.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class Zcs_compress {

    /**
     * C = compressed-column form of a triplet matrix T. The columns of C are
     * not sorted, and duplicate entries may be present in C.
     *
     * @param T
     *            triplet matrix
     * @return C if successful, null on error
     */
    public static Zcs cs_compress(Zcs T) {
        int m, n, nz, p, k, Cp[], Ci[], w[], Ti[], Tj[];
        Complex Cx[], Tx[];
        Zcs C;
        if (!Zcs_util.CS_TRIPLET(T))
            return (null); /* check inputs */
        m = T.m;
        n = T.n;
        Ti = T.i;
        Tj = T.p;
        Tx = T.x;
        nz = T.nz;
        C = Zcs_util.cs_spalloc(m, n, nz, Tx != null, false); /* allocate result */
        w = new int[n]; /* get workspace */
        Cp = C.p;
        Ci = C.i;
        Cx = C.x;
        for (k = 0; k < nz; k++)
            w[Tj[k]]++; /* column counts */
        Zcs_cumsum.cs_cumsum(Cp, w, n); /* column pointers */
        for (k = 0; k < nz; k++) {
            Ci[p = w[Tj[k]]++] = Ti[k]; /* A(i,j) is the pth entry in C */
            if (Cx != null)
                Cx[p] = Tx[k];
        }
        return C;
    }
}
