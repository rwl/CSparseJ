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
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex;

/**
 * Drop small entries from a sparse matrix.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class DZcs_droptol {

    private static class Cs_tol implements DZcs_ifkeep {
        @Override
        public boolean fkeep(int i, int j, double[] aij, Object other) {
            return (DZcs_complex.cs_cabs(aij) > (DZcs_complex.cs_cabs((double[]) other)));
        }
    }

    /**
     * Removes entries from a matrix with absolute value <= tol.
     *
     * @param A
     *            column-compressed matrix
     * @param tol
     *            drop tolerance
     * @return nz, new number of entries in A, -1 on error
     */
    public static int cs_droptol(DZcs A, double tol) {
        return (DZcs_fkeep.cs_fkeep(A, new Cs_tol(), tol)); /* keep all large entries */
    }
}
