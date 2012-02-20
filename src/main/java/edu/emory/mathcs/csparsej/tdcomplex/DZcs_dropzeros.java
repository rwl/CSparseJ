/*
 * CXSparse: a Concise Sparse matrix package.
 * Copyright (C) 2006-2011, Timothy A. Davis.
 * Copyright (C) 2011-2012, Richard W. Lincoln.
 * http://www.cise.ufl.edu/research/sparse/CXSparse
 *
 * -------------------------------------------------------------------------
 *
 * CXSparseJ is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * CXSparseJ is distributed in the hope that it will be useful,
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
 * Drop zeros from a sparse matrix.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class DZcs_dropzeros {

    private static class Cs_nonzero implements DZcs_ifkeep {
        @Override
        public boolean fkeep(int i, int j, double[] aij, Object other) {
            return (DZcs_complex.cs_creal(aij) != 0 && DZcs_complex.cs_cimag(aij) != 0);
        }
    }

    /**
     * Removes numerically zero entries from a matrix.
     *
     * @param A
     *            column-compressed matrix
     * @return nz, new number of entries in A, -1 on error
     */
    public static int cs_dropzeros(DZcs A) {
        return (DZcs_fkeep.cs_fkeep(A, new Cs_nonzero(), null)); /* keep all nonzero entries */
    }

}
