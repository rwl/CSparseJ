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

import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcsa;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex;

/**
 * Compute Householder reflection.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class DZcs_house {

    /**
     * Compute a Householder reflection [v,beta,s]=house(x), overwrite x with v,
     * where (I-beta*v*v')*x = s*e1 and e1 = [1 0 ... 0]'.
     * Note that this CXSparseJ version is different than CXSparseJ.  See Higham,
     * Accuracy & Stability of Num Algorithms, 2nd ed, 2002, page 357.
     *
     * @param x
     *            x on output, v on input
     * @param beta
     *            scalar beta
     * @param n
     *            the length of x
     * @return norm2(x), -1 on error
     */
    public static double[] cs_house(DZcsa x, double[] beta, int n) {
        double[] s = DZcs_complex.cs_czero();
        int i;
        if (x == null)
            return new double[] {-1.0, 0.0}; /* check inputs */
        for (i = 1; i < n; i++)
            s = DZcs_complex.cs_cplus(s, DZcs_complex.cs_cmult(x.get(i), DZcs_complex.cs_conj(x.get(i))));
        s = DZcs_complex.cs_csqrt(s);
        if (DZcs_complex.cs_cequal(s, DZcs_complex.cs_czero())) {
            beta[0] = 0.0;
            x.set(0, DZcs_complex.cs_cone());
        } else {
            /* s = sign(x[0]) * norm (x) ; */
            if (!DZcs_complex.cs_cequal(x.get(0), DZcs_complex.cs_czero())) {
                s = DZcs_complex.cs_cmult(s, DZcs_complex.cs_cdiv(x.get(0), new double[] {DZcs_complex.cs_cabs(x.get(0)), 0.0}));
            }
            x.set(0, DZcs_complex.cs_cplus(x.get(0), s));
            beta[0] = 1 / DZcs_complex.cs_creal(DZcs_complex.cs_cmult(DZcs_complex.cs_conj(s), x.get(0)));
        }
        return (s);
    }
}
