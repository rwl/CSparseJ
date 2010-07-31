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

/**
 * Compute Householder reflection.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class Zcs_house {

    /**
     * Compute a Householder reflection [v,beta,s]=house(x), overwrite x with v,
     * where (I-beta*v*v')*x = s*e1 and e1 = [1 0 ... 0]'.
     * Note that this CXSparseJ version is different than CSparseJ.  See Higham,
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
    public static Complex cs_house(Complex[] x, double[] beta, int n) {
        Complex s = Complex.ZERO;
        int i;
        if (x == null)
            return Complex.ONE.multiply(-1.0); /* check inputs */
        for (i = 1; i < n; i++)
            s = s.add(x[i].multiply(x[i].conjugate()));
        s = s.sqrt();
        if (s.equals(Complex.ZERO)) {
            beta[0] = 0.0;
            x[0] = Complex.ONE;
        } else {
            /* s = sign(x[0]) * norm (x) ; */
            if (!x[0].equals(Complex.ZERO)) {
                s = s.multiply(x[0].divide(new Complex(x[0].abs(), 0.0)));
            }
            x[0] = x[0].add(s);
            beta[0] = 1 / s.conjugate().multiply(x[0]).getReal();
        }
        return (s);
    }
}
