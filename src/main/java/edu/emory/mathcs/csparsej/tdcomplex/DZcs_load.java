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

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcs;

/**
 * Load a sparse matrix from a file.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class DZcs_load {

    /**
     * Loads a triplet matrix T from a file. Each line of the file contains
     * four values: a row index i, a column index j, a real value aij, and an
     * imaginary value aij.
     * The file is zero-based.
     *
     * @param fileName
     *            file name
     * @return T if successful, null on error
     */
    public static DZcs cs_load(String fileName) {
        int i, j;
        double[] x;
        DZcs T;
        BufferedReader in;
        try {
            in = new BufferedReader(new FileReader(fileName));
        } catch (FileNotFoundException e1) {
            return (null);
        }
        T = DZcs_util.cs_spalloc(0, 0, 1, true, true); /* allocate result */
        String line;
        try {
            while ((line = in.readLine()) != null) {
                String[] tokens = line.trim().split("\\s+");
                if (tokens.length != 4) {
                    return null;
                }
                i = Integer.parseInt(tokens[0]);
                j = Integer.parseInt(tokens[1]);
                double x_re = Double.parseDouble(tokens[2]);
                double x_im = Double.parseDouble(tokens[3]);
                x = new double[] {x_re, x_im};
                if (!DZcs_entry.cs_entry(T, i, j, x))
                    return (null);
            }
        } catch (IOException e) {
            return (null);
        }
        return (T);
    }
}
