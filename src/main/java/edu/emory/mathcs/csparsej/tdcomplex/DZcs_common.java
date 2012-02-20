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

/**
 * Common data structures.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class DZcs_common {

    public static final int CS_VER = 1; /* CSparseJ Version 1.0.0 */
    public static final int CS_SUBVER = 0;
    public static final int CS_SUBSUB = 0;
    public static final String CS_DATE = "June 13, 2009"; /* CSparseJ release date */
    public static final String CS_COPYRIGHT = "Copyright (C) Timothy A. Davis, 2006-2009";

    /**
     *
     * Complex array.
     *
     */
    public static class DZcsa {

        /**
         * numerical values
         */
        public double[] x;

        /**
         *
         */
        public DZcsa() {

        }

        /**
         * Constructs an array of the given length.
         */
        public DZcsa(int len) {
            this.x = new double[2 * len];
        }

        /**
         *
         * @param idx
         * @return
         */
        public double[] get(int idx) {
            return new double[] {x[idx], x[idx + 1]};
        }

        /**
         *
         * @param idx
         * @param val
         */
        public void set(int idx, double[] val) {
            x[idx] = val[0];
            x[idx + 1] = val[1];
        }
    }

    /**
     *
     * Complex matrix in compressed-column or triplet form.
     *
     */
    public static class DZcs {

        /**
         * maximum number of entries
         */
        public int nzmax;

        /**
         * number of rows
         */
        public int m;

        /**
         * number of columns
         */
        public int n;

        /**
         * column pointers (size n+1) or col indices (size nzmax)
         */
        public int[] p;

        /**
         * row indices, size nzmax
         */
        public int[] i;

        /**
         * numerical values, size 2 * nzmax
         */
        public double[] x;

        /**
         * # of entries in triplet matrix, -1 for compressed-col
         */
        public int nz;

        public DZcs() {

        }

        public double[] get(int idx) {
            return new double[] {x[idx], x[idx + 1]};
        }

        public void set(int idx, double[] val) {
            x[idx] = val[0];
            x[idx + 1] = val[1];
        }

//        public synchronized double[] getQuick(int row, int column) {
//            //        int k = cern.colt.Sorting.binarySearchFromTo(dcs.i, row, dcs.p[column], dcs.p[column + 1] - 1);
//            int k = searchFromTo(rowIndexes, row, columnPointers[column], columnPointers[column + 1] - 1);
//            double[] v = new double[2];
//            if (k >= 0) {
//                v[0] = values[2 * k];
//                v[1] = values[2 * k + 1];
//            }
//            return v;
//        }
//
//        public synchronized void setQuick(int row, int column, double[] value) {
//            //        int k = cern.colt.Sorting.binarySearchFromTo(dcs.i, row, dcs.p[column], dcs.p[column + 1] - 1);
//            int k = searchFromTo(rowIndexes, row, columnPointers[column], columnPointers[column + 1] - 1);
//
//            if (k >= 0) { // found
//                if (value[0] == 0 && value[1] == 0) {
//                    remove(column, k);
//                } else {
//                    values[2 * k] = value[0];
//                    values[2 * k + 1] = value[1];
//                }
//                return;
//            }
//
//            if (value[0] != 0 || value[1] != 0) {
//                k = -k - 1;
//                insert(row, column, k, value);
//            }
//        }
//
//        protected void insert(int row, int column, int index, double[] value) {
//            IntArrayList rowIndexesList = new IntArrayList(rowIndexes);
//            rowIndexesList.setSizeRaw(columnPointers[columns]);
//            DoubleArrayList valuesList = new DoubleArrayList(values);
//            valuesList.setSizeRaw(2 * columnPointers[columns]);
//            rowIndexesList.beforeInsert(index, row);
//            valuesList.beforeInsert(2 * index, value[0]);
//            valuesList.beforeInsert(2 * index + 1, value[1]);
//            for (int i = columnPointers.length; --i > column;)
//                columnPointers[i]++;
//            rowIndexes = rowIndexesList.elements();
//            values = valuesList.elements();
//        }
//
//        private static int searchFromTo(int[] list, int key, int from, int to) {
//            while (from <= to) {
//                if (list[from] == key) {
//                    return from;
//                } else {
//                    from++;
//                    continue;
//                }
//            }
//            return -(from + 1); // key not found.
//        }

    };

    /**
     *
     * Output of symbolic Cholesky, LU, or QR analysis.
     *
     */
    public static class DZcss {
        /**
         * inverse row perm. for QR, fill red. perm for Chol
         */
        public int[] pinv;

        /**
         * fill-reducing column permutation for LU and QR
         */
        public int[] q;

        /**
         * elimination tree for Cholesky and QR
         */
        public int[] parent;

        /**
         * column pointers for Cholesky, row counts for QR
         */
        public int[] cp;

        /**
         * leftmost[i] = min(find(A(i,:))), for QR
         */
        public int[] leftmost;

        /**
         * # of rows for QR, after adding fictitious rows
         */
        public int m2;

        /**
         * # entries in L for LU or Cholesky; in V for QR
         */
        public int lnz;

        /**
         * # entries in U for LU; in R for QR
         */
        public int unz;

        public DZcss() {
        }
    };

    /**
     *
     * Output of numeric Cholesky, LU, or QR factorization
     *
     */
    public static class DZcsn {
        /**
         * L for LU and Cholesky, V for QR
         */
        public DZcs L;

        /**
         * U for LU, R for QR, not used for Cholesky
         */
        public DZcs U;

        /**
         * partial pivoting for LU
         */
        public int[] pinv;

        /**
         * beta [0..n-1] for QR
         */
        public double[] B;

        public DZcsn() {
        }

    };

    /**
     *
     * Output of Dulmage-Mendelsohn decomposition.
     *
     */
    public static class DZcsd {

        /**
         * size m, row permutation
         */
        public int[] p;

        /**
         * size n, column permutation
         */
        public int[] q;

        /**
         * size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q)
         */
        public int[] r;

        /**
         * size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q)
         */
        public int[] s;

        /**
         * # of blocks in fine dmperm decomposition
         */
        public int nb;

        /**
         * coarse row decomposition
         */
        public int[] rr;

        /**
         * coarse column decomposition
         */
        public int[] cc;
    };
}
