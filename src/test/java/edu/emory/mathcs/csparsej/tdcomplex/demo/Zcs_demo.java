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

package edu.emory.mathcs.csparsej.tdcomplex.demo;

import java.util.Random;

import org.apache.commons.math.complex.Complex;

import edu.emory.mathcs.csparsej.tdcomplex.DZcs_add;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_chol;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_cholsol;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_compress;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_dmperm;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_droptol;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_dropzeros;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_dupl;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_fkeep;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_gaxpy;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_ifkeep;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_ipvec;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_load;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_lsolve;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_ltsolve;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_lusol;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_multiply;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_norm;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_permute;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_pinv;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_pvec;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_qrsol;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_schol;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_transpose;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_updown;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_util;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcs;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcsd;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcsn;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcss;

/**
 * Support routines for Zcs_demo*.java
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class Zcs_demo {

    /**
     *
     * A structure for a demo problem.
     *
     */
    public static class Zproblem {
        public DZcs A;
        public DZcs C;
        public int sym;
        public Complex[] x;
        public Complex[] b;
        public Complex[] resid;

        public Zproblem() {

        }
    };

    /* 1 if A is square & upper tri., -1 if square & lower tri., 0 otherwise */
    private static int is_sym(DZcs A) {
        int j, p, n = A.n, m = A.m, Ap[] = A.p, Ai[] = A.i;
        boolean is_upper, is_lower;
        if (m != n)
            return (0);
        is_upper = true;
        is_lower = true;
        for (j = 0; j < n; j++) {
            for (p = Ap[j]; p < Ap[j + 1]; p++) {
                if (Ai[p] > j)
                    is_upper = false;
                if (Ai[p] < j)
                    is_lower = false;
            }
        }
        return (is_upper ? 1 : (is_lower ? -1 : 0));
    }

    /* true for off-diagonal entries */
    private static class Dropdiag implements DZcs_ifkeep {

        @Override
        public boolean fkeep(int i, int j, Complex aij, Object other) {
            return (i != j);
        }

    }

    /* C = A + triu(A,1)' */
    private static DZcs make_sym(DZcs A) {
        DZcs AT, C;
        AT = DZcs_transpose.cs_transpose(A, true); /* AT = A' */
        DZcs_fkeep.cs_fkeep(AT, new Dropdiag(), null); /* drop diagonal entries from AT */
        C = DZcs_add.cs_add(A, AT, Complex.ONE, Complex.ONE); /* C = A+AT */
        return (C);
    }

    /* create a right-hand side */
    private static void rhs(Complex[] x, Complex[] b, int m) {
        int i;
        for (i = 0; i < m; i++)
            b[i] = new Complex(1 + ((double) i) / m, 0.0);
        for (i = 0; i < m; i++)
            x[i] = b[i];
    }

    /* infinity-norm of x */
    private static double norm(Complex[] x, int n) {
        int i;
        double normx = 0;
        for (i = 0; i < n; i++)
            normx = Math.max(normx, x[i].abs());
        return (normx);
    }

    /* compute residual, norm(A*x-b,inf) / (norm(A,1)*norm(x,inf) + norm(b,inf)) */
    private static void print_resid(boolean ok, DZcs A, Complex[] x, Complex[] b, Complex[] resid) {
        int i, m, n;
        if (!ok) {
            System.out.print("    (failed)\n");
            return;
        }
        m = A.m;
        n = A.n;
        for (i = 0; i < m; i++)
            resid[i] = b[i].multiply(-1.0); /* resid = -b */
        DZcs_gaxpy.cs_gaxpy(A, x, resid); /* resid = resid + A*x  */
        System.out.print(String.format("resid: %8.2e\n", norm(resid, m)
                / ((n == 0) ? 1 : (DZcs_norm.cs_norm(A) * norm(x, n) + norm(b, m)))));
    }

    private static double tic() {
        return System.nanoTime();
    }

    private static double toc(double t) {
        double s = tic();
        return (Math.max(0, s - t)) / 1000000.0;
    }

    private static void print_order(int order) {
        switch (order) {
        case 0:
            System.out.print("natural    ");
            break;
        case 1:
            System.out.print("amd(A+A')  ");
            break;
        case 2:
            System.out.print("amd(S'*S)  ");
            break;
        case 3:
            System.out.print("amd(A'*A)  ");
            break;
        }
    }

    /**
     * Reads a problem from a file.
     *
     * @param fileName
     *            file name
     * @param tol
     *            drop tolerance
     * @return problem
     */
    public static Zproblem get_problem(String fileName, double tol) {
        DZcs T, A, C;
        int sym, m, n, mn, nz1, nz2;
        Zproblem Prob;
        Prob = new Zproblem();
        T = DZcs_load.cs_load(fileName); /* load triplet matrix T from a file */
        Prob.A = A = DZcs_compress.cs_compress(T); /* A = compressed-column form of T */
        T = null; /* clear T */
        DZcs_dupl.cs_dupl(A);
        Prob.sym = sym = is_sym(A); /* determine if A is symmetric */
        m = A.m;
        n = A.n;
        mn = Math.max(m, n);
        nz1 = A.p[n];
        DZcs_dropzeros.cs_dropzeros(A); /* drop zero entries */
        nz2 = A.p[n];
        if (tol > 0)
            DZcs_droptol.cs_droptol(A, tol); /* drop tiny entries (just to test) */
        Prob.C = C = sym != 0 ? make_sym(A) : A; /* C = A + triu(A,1)', or C=A */
        if (C == null)
            return (null);
        System.out.print(String.format("\n--- Matrix: %d-by-%d, nnz: %d (sym: %d: nnz %d), norm: %8.2e\n", m, n,
                A.p[n], sym, sym != 0 ? C.p[n] : 0, DZcs_norm.cs_norm(C)));
        if (nz1 != nz2)
            System.out.print(String.format("zero entries dropped: %d\n", nz1 - nz2));
        if (nz2 != A.p[n])
            System.out.print(String.format("tiny entries dropped: %d\n", nz2 - A.p[n]));
        Prob.b = new Complex[mn];
        Prob.x = new Complex[mn];
        Prob.resid = new Complex[mn];
        return Prob;
    }

    /**
     * Solves a linear system using Cholesky, LU, and QR, with various
     * orderings.
     *
     * @param Prob
     *            problem
     * @return true if successful, false on error
     */
    public static boolean demo2(Zproblem Prob) {
        DZcs A, C;
        Complex b[], x[], resid[];
        double t, tol;
        int k, m, n, order, nb, ns, r[], s[], rr[], sprank;
        boolean ok;
        DZcsd D;
        if (Prob == null)
            return (false);
        A = Prob.A;
        C = Prob.C;
        b = Prob.b;
        x = Prob.x;
        resid = Prob.resid;
        m = A.m;
        n = A.n;
        tol = Prob.sym != 0 ? 0.001 : 1; /* partial pivoting tolerance */
        D = DZcs_dmperm.cs_dmperm(C, 1); /* randomized dmperm analysis */
        if (D == null)
            return (false);
        nb = D.nb;
        r = D.r;
        s = D.s;
        rr = D.rr;
        sprank = rr[3];
        for (ns = 0, k = 0; k < nb; k++) {
            if ((r[k + 1] == r[k] + 1) && (s[k + 1] == s[k] + 1)) {
                ns++;
            }
        }
        System.out.print(String.format("blocks: %d singletons: %d structural rank: %d\n", nb, ns, sprank));
        D = null;
        for (order = 0; order <= 3; order += 3) /* natural and amd(A'*A) */
        {
            if (order == 0 && m > 1000)
                continue;
            System.out.print("QR   ");
            print_order(order);
            rhs(x, b, m); /* compute right-hand side */
            t = tic();
            ok = DZcs_qrsol.cs_qrsol(order, C, x); /* min norm(Ax-b) with QR */
            System.out.print(String.format("time: %8.2f ms ", toc(t)));
            print_resid(ok, C, x, b, resid); /* print residual */
        }
        if (m != n || sprank < n)
            return (true); /* return if rect. or singular*/
        for (order = 0; order <= 3; order++) /* try all orderings */
        {
            if (order == 0 && m > 1000)
                continue;
            System.out.print("LU   ");
            print_order(order);
            rhs(x, b, m); /* compute right-hand side */
            t = tic();
            ok = DZcs_lusol.cs_lusol(order, C, x, tol); /* solve Ax=b with LU */
            System.out.print(String.format("time: %8.2f ms ", toc(t)));
            print_resid(ok, C, x, b, resid); /* print residual */
        }
        if (Prob.sym == 0)
            return (true);
        for (order = 0; order <= 1; order++) /* natural and amd(A+A') */
        {
            if (order == 0 && m > 1000)
                continue;
            System.out.print("Chol ");
            print_order(order);
            rhs(x, b, m); /* compute right-hand side */
            t = tic();
            ok = DZcs_cholsol.cs_cholsol(order, C, x); /* solve Ax=b with Cholesky */
            System.out.print(String.format("time: %8.2f ms ", toc(t)));
            print_resid(ok, C, x, b, resid); /* print residual */
        }
        return (true);
    }

    /**
     * Cholesky update/downdate
     *
     * @param Prob
     *            problem
     * @return true if successful, false on error
     */
    public static boolean demo3(Zproblem Prob) {
        DZcs A, C, W = null, WW, WT, E = null, W2;
        int n, k, Li[], Lp[], Wi[], Wp[], p1, p2, p[] = null;
        boolean ok;
        Complex b[], x[], resid[], y[] = null, Lx[], Wx[], s;
        double t, t1;
        DZcss S = null;
        DZcsn N = null;
        if (Prob == null || Prob.sym == 0 || Prob.A.n == 0)
            return (false);
        A = Prob.A;
        C = Prob.C;
        b = Prob.b;
        x = Prob.x;
        resid = Prob.resid;
        n = A.n;
        if (Prob.sym == 0 || n == 0)
            return (true);
        rhs(x, b, n); /* compute right-hand side */
        System.out.print("\nchol then update/downdate ");
        print_order(1);
        y = new Complex[n];
        t = tic();
        S = DZcs_schol.cs_schol(1, C); /* symbolic Chol, amd(A+A') */
        System.out.print(String.format("\nsymbolic chol time %8.2f ms\n", toc(t)));
        t = tic();
        N = DZcs_chol.cs_chol(C, S); /* numeric Cholesky */
        System.out.print(String.format("numeric  chol time %8.2f ms\n", toc(t)));
        if (S == null || N == null)
            return (false);
        t = tic();
        DZcs_ipvec.cs_ipvec(S.pinv, b, y, n); /* y = P*b */
        DZcs_lsolve.cs_lsolve(N.L, y); /* y = L\y */
        DZcs_ltsolve.cs_ltsolve(N.L, y); /* y = L'\y */
        DZcs_pvec.cs_pvec(S.pinv, y, x, n); /* x = P'*y */
        System.out.print(String.format("solve    chol time %8.2f ms\n", toc(t)));
        System.out.println("original: ");
        print_resid(true, C, x, b, resid); /* print residual */
        k = n / 2; /* construct W  */
        W = DZcs_util.cs_spalloc(n, 1, n, true, false);
        Lp = N.L.p;
        Li = N.L.i;
        Lx = N.L.x;
        Wp = W.p;
        Wi = W.i;
        Wx = W.x;
        Wp[0] = 0;
        p1 = Lp[k];
        Wp[1] = Lp[k + 1] - p1;
        s = Lx[p1];
        Random r = new Random(1);
        for (; p1 < Lp[k + 1]; p1++) {
            p2 = p1 - Lp[k];
            Wi[p2] = Li[p1];
            Wx[p2] = s.multiply(r.nextDouble());
        }
        t = tic();
        ok = DZcs_updown.cs_updown(N.L, +1, W, S.parent); /* update: L*L'+W*W' */
        t1 = toc(t);
        System.out.print(String.format("update:   time: %8.2f ms\n", t1));
        if (!ok)
            return (false);
        t = tic();
        DZcs_ipvec.cs_ipvec(S.pinv, b, y, n); /* y = P*b */
        DZcs_lsolve.cs_lsolve(N.L, y); /* y = L\y */
        DZcs_ltsolve.cs_ltsolve(N.L, y); /* y = L'\y */
        DZcs_pvec.cs_pvec(S.pinv, y, x, n); /* x = P'*y */
        t = toc(t);
        p = DZcs_pinv.cs_pinv(S.pinv, n);
        W2 = DZcs_permute.cs_permute(W, p, null, true); /* E = C + (P'W)*(P'W)' */
        WT = DZcs_transpose.cs_transpose(W2, true);
        WW = DZcs_multiply.cs_multiply(W2, WT);
        WT = null;
        W2 = null;
        E = DZcs_add.cs_add(C, WW, Complex.ONE, Complex.ONE);
        WW = null;
        if (E == null || p == null)
            return (false);
        System.out.print(String.format("update:   time: %8.2f ms(incl solve) ", t1 + t));
        print_resid(true, E, x, b, resid); /* print residual */
        N = null; /* clear N */
        t = tic();
        N = DZcs_chol.cs_chol(E, S); /* numeric Cholesky */
        if (N == null)
            return (false);
        DZcs_ipvec.cs_ipvec(S.pinv, b, y, n); /* y = P*b */
        DZcs_lsolve.cs_lsolve(N.L, y); /* y = L\y */
        DZcs_ltsolve.cs_ltsolve(N.L, y); /* y = L'\y */
        DZcs_pvec.cs_pvec(S.pinv, y, x, n); /* x = P'*y */
        t = toc(t);
        System.out.print(String.format("rechol:   time: %8.2f ms(incl solve) ", t));
        print_resid(true, E, x, b, resid); /* print residual */
        t = tic();
        ok = DZcs_updown.cs_updown(N.L, -1, W, S.parent); /* downdate: L*L'-W*W' */
        t1 = toc(t);
        if (!ok)
            return (false);
        System.out.print(String.format("downdate: time: %8.2f\n", t1));
        t = tic();
        DZcs_ipvec.cs_ipvec(S.pinv, b, y, n); /* y = P*b */
        DZcs_lsolve.cs_lsolve(N.L, y); /* y = L\y */
        DZcs_ltsolve.cs_ltsolve(N.L, y); /* y = L'\y */
        DZcs_pvec.cs_pvec(S.pinv, y, x, n); /* x = P'*y */
        t = toc(t);
        System.out.print(String.format("downdate: time: %8.2f ms(incl solve) ", t1 + t));
        print_resid(true, C, x, b, resid); /* print residual */
        return (true);
    }

}
