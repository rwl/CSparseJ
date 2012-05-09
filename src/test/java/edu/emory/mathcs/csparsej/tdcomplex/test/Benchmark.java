package edu.emory.mathcs.csparsej.tdcomplex.test;

import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_qrsol.cs_qrsol;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_lusol.cs_lusol;

import java.io.InputStream;

import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcs;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcsa;

public class Benchmark extends DZcs_test {

	private static final int N = 100 ;

	private static final int ORDER = 0 ;

	public static void main(String[] args) {

		int m ;
		double tol ;
		DZcs A, C ;
		DZcsa b, x ;

		InputStream in ;
		DZproblem prob ;

		in = get_stream (C_IBM32B) ;
		prob = get_problem (in, DROP_TOL) ;

		A = prob.A ; C = prob.C ; b = prob.b ; x = prob.x ;
		m = A.m ;

		/* partial pivoting tolerance */
		tol = prob.sym != 0 ? 0.001 : 1 ;

		rhs (x, b, m) ;

		System.out.println("CSparseJ: c_ibm32a") ;

		benchmark(C, x, b, m, tol, ORDER) ;
	}

	private static void benchmark(DZcs C, DZcsa x, DZcsa b, int m, double tol, int order) {

		double t ;
		double [] t_lu = new double [N] ;
		double [] t_qr = new double [N] ;

		for (int i = 0; i < N; i++) {
			t = tic() ;
			cs_lusol (order, C, x, tol) ;
			t_lu [i] = toc (t) ;

			System.arraycopy(b.x, 0, x.x, 0, m) ;

			t = tic() ;
			cs_qrsol (order, C, x) ;
			t_qr [i] = toc (t) ;

			System.arraycopy(b.x, 0, x.x, 0, m) ;
		}

		System.out.printf("LU - min: %6.2f, max: %6.2f, avg: %6.2f",
				min (t_lu), max (t_lu), avg (t_lu)) ;
		System.out.println() ;

		System.out.printf("QR - min: %6.2f, max: %6.2f, avg: %6.2f",
				min (t_qr), max (t_qr), avg (t_qr)) ;
		System.out.println() ;
	}

	private static double min(double[] tt) {
		int l = tt.length ;
		assert l > 0 ;

		double tmin = tt [0] ;

		for (int i = 1 ; i < l ; i++) {
			if (tt [i] < tmin) tmin = tt [i] ;
		}

		return tmin ;
	}

	private static double max(double[] tt) {
		int l = tt.length ;
		assert l > 0 ;

		double tmax = tt [0] ;

		for (int i = 1 ; i < l ; i++) {
			if (tt [i] > tmax) tmax = tt [i] ;
		}

		return tmax ;
	}

	private static double avg(double[] tt) {
		int l = tt.length ;
		assert l > 0 ;

		double sum = 0.0 ;
		for (int i = 0; i < l; i++) {
			sum = sum + tt [i] ;
		}

		return sum / l ;
	}

}
