package edu.emory.mathcs.csparsej.tdcomplex.test;

import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_qrsol.cs_qrsol;

import java.io.InputStream;

import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcs;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcsa;

public class Main extends DZcs_test {

	public static void main(String[] args) {

		int m ;
		double t ;
		DZcs A, C ;
		DZcsa b, x ;

		InputStream in = get_stream (C_IBM32B) ;
		DZproblem prob = get_problem (in, DROP_TOL) ;

		A = prob.A ; C = prob.C ; b = prob.b ; x = prob.x ;
		m = A.m ;

		rhs (x, b, m) ;

		t = tic() ;

		cs_qrsol (0, C, x) ;

		System.out.printf("time: %8.2f ms ", toc (t)) ;

	}

}
