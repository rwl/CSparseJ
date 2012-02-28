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

package edu.emory.mathcs.csparsej.tdcomplex.test ;

import edu.emory.mathcs.csparsej.tdcomplex.test.DZcs_test.DZproblem;

/**
 * Read a matrix, solve a linear system, update/downdate.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class DZcs_test3 extends DZcs_test {

	public static void main(String[] args)
	{
		DZproblem Prob = null ;
		if (args.length == 0)
		{
			throw new IllegalArgumentException("Usage: java edu.emory.mathcs.csparsej.tdcomplex.demo.DZcs_demo3 fileName") ;
		}
//		Prob = get_problem (args[0], 1e-14) ;
//		demo3 (Prob) ;
	}

}
