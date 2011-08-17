/*
 * Based on R signrank.c...
 *
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1999-2007  The R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double dsignrank(double x, double n, int give_log)
 *    double psignrank(double x, double n, int lower_tail, int log_p)
 *    double qsignrank(double x, double n, int lower_tail, int log_p)
 *    double rsignrank(double n)
 *
 *  DESCRIPTION
 *
 *    dsignrank	   The density of the Wilcoxon Signed Rank distribution.
 *    psignrank	   The distribution function of the Wilcoxon Signed Rank
 *		   distribution.
 *    qsignrank	   The quantile function of the Wilcoxon Signed Rank
 *		   distribution.
 *    rsignrank	   Random variates from the Wilcoxon Signed Rank
 *		   distribution.
 */

package data.t1d;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import utils.MathsTools;

import cern.colt.list.DoubleArrayList;

public class Wilcoxon {
	private static class Datum {
		private final double x;
		private final double ax;
		
		public Datum(double x) {
			this.x = x;
			this.ax = Math.abs(x);
		}
	}
	
	public static int wstat(DoubleArrayList dal) {
		List<Datum> dl = new ArrayList<Datum>();
		for (int i = 0; i < dal.size(); ++i) {
			dl.add(new Datum(dal.get(i)));
		}
		
		Collections.sort(dl, new Comparator<Datum>() {
			public int compare(Datum o1, Datum o2) {
				return MathsTools.sign(o1.ax - o2.ax);
			}
		});
		
		int tot = 0;
		for (int i = 0; i < dl.size(); ++i) {
			if (dl.get(i).x > 0) {
				tot += (i + 1);
			}
		}
			
		return tot;
	}
	
	public static double pTwoSided(double x, int n) {
		double p;
		if (x > (n * (n+1) / 4)) {
			p = psignrank(x - 1, n, false, false);
		} else {
			p = psignrank(x, n, false, false);
		}
		return Math.min(2 * p, 1);
	}
	
	public static void main(String[] args) {
		for (int i = 1; i <= 120; ++i) {
			System.out.println(pTwoSided(i, 15));
		}
	}
	
	private static final double M_LN2 = 0.693147180559945309417232121458;
	
	// Quick'n'dirty port of R code.
	// Horribly, ghastlily, non-thread-safe!!!!
	
	private static double[] w;
	private static int allocated_n;

	static void w_init_maybe(int n)
	{
	    int u, c;

	    u = n * (n + 1) / 2;
	    c = (u / 2);

	    if (w != null) {
	        if(n != allocated_n) {
	            w = null;
	        }
	        else return;
	    }

	    if(w == null) {
	        w = new double[c + 1];
	        allocated_n = n;
	    }
	}
	
	static double
	csignrank(int k, int n)
	{
	    int c, u, j;

	    u = n * (n + 1) / 2;
	    c = (u / 2);

	    if (k < 0 || k > u)
	        return 0;
	    if (k > c)
	        k = u - k;

	    if (n == 1)
	        return 1.;
	    if (w[0] == 1.)
	        return w[k];

	    w[0] = w[1] = 1.;
	    for(j=2; j < n+1; ++j) {
	        // int i, end = imin2(j*(j+1)/2, c);
	    	int i, end = Math.min(j*(j+1)/2, c);
	        for(i=end; i >= j; --i) {
	            w[i] += w[i-j];
	        }
	    }

	    return w[k];
	}

	static double dsignrank(double x, double n, int give_log)
	{
	    double d;

	    /* NaNs propagated correctly */
	    if (Double.isNaN(x) || Double.isNaN(n)) return(x + n);

	    n = Math.floor(n + 0.5);
	    if (n <= 0)
	        return Double.NaN;

	    if (Math.abs(x - Math.floor(x + 0.5)) > 1e-7)
	        return 0.0; // return(R_D__0);
	    x = Math.floor(x + 0.5);
	    if ((x < 0) || (x > (n * (n + 1) / 2)))
	        return 0.0; // return(R_D__0);

	    w_init_maybe((int) n);
	    // d = R_D_exp(log(csignrank(x, n)) - n * M_LN2);
	    d = Math.exp(Math.log(csignrank((int) x, (int) n)) - n * M_LN2);
	    
	    return(d);
	}

	static double psignrank(double x, double n, boolean lower_tail, boolean log_p)
	{
	    int i;
	    double f, p;

	    if (Double.isNaN(x) || Double.isNaN(n))
	      return(x + n);

	    if (Double.isInfinite(n))
	    	return Double.NaN; 
	    n = Math.floor(n + 0.5);
	    if (n <= 0) 
	    	return Double.NaN;

	    x = Math.floor(x + 1e-7);
	    if (x < 0.0)
	        return 0;   // R_DT
	    if (x >= n * (n + 1) / 2)
	        return 1.0; // R_DT

	    w_init_maybe((int) n);
	    f = Math.exp(- n * M_LN2);
	    p = 0;
	    if (x <= (n * (n + 1) / 4)) {
	        for (i = 0; i <= x; i++)
	            p += csignrank((int) i, (int) n) * f;
	    }
	    else {
	        x = n * (n + 1) / 2 - x;
	        for (i = 0; i < x; i++)
	            p += csignrank((int) i, (int) n) * f;
	        lower_tail = !lower_tail; /* p = 1 - p; */
	    }

	    // return(R_DT_val(p));
	    return p;
	} /* psignrank() */

	static double qsignrank(double x, double n, boolean lower_tail, boolean log_p)
	{
	    double f, p, q;

	    if (Double.isNaN(x) || Double.isNaN(n))
	        return(x + n);
	    
	    
	    if (Double.isInfinite(x) || Double.isInfinite(n)) 
	    	return Double.NaN;
	    // if (!R_FINITE(x) || !R_FINITE(n))
	    //   ML_ERR_return_NAN;
	    if (x < 0 || x >1) {
	    	return Double.NaN;
	    }

	    n = Math.floor(n + 0.5);
	    if (n <= 0)
	        return Double.NaN;

	    if (x == 0)
	        return(0);
	    if (x == 1)
	        return(n * (n + 1) / 2);

	    if(log_p || !lower_tail) {
	    	throw new RuntimeException("FIXME");
	        // x = R_DT_qIv(x); /* lower_tail,non-log "p" */
	    }

	    w_init_maybe((int) n);
	    f = Math.exp(- n * M_LN2);
	    p = 0;
	    q = 0;
	    if (x <= 0.5) {
	        x = x - 10 * Double.MIN_VALUE;
	        for (;;) {
	            p += csignrank((int) q, (int) n) * f;
	            if (p >= x)
	                break;
	            q++;
	        }
	    }
	    else {
	        x = 1 - x + 10 * Double.MIN_VALUE;
	        for (;;) {
	            p += csignrank((int) q, (int) n) * f;
	            if (p > x) {
	                q = n * (n + 1) / 2 - q;
	                break;
	            }
	            q++;
	        }
	    }

	    return(q);
	}

	double rsignrank(double n)
	{
	    int i, k;
	    double r;

	    /* NaNs propagated correctly */
	    if (Double.isNaN(n)) return(n);

	    n = Math.floor(n + 0.5);
	    if (n < 0) 
	    	return Double.NaN;
	    if (n == 0)
	        return(0);
	    r = 0.0;
	    k = (int) n;
	    for (i = 0; i < k; ) {
	        r += (++i) * Math.floor(Math.random() + 0.5);
	    }
	    return(r);
	}
	
	
}
