package data.t1d;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.StringTokenizer;

import utils.CollectTools;

import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

import cern.colt.list.DoubleArrayList;
import cern.jet.stat.Descriptive;
import cern.jet.stat.Probability;

import utils.AbstractLineProcessor;

@App(overview="Call differences in methylation array data using sign-rank method", generateStub=true)
public class PairwiseDMRCallerW extends AbstractLineProcessor {
	private int[] fg;
	private int[] bg;
	private int offset = 0;
	
	@Option(help="Offset to apply to column indices (default=0, other values not recommended)", optional=true)
	public void setOffset(int i) {
		this.offset = i;
	}
	
	@Option(help="List of zero-based column indices for foreground samples")
	public void setFg(String s) {
		fg = uncut(new String[] {s});
	}
	
	@Option(help="List of zero-based column indices for background samples")
	public void setBg(String s) {
		bg = uncut(new String[] {s});
	}
	
	private int[] uncut(String[] s) {
		List<Integer> il = new ArrayList<Integer>();
		for (String k : s) {
			StringTokenizer toke = new StringTokenizer(k, ",");
			while (toke.hasMoreTokens()) {
				String t = toke.nextToken();
				int c = t.indexOf('-');
				if (c > 0) {
					int min = Integer.parseInt(t.substring(0, c));
					int max = Integer.parseInt(t.substring(c + 1));
					for (int i = min; i <= max; ++i) {
						il.add(i);
					}
				} else {
					il.add(Integer.parseInt(t));
				}
			}
		}
		// Collections.sort(il);
		return CollectTools.toIntArray(il);
	}
	
	
	public void pre() {
		if (fg.length != bg.length) {
			throw new IllegalArgumentException("FG and BG samples must be paired!");
		}
		setProcessComments(true);
	}
	
	public void processTokens(String[] args) 
	{
		System.out.print(args[0]);
		DoubleArrayList dal = new DoubleArrayList();
		for (int i = 0; i < fg.length; ++i) {
			double diff = Double.parseDouble(args[fg[i] + offset]) - Double.parseDouble(args[bg[i] + offset]);
			if (!Double.isNaN(diff)) {
				dal.add(diff);
			}
		}
		
		try {
			int n = dal.size();
			double w = Wilcoxon.wstat(dal);
			double p = Wilcoxon.pTwoSided(w, dal.size());
			double dom = Descriptive.mean(dal);
			
			System.out.printf("\t%g\t%g\t%g%n", dom, p, w);
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}
}
