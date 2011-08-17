package data.t1d;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import utils.MathsTools;

import org.bjv2.util.cli.App;

import utils.AbstractLineProcessor;

@App(overview="Quantile normalization", generateStub=true)
public class MatrixQuantile extends AbstractLineProcessor {
	private static class RankedDatum {
		public final double x;
		public int index;
		public int rank;
		
		public RankedDatum(double x, int index) {
			this.x = x;
			this.index = index;
		}
	}
	
	private List<String> keys = new ArrayList<String>();
	private List<List<RankedDatum>> data;
	
	public void pre() {
		setProcessComments(true);
	}
	
	public void processTokens(String[] t) {
		if (data == null) {
			data = new ArrayList<List<RankedDatum>>();
			for (int i = 1; i < t.length; ++i) {
				data.add(new ArrayList<RankedDatum>());
			}
		}
		
		keys.add(t[0]);
		for (int i = 0; i < data.size(); ++i) {
			data.get(i).add(new RankedDatum(Double.parseDouble(t[i + 1]), data.get(i).size()));
		}
	}

	public void post() {
		for (List<RankedDatum> dl : data) {
			Collections.sort(dl, new Comparator<RankedDatum>() {
				public int compare(RankedDatum o1, RankedDatum o2) {
					return MathsTools.sign(o1.x - o2.x);
				}
			});
			for (int i = 0; i < dl.size(); ++i) {
				dl.get(i).rank = i;
			}
		}
		
		List<Double> template = new ArrayList<Double>();
		for (int i = 0; i < data.get(0).size(); ++i) {
			double tot = 0;
			int cnt = 0;
			for (int d = 0; d < data.size(); ++d) {
				tot += data.get(d).get(i).x;
				++cnt;
			}
			template.add(tot / cnt);
		}
		
		for (List<RankedDatum> dl : data) {
			Collections.sort(dl, new Comparator<RankedDatum>() {
				public int compare(RankedDatum o1, RankedDatum o2) {
					return o1.index - o2.index;
				}
			});
		}
		
		for (int i = 0; i < template.size(); ++i) {
			System.out.print(keys.get(i));
			for (int d = 0; d < data.size(); ++d) {
				System.out.print('\t');
				System.out.print(template.get(data.get(d).get(i).rank).doubleValue());
			}
			System.out.println();
		}
	}
}
