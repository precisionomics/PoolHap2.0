package PoolHap; 

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Scanner;
import java.util.concurrent.TimeUnit;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class SATVertexColoring {

	public SATVertexColoring(String name,
			GraphUV G, int[] clique)
		{
			assert clique.length <= G.numVertices();
			assert clique.length > 0;
			
			this.name = name;
			this.G = G;
			this.clique = clique;
			this.colors = new int[G.numVertices()];
			
			this.timeout = Long.parseLong(System.getProperty(
				"reads.vc.sat.timeout", "3600"));
		}
		
		public int[] getColors()
		{
			return colors;
		}
		
		public boolean run(int n, String method, String minisat_path)
			throws Exception
		{
			if (method.equals("l2md"))
			{
				runTree(n, false, minisat_path);
				return runL2MD(n, true, minisat_path);
			}
			else if (method.equals("tree"))
			{
				return runTree(n, true, minisat_path);
			}
			else if (method.equals("r2md"))
			{
				return runR2MD(n, true, minisat_path);
			}
			else if (method.equals("md"))
			{
				return runMD(n, true, minisat_path);
			}
			else if (method.equals("nmd"))
			{
				return runNMD(n, true, minisat_path);
			}
			else
			{
				throw new Exception("No such method "+method);
			}
		}
		
		private boolean runL2MD(int n, boolean solve, String minisat_path)
			throws Exception
		{
			assert (clique == null) || (n >= clique.length);

			// linear-2 has 3 branches
			final int mdVars = (n+2)/3;
			
			// for muldirect
			int[] oneOr = new int[2+mdVars];
			Arrays.fill(oneOr, 1);
			oneOr[0] = oneOr[1] = 0;
			
			// generate mdVars*3 colors
			// will restrict to n using clauses
			int[][] colorAnds = new int[mdVars*3][2+mdVars];
			for (int i = 0; i < mdVars; ++i)
			{
				// branch 0x
				Arrays.fill(colorAnds[i], 0);
				colorAnds[i][0] = -1;
				colorAnds[i][i+2] = 1;
				
				// branch 10
				Arrays.fill(colorAnds[i+mdVars], 0);
				colorAnds[i+mdVars][0] = 1;
				colorAnds[i+mdVars][1] = -1;
				colorAnds[i+mdVars][i+2] = 1;
				
				// branch 11
				Arrays.fill(colorAnds[i+2*mdVars], 0);
				colorAnds[i+2*mdVars][0] = 1;
				colorAnds[i+2*mdVars][1] = 1;
				colorAnds[i+2*mdVars][i+2] = 1;
			}
			
			String cnf = writeCNF("l2md", n, 2+mdVars, colorAnds, oneOr);
			
			if (solve)
				return runSat(cnf, n, 2+mdVars, colorAnds, minisat_path);
			else
				return true;
		}

		private boolean runTree(int n, boolean solve, String minisat_path)
			throws Exception
		{
			assert (clique == null) || (n >= clique.length);

			// as many as log2 variables
			int vars = 1;
			for (; (1 << vars) < n; ++vars)
				;
			
			// generate n colors
			int[][] colorAnds = new int[n][vars];
			encodeTreeColors(0, 0, n, colorAnds);
			
			String cnf = writeCNF("tree", n, vars, colorAnds, null);
			
			if (solve)
				return runSat(cnf, n, vars, colorAnds, minisat_path);
			else
				return true;
		}
		
		private void encodeTreeColors(
			int cur, int begin, int n, int[][] colorAnds)
		{
			assert n > 1;
			
			int next0 = n/2;
			int next1 = n-next0;
			
			for (int i = begin; i < begin+next0; ++i)
				colorAnds[i][cur] = -1;
			for (int i = begin+next0; i < begin+n; ++i)
				colorAnds[i][cur] = 1;
			
			if (next0 > 1)
				encodeTreeColors(cur+1,
					begin, next0, colorAnds);
			if (next1 > 1)
				encodeTreeColors(cur+1,
					begin+next0, next1, colorAnds);
		}
		
		private boolean runR2MD(int n, boolean solve, String minisat_path)
			throws Exception
		{
			assert (clique == null) || (n >= clique.length);

			// reduced-2 has 3 branches
			final int mdVars = (n+2)/3;
			
			// for muldirect
			int[] oneOr = new int[2+mdVars];
			Arrays.fill(oneOr, 1);
			oneOr[0] = oneOr[1] = 0;
			
			// generate mdVars*3 colors
			// will restrict to n using clauses
			int[][] colorAnds = new int[mdVars*3][2+mdVars];
			for (int i = 0; i < mdVars; ++i)
			{
				// branch 1x
				Arrays.fill(colorAnds[i], 0);
				colorAnds[i][0] = 1;
				colorAnds[i][i+2] = 1;
				
				// branch x1
				Arrays.fill(colorAnds[i+mdVars], 0);
				colorAnds[i+mdVars][1] = 1;
				colorAnds[i+mdVars][i+2] = 1;
				
				// branch 00
				Arrays.fill(colorAnds[i+2*mdVars], 0);
				colorAnds[i+2*mdVars][0] = -1;
				colorAnds[i+2*mdVars][1] = -1;
				colorAnds[i+2*mdVars][i+2] = 1;
			}
			
			String cnf = writeCNF("r2md", n, 2+mdVars, colorAnds, oneOr);
			
			if (solve)
				return runSat(cnf, n, 2+mdVars, colorAnds, minisat_path);
			else
				return true;
		}

		private boolean runMD(int n, boolean solve, String minisat_path)
			throws Exception
		{
			assert (clique == null) || (n >= clique.length);

			int[] oneOr = new int[n];
			Arrays.fill(oneOr, 1);
			
			// generate mdVars*3 colors
			// will restrict to n using clauses
			int[][] colorAnds = new int[n][n];
			for (int i = 0; i < n; ++i)
			{
				Arrays.fill(colorAnds[i], 0);
				colorAnds[i][i] = 1;
			}
			
			String cnf = writeCNF("md", n, n, colorAnds, oneOr);
			
			if (solve)
				return runSat(cnf, n, n, colorAnds, minisat_path);
			else
				return true;
		}

		private boolean runNMD(int n, boolean solve, String minisat_path)
			throws Exception
		{
			assert (clique == null) || (n >= clique.length);

			int[] oneOr = new int[n];
			Arrays.fill(oneOr, -1);
			
			// generate mdVars*3 colors
			// will restrict to n using clauses
			int[][] colorAnds = new int[n][n];
			for (int i = 0; i < n; ++i)
			{
				Arrays.fill(colorAnds[i], 0);
				colorAnds[i][i] = -1;
			}
			
			String cnf = writeCNF("nmd", n, n, colorAnds, oneOr);
			
			if (solve)
				return runSat(cnf, n, n, colorAnds, minisat_path);
			else
				return true;
		}
		
		private boolean runSat(String cnf,
			int n, int l, int[][] colorAnds, String minisat_path)
			throws Exception
		{
			Timer t = new Timer();
			
			System.out.flush();
			System.err.flush();
			
			try
			{
				Process p = new ProcessBuilder(
					minisat_path, cnf, cnf+".output")
					.inheritIO().start();
				
				if (!p.waitFor(timeout, TimeUnit.SECONDS))
					p.destroyForcibly();
			}
			catch (Exception e)
			{
				logger.info(
					"SAT: minisat failed: {}",
					e.toString());
			}
		
			System.out.flush();
			System.err.flush();
		
			boolean ret = readOutput(
				new File(cnf+".output"), n, l, colorAnds);
			
			System.out.printf(
				"@sat %s, real_sat %.3f%n",
				ret? "sat": "unsat", t.now());
			
			logger.info(
				"SAT done: t {}, {}",
				t.now(), ret? "sat": "unsat");

			return ret;
		}
		
		private String writeCNF(String method,
			int n, int l, int[][] colorAnds, int[] oneOr)
			throws Exception
		{
			assert colorAnds.length >= n;
			assert (oneOr == null) || (oneOr.length == l);
			
			String cnf = String.format("%s.%d.%s", name, n, method);
			
			logger.info(
				"SAT {}: V/E {}/{}, n/l {}/{}",
				method, G.numVertices(), G.numEdges(), n, l);
			
			for (int c = 0; c < colorAnds.length; ++c)
			{
				assert colorAnds[c].length == l;
				logger.debug("SAT {} color: c {}, {}",
					method, c, Arrays.toString(colorAnds[c]));
			}
			
			File file = new File(cnf);
			
			try(
				FileOutputStream f = new FileOutputStream(file);
				OutputStreamWriter osr = new OutputStreamWriter(f, "UTF-8");
				PrintWriter pw = new PrintWriter(osr))
			{
				pw.printf("c %s%n", file.getName());
				pw.println("c");

				for (int v = 0; v < G.numVertices(); ++v)
				{
					pw.printf("c %d:", v);
					for (int i = 0; i < l; ++i)
						pw.printf(" %d", v*l+i+1);
					pw.println("");
				}
				pw.println("c");
				
				int numClique = (clique == null)? 0: l*clique.length;
				int numAtLeastOne = (oneOr == null)? 0: G.numVertices();
				int numRestrict = (colorAnds.length-n)*G.numVertices();
				int numConflict = n*G.numEdges();
				
				pw.printf("p cnf %d %d%n",
					l*G.numVertices(),
					numClique+numAtLeastOne+numRestrict+numConflict);

				// assign color to known clique
				if (clique != null)
					for (int c = 0; c < clique.length; ++c)
						assignColor(pw, clique[c], colorAnds[c]);
				
				for (int v = 0; v < G.numVertices(); ++v)
				{
					// write at least one
					if (oneOr != null)
						atLeastOne(pw, v, oneOr);
					
					// restrict color
					for (int c = n; c < colorAnds.length; ++c)
					{
						writeNot(pw, v, colorAnds[c]);
						pw.println("0");
					}
				}
				
				// conflict along edges 
				for (int e = 0; e < G.numEdges(); ++e)
					for (int c = 0; c < n; ++c)
					{
						writeNot(pw, G.fromVertex(e), colorAnds[c]);
						writeNot(pw, G.toVertex(e), colorAnds[c]);
						pw.println("0");
					}
			}
			
			logger.info(
				"SAT {}: write {}, timeout {}",
				method, cnf, timeout);
			
			return cnf;
		}
		
		private void atLeastOne(PrintWriter pw, int v, int[] oneOr)
		{
			for (int i = 0; i < oneOr.length; ++i)
			{
				if (oneOr[i] == 0)
					continue;
				
				int var = v*oneOr.length+i+1;
				
				pw.print((oneOr[i] > 0)? var: -var);
				pw.print(" ");
			}
			pw.println("0");
		}
		
		private void writeNot(PrintWriter pw, int v, int[] andTerm)
		{
			for (int i = 0; i < andTerm.length; ++i)
			{
				if (andTerm[i] == 0)
					continue;
				
				int var = v*andTerm.length+i+1;
				
				pw.print((andTerm[i] > 0)? -var: var);
				pw.print(" ");
			}
		}

		private void assignColor(PrintWriter pw,
			int v, int[] colorAnd)
		{
			for (int i = 0; i < colorAnd.length; ++i)
			{
				if (colorAnd[i] == 0)
					continue;

				int var = v*colorAnd.length+i+1;
				
				pw.print((colorAnd[i] > 0)? var: -var);
				pw.println(" 0");
			}
		}
		
		private boolean readOutput(File file,
			int n, int l, int[][] colorAnds)
		{
			boolean[] vars = new boolean[l*G.numVertices()];
			
			try(Scanner s = new Scanner(file))
			{
				String sat = s.nextLine();
				
				if (!sat.startsWith("SAT"))
					return false;
				
				for (int i = 0; i < vars.length; ++i)
				{
					int var = s.nextInt();
					
					if ((var != (i+1)) && (var != (-i-1)))
						throw new Exception("Invalid var "+var);
					
					vars[i] = var > 0;
				}
				
				for (int v = 0; v < G.numVertices(); ++v)
				{
					colors[v] = findColor(v, n, colorAnds,
						Arrays.copyOfRange(vars, v*l, v*l+l));
					
					logger.debug("v {}, color {}", v, colors[v]);
				}
				
				return true;
			}
			catch (Exception e)
			{
				logger.info(
					"SAT: read output failed: {}",
					e.toString());
				
				return false;
			}
		}
		
		private int findColor(int v, int n,
			int[][] colorAnds, boolean[] vars)
			throws Exception
		{
			for (int c = 0; c < n; ++c)
			{
				if (matchColor(vars, colorAnds[c]))
					return c;
			}
			
			logger.debug("No color: v {}, {}",
				v, Arrays.toString(vars));
			
			throw new Exception("No color match");
		}
		
		private boolean matchColor(boolean[] vars, int[] colorAnd)
		{
			assert vars.length == colorAnd.length;
			
			for (int i = 0; i < colorAnd.length; ++i)
			{
				if (colorAnd[i] == 0)
					continue;
				
				if ((colorAnd[i] > 0) && !vars[i])
					return false;
				if ((colorAnd[i] < 0) && vars[i])
					return false;
			}
			
			return true;
		}
		
		private final String name;
		private final GraphUV G;
		private final int[] clique, colors;
		private final long timeout;
		
		private static final Logger logger
			= LoggerFactory.getLogger(SATVertexColoring.class);

}
