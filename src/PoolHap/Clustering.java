package PoolHap;
import java.util.regex.Pattern;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.mllib.clustering.KMeans;
import org.apache.spark.mllib.clustering.KMeansModel;
import org.apache.spark.mllib.linalg.Vector;
import org.apache.spark.mllib.linalg.Vectors;
import scala.util.Random;
import weka.core.Instances;
/**
 * Example using MLlib KMeans from Java.
 */




public final class Clustering {
	
	public int getDistance(String str1, String str2) {
		int distance;
		if (str1.length() != str2.length()) {
			distance = -1;
		} else {
			distance = 0;
			for (int i = 0; i < str1.length(); i++) {
				if (str1.charAt(i) != str2.charAt(i)) {
					distance++;
				}
			}
		}
		return distance;
	}

	  private static class ParsePoint implements Function<String, Vector> {
	    private static final Pattern SPACE = Pattern.compile(" ");
	
	@Override
	public Vector call(String line) {
		String[] tok = SPACE.split(line);
		double[] point = new double[tok.length];
	  for (int i = 0; i < tok.length; ++i) {
	//    	  System.out.println("#########"+tok[i]);
	        point[i] = Double.parseDouble(tok[i]);
	      }
	      return Vectors.dense(point);
	    }
	  }
	
	  public static void main(String[] args) {
	
	  String inputFile = "/home/chencao/Desktop/kmeans.txt";
	int k = 2; // two clusters
	int iterations = 20;
	int runs = 1;
	
	JavaSparkContext sc = new JavaSparkContext("local", "JavaKMeans");
	JavaRDD<String> lines = sc.textFile(inputFile);
	
	JavaRDD<Vector> points = lines.map(new ParsePoint());
	
	KMeansModel model = KMeans.train(points.rdd(), k, iterations, runs, KMeans.K_MEANS_PARALLEL());
	
	System.out.println("Cluster centers:");
	for (Vector center : model.clusterCenters()) {
	  System.out.println(" " + center);
	}
	double cost = model.computeCost(points.rdd());
	System.out.println("Cost: " + cost);
	
	    sc.stop();
	}
}


