package Viral_Reconstructions_Tools;

import org.apache.hadoop.hive.ql.parse.HiveParser_IdentifiersParser.identifier_return;

public class test {

	public static void main(String[] args) {
		String d0String= "9.0E-4";
		double d=Double.parseDouble(d0String);
		if(d>0.01) {
			System.out.println(d);
		}

	}

}
