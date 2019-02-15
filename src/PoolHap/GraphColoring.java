package PoolHap;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;

public class GraphColoring
{
	public static void gc(String VefFile, HashMap<Integer, Integer> pos_dict, String OutputFile)
    throws IOException
  {
    int n_pos = pos_dict.size();
    Vector<Integer> readindex_arr_tmp = new Vector<Integer>();
    Vector<String> readinfo_arr_tmp = new Vector<String>();
    Vector<Integer> readindex_arr = new Vector<Integer>();
    Vector<String> readinfo_arr = new Vector<String>();
    HashMap<String, Integer> geno_dict = new HashMap<String, Integer>();
    
    int count = 0;
    FileReader vf = new FileReader(VefFile);
    BufferedReader bufferedreader = new BufferedReader(vf);
    String line = "";
    int max_num_geno = 32878;
    while ((line = bufferedreader.readLine()) != null)
    {
      line = line.replace("\r", "");
      String[] line_arr = line.split("\t");
      if (line_arr[1].contains("="))
      {
        String tmp_geno = line_arr[1];
        if (!geno_dict.containsKey(tmp_geno))
        {
          geno_dict.put(tmp_geno, Integer.valueOf(1));
          readindex_arr_tmp.add(Integer.valueOf(count));
          
          count++;
          readinfo_arr_tmp.add(line_arr[1]);
        }
        else
        {
          int tmp_num = ((Integer)geno_dict.get(tmp_geno)).intValue();
          geno_dict.remove(tmp_geno);
          geno_dict.put(tmp_geno, Integer.valueOf(tmp_num + 1));
          if (((Integer)geno_dict.get(tmp_geno)).intValue() < max_num_geno)
          {
            readindex_arr_tmp.add(Integer.valueOf(count));
            
            count++;
            readinfo_arr_tmp.add(line_arr[1]);
          }
        }
      }
    }
    bufferedreader.close();
    
    int[] index_arr_tmp = new int[readindex_arr_tmp.size()];
    for (int k = 0; k < index_arr_tmp.length; k++) {
      index_arr_tmp[k] = k;
    }
    ArrayList<Integer> list = new ArrayList<Integer>();
    for (int i = 0; i < index_arr_tmp.length; i++) {
      list.add(Integer.valueOf(index_arr_tmp[i]));
    }
    int[] index_arr = new int[readindex_arr_tmp.size()];
    Collections.shuffle(list);
    Iterator<Integer> ite = list.iterator();
    int tmp_i = 0;
    while (ite.hasNext())
    {
      index_arr[tmp_i] = ite.next();
      tmp_i++;
    }
    for (int i = 0; i < readindex_arr_tmp.size(); i++)
    {
      readinfo_arr.add((String)readinfo_arr_tmp.get(index_arr[i]));
      readindex_arr.add((Integer)readindex_arr_tmp.get(index_arr[i]));
    }
    ArrayList<ArrayList<Integer>> read_pos_2D_arr = new ArrayList<ArrayList<Integer>>();
    ArrayList<ArrayList<String>> read_geno_2D_arr = new ArrayList<ArrayList<String>>();
    for (int i = 0; i < readinfo_arr.size(); i++)
    {
      String tmp_str = (String)readinfo_arr.get(i);
      String[] tmp_str_arr = tmp_str.split(";");
      ArrayList<String> i_geno_arr = new ArrayList<String>();
      ArrayList<Integer> i_pos_arr = new ArrayList<Integer>();
      for (int k = 0; k < tmp_str_arr.length; k++)
      {
        String[] tmp2_arr = tmp_str_arr[k].split("=");
        int pos = Integer.parseInt(tmp2_arr[0]);
        i_pos_arr.add((Integer)pos_dict.get(Integer.valueOf(pos)));
        i_geno_arr.add(tmp2_arr[1]);
      }
      read_pos_2D_arr.add(i_pos_arr);
      read_geno_2D_arr.add(i_geno_arr);
    }
    ArrayList<ArrayList<Integer>> readadj_arr = new ArrayList<ArrayList<Integer>>();
    for (int i = 0; i < read_pos_2D_arr.size(); i++)	// This loop takes the longest because there are so many reads.
    {
      // if (i % 1000 == 0) System.out.println(i + " read-genotypes evaluated for confidence."); 
      ArrayList<Integer> adj_arr = new ArrayList<Integer>();
      for (int j = 0; j < read_pos_2D_arr.size(); j++) {
        if (i != j)
        {
          boolean IsConnect = false;
          for (int k = 0; k < (read_pos_2D_arr.get(i)).size(); k++) {
            for (int l = 0; l < (read_pos_2D_arr.get(j)).size(); l++) {
              if (((read_pos_2D_arr.get(i)).get(k) == (read_pos_2D_arr.get(j)).get(l)) && 
                (!(read_geno_2D_arr.get(i)).get(k).equals((read_geno_2D_arr.get(j)).get(l)))) {
                IsConnect = true;
              }
            }
          }
          if (IsConnect) {
            adj_arr.add((Integer)readindex_arr.get(j));
          }
        }
      }
      readadj_arr.add(adj_arr);
    }
    System.out.println("Finished assessing read genotypes.");
    int max_color = readindex_arr.size();
    ArrayList<Integer> read_color_arr = new ArrayList<Integer>();
    ArrayList<HashSet<Integer>> nb_color_arr = new ArrayList<HashSet<Integer>>();
    for (int i = 0; i < readadj_arr.size(); i++)
    {
      nb_color_arr.add(new HashSet<Integer>());
      read_color_arr.add(-1);
    }
    read_color_arr.set(0,0);
    for (int i = 0; i < (readadj_arr.get(0)).size(); i++)
    {
      int index = readadj_arr.get(0).get(i);
      nb_color_arr.get(index).add(0);
    }
    int real_max_color = 0;
    ArrayList<HashSet<String>> color_geno_set_arr = new ArrayList<HashSet<String>>();
    color_geno_set_arr.add(new HashSet<String>());
    color_geno_set_arr.get(0).add(readinfo_arr.get(0));
    for (;;)
    {
      int max_nb_color = -1;
      int index = -1;
      for (int i = 0; i < readadj_arr.size(); i++) {
        if (((read_color_arr.get(i)).intValue() == -1) && 
          ((nb_color_arr.get(i)).size() > max_nb_color))
        {
          index = i;
          max_nb_color = (nb_color_arr.get(i)).size();
        }
      }
      if (index == -1) {
        break;
      }
      int color = -1;
      for (int i = 0; i < max_color; i++) {
        if (!(nb_color_arr.get(index)).contains(i)) {
          if (i < color_geno_set_arr.size())
          {
            if (!(color_geno_set_arr.get(i)).contains(readinfo_arr.get(index)))
            {
              color = i;
              (color_geno_set_arr.get(i)).add(readinfo_arr.get(index));
              break;
            }
          }
          else
          {
            color_geno_set_arr.add(new HashSet<String>());
            (color_geno_set_arr.get(color_geno_set_arr.size() - 1)).add(readinfo_arr.get(index));
            color = i;
            break;
          }
        }
      }
      if (color > real_max_color) {
        real_max_color = color;
      }
      read_color_arr.set(index, color);
      for (int i = 0; i < (readadj_arr.get(index)).size(); i++) {
        nb_color_arr.get(readadj_arr.get(index).get(i)).add(color);
      }
    }
    System.out.println("Finished assessing possible haplotypes.");
    String null_ref = "*";
    String conf_ref = "";
    for (int i = 1; i < n_pos; i++)
    {
      null_ref = null_ref + "*";
      conf_ref = conf_ref + "?";
    }
    String[] ref_arr = new String[real_max_color + 1];
    String[] conf_arr = new String[real_max_color + 1];
    for (int i = 0; i <= real_max_color; i++)
    {
      ref_arr[i] = null_ref;
      conf_arr[i] = conf_ref;
    }
    for (int i = 0; i < read_color_arr.size(); i++)
    {
      int i_color = ((Integer)read_color_arr.get(i)).intValue();
      for (int j = 0; j < (read_pos_2D_arr.get(i)).size(); j++)
      {
        int p = read_pos_2D_arr.get(i).get(j);
        p++;
        ref_arr[i_color] = 
          (ref_arr[i_color].substring(0, p - 1) + (read_geno_2D_arr.get(i)).get(j).toString() + ref_arr[i_color].substring(p, ref_arr[i_color].length()));
      }
      if ((read_pos_2D_arr.get(i)).size() > 1) {
        for (int j = 0; j < (read_pos_2D_arr.get(i)).size() - 1; j++) {
          if ((Math.abs(read_pos_2D_arr.get(i).get(j + 1)) - 
            read_pos_2D_arr.get(i).get(j)) == 1)
          {
            int p = read_pos_2D_arr.get(i).get(j);
            p++;
            conf_arr[i_color] = 
              (conf_arr[i_color].substring(0, p - 1) + "-" + conf_arr[i_color].substring(p, conf_arr[i_color].length()));
          }
        }
      }
    }
    System.out.println("Finished assessing confidences.");
    HashMap<String, Integer> output_ref_arr = new HashMap<String, Integer>();
    HashMap<String, String> conf_ref_arr = new HashMap<String, String>();
    for (int i = 0; i <= real_max_color; i++) {
      if (!ref_arr[i].contains("*")) {
        if (output_ref_arr.containsKey(ref_arr[i]))
        {
          output_ref_arr.put(ref_arr[i], output_ref_arr.get(ref_arr[i]) + 1);
          conf_ref_arr.put(ref_arr[i], conf_arr[i]);
        }
        else
        {
          output_ref_arr.put(ref_arr[i], 1);
          conf_ref_arr.put(ref_arr[i], conf_arr[i]);
        }
      }
    }
    FileWriter mydata = new FileWriter(OutputFile, false);
    PrintWriter pw = new PrintWriter(mydata);
    for (String key : output_ref_arr.keySet())
    {
      String x = key;
      String a = x;
      String b = conf_ref_arr.get(x);
      String c = "";
      for (int i = 0; i < b.length(); i++) {
        c = c + a.substring(i, i + 1) + b.substring(i, i + 1);
      }
      c = c + a.substring(a.length() - 1, a.length());
      pw.write(c + "\t" + output_ref_arr.get(x) + "\n");
    }
    pw.flush();
    pw.close();
  }
}
