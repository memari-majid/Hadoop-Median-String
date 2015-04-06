package cs591findingmotifproblem;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

public class FindingMotif
{
  public final char[] DNAChar = { 'a', 'c', 'g', 't' };
  public String[] DNA;
  public int motifLength;
  public int SeqNum;
  
  public FindingMotif(String filename, int SeqNum, int motifLength)
    throws FileNotFoundException
  {
    Scanner read = new Scanner(new File(filename));
    this.DNA = new String[SeqNum];
    this.SeqNum = SeqNum;
    this.motifLength = motifLength;
    int i = 0;
    while ((read.hasNextLine()) && (i < SeqNum))
    {
      this.DNA[i] = read.nextLine();
      i++;
    }
  }
  
  public Pair<String, Double> findingMotifAllalgorithms(int algorithm)
  {
    if (algorithm == 0)
    {
      long starttime = System.nanoTime();
      String motif = FromfinalPositionsToMotif(BruteForceMotifSearchAgain(this.DNA, this.SeqNum, this.DNA[0].length(), this.motifLength), this.DNA, this.motifLength);
      long endtime = System.nanoTime();
      double time = (endtime - starttime) / 1000000000.0D;
      return new Pair(motif, Double.valueOf(time));
    }
    if (algorithm == 1)
    {
      long starttime = System.nanoTime();
      String motif = FromfinalPositionsToMotif(BranchAndBoundMotifSearch(this.DNA, this.SeqNum, this.DNA[0].length(), this.motifLength), this.DNA, this.motifLength);
      long endtime = System.nanoTime();
      double time = (endtime - starttime) / 1000000000.0D;
      return new Pair(motif, Double.valueOf(time));
    }
    if (algorithm == 2)
    {
      long starttime = System.nanoTime();
      String motif = BruteForceMedianSearch(this.DNA, this.SeqNum, this.DNA[0].length(), this.motifLength);
      long endtime = System.nanoTime();
      double time = (endtime - starttime) / 1000000000.0D;
      return new Pair(motif, Double.valueOf(time));
    }
    if (algorithm == 3)
    {
      long starttime = System.nanoTime();
      String motif = BranchAndBoundMedianSearch(this.DNA, this.SeqNum, this.DNA[0].length(), this.motifLength);
      long endtime = System.nanoTime();
      double time = (endtime - starttime) / 1000000000.0D;
      return new Pair(motif, Double.valueOf(time));
    }
    long starttime = System.nanoTime();
    String motif = TighterBoundForBranchAndBoundMedianSearch(this.DNA, this.SeqNum, this.DNA[0].length(), this.motifLength);
    long endtime = System.nanoTime();
    double time = (endtime - starttime) / 1000000000.0D;
    return new Pair(motif, Double.valueOf(time));
  }
  
  public String TighterBoundForBranchAndBoundMedianSearch(String[] DNA, int t, int n, int l)
  {
    int[] s = new int[l];
    
    int bestDistance = 2147483647;
    int i = -1;
    int distance = 0;int optimisticSuffixDistance = 0;int optimisticPrefixDistance = 0;
    String bestWord = new String();
    
    int[] bestSubStrings = new int[l - 1];
    for (int j = 0; j < bestSubStrings.length; j++) {
      bestSubStrings[j] = 2147483647;
    }
    while (i > -2) {
      if (i < l - 1)
      {
        if (i != -1)
        {
          String prefix = NumbersToWord(s, i + 1);
          optimisticPrefixDistance = TotalDistance(prefix, DNA, t, n - l + 1);
          if (optimisticPrefixDistance < bestSubStrings[i]) {
            bestSubStrings[i] = optimisticPrefixDistance;
          }
          if (l - 1 - i < i) {
            optimisticSuffixDistance = bestSubStrings[(l - 1 - i)];
          } else {
            optimisticSuffixDistance = 0;
          }
        }
        if (optimisticPrefixDistance + optimisticSuffixDistance >= bestDistance)
        {
          Pair<int[], Integer> p = ByPass(s, i, l, 3);
          i = ((Integer)p.u).intValue();
          System.arraycopy(p.t, 0, s, 0, l);
        }
        else
        {
          Pair<int[], Integer> p = NextVertex(s, i, l, 3);
          i = ((Integer)p.u).intValue();
          System.arraycopy(p.t, 0, s, 0, l);
        }
      }
      else
      {
        String word = NumbersToWord(s, l);
        distance = TotalDistance(word, DNA, t, n - l + 1);
        if (distance < bestDistance)
        {
          bestDistance = distance;
          bestWord = word;
        }
        Pair<int[], Integer> p = NextVertex(s, i, l, 3);
        i = ((Integer)p.u).intValue();
        System.arraycopy(p.t, 0, s, 0, l);
      }
    }
    return bestWord;
  }
  
  public String BranchAndBoundMedianSearch(String[] DNA, int t, int n, int l)
  {
    int[] s = new int[l];
    
    int bestDistance = 2147483647;
    int i = -1;
    int distance = 0;int optimisticDistance = 0;
    String bestWord = new String();
    while (i > -2) {
      if (i < l - 1)
      {
        if (i != -1)
        {
          String prefix = NumbersToWord(s, i + 1);
          optimisticDistance = TotalDistance(prefix, DNA, t, n - l + 1);
        }
        if (optimisticDistance > bestDistance)
        {
          Pair<int[], Integer> p = ByPass(s, i, l, 3);
          i = ((Integer)p.u).intValue();
          System.arraycopy(p.t, 0, s, 0, l);
        }
        else
        {
          Pair<int[], Integer> p = NextVertex(s, i, l, 3);
          i = ((Integer)p.u).intValue();
          System.arraycopy(p.t, 0, s, 0, l);
        }
      }
      else
      {
        String word = NumbersToWord(s, l);
        distance = TotalDistance(word, DNA, t, n - l + 1);
        if (distance < bestDistance)
        {
          bestDistance = distance;
          bestWord = word;
        }
        Pair<int[], Integer> p = NextVertex(s, i, l, 3);
        i = ((Integer)p.u).intValue();
        System.arraycopy(p.t, 0, s, 0, l);
      }
    }
    return bestWord;
  }
  
  public String BruteForceMedianSearch(String[] DNA, int t, int n, int l)
  {
    int[] s = new int[l];
    

    int bestDistance = TotalDistance(NumbersToWord(s, l), DNA, t, n - l + 1);
    String bestWord = NumbersToWord(s, l);
    do
    {
      s = Nextleaf(s, l, 3);
      String word = NumbersToWord(s, l);
      int distance = TotalDistance(word, DNA, t, n - l + 1);
      if (distance < bestDistance)
      {
        bestDistance = distance;
        bestWord = word;
      }
    } while (!isZeroArray(s));
    return bestWord;
  }
  
  public int[] BranchAndBoundMotifSearch(String[] DNA, int t, int n, int l)
  {
    int[] s = new int[t];
    int[] bestMotif = new int[t];
    int bestScore = 0;
    int i = -1;
    int score = 0;int optimisticScore = 0;
    while (i > -2) {
      if (i < t - 1)
      {
        if (i != -1) {
          optimisticScore = Score(s, DNA, l, i + 1) + (t - (i + 1)) * l;
        }
        if (optimisticScore < bestScore)
        {
          Pair<int[], Integer> p = ByPass(s, i, t, n - l);
          i = ((Integer)p.u).intValue();
          System.arraycopy(p.t, 0, s, 0, t);
        }
        else
        {
          Pair<int[], Integer> p = NextVertex(s, i, t, n - l);
          i = ((Integer)p.u).intValue();
          System.arraycopy(p.t, 0, s, 0, t);
        }
      }
      else
      {
        score = Score(s, DNA, l, t);
        if (score > bestScore)
        {
          bestScore = score;
          System.arraycopy(s, 0, bestMotif, 0, t);
        }
        Pair<int[], Integer> p = NextVertex(s, i, t, n - l);
        i = ((Integer)p.u).intValue();
        System.arraycopy(p.t, 0, s, 0, t);
      }
    }
    return bestMotif;
  }
  
  public int[] BruteForceMotifSearchAgain(String[] DNA, int t, int n, int l)
  {
    int[] s = new int[t];
    int bestScore = Score(s, DNA, l, t);
    int[] bestMotif = new int[t];
    System.arraycopy(s, 0, bestMotif, 0, t);
    do
    {
      s = Nextleaf(s, t, n - l);
      int var = Score(s, DNA, l, t);
      if (var > bestScore)
      {
        bestScore = var;
        System.arraycopy(s, 0, bestMotif, 0, t);
      }
    } while (!isZeroArray(s));
    return bestMotif;
  }
  
  public int TotalDistance(String word, String[] DNA, int t, int k)
  {
    int totalDistance = 0;
    for (int i = 0; i < t; i++)
    {
      int minDistance = 2147483647;
      for (int j = 0; j < k; j++)
      {
        int distance = HammingDistance(word, DNA[i].substring(j, j + word.length()));
        if (distance < minDistance) {
          minDistance = distance;
        }
      }
      totalDistance += minDistance;
    }
    return totalDistance;
  }
  
  public int HammingDistance(String word1, String word2)
  {
    int distance = 0;
    for (int i = 0; i < word1.length(); i++) {
      if (word1.charAt(i) != word2.charAt(i)) {
        distance++;
      }
    }
    return distance;
  }
  
  public String NumbersToWord(int[] s, int upto)
  {
    StringBuffer word = new StringBuffer();
    for (int i = 0; i < upto; i++)
    {
      char c = s[i] == 2 ? this.DNAChar[2] : s[i] == 1 ? this.DNAChar[1] : s[i] == 0 ? this.DNAChar[0] : this.DNAChar[3];
      word.append(c);
    }
    return word.toString();
  }
  
  public String printSeq()
  {
    StringBuffer result = new StringBuffer();
    for (int i = 0; i < this.DNA.length; i++) {
      result.append(this.DNA[i]).append("\n");
    }
    return result.toString();
  }
  
  public String FromfinalPositionsToMotif(int[] s, String[] DNA, int l)
  {
    StringBuffer motifString = new StringBuffer();
    int[] s1 = new int[s.length];
    System.arraycopy(s, 0, s1, 0, s.length);
    for (int i = 0; i < l; i++)
    {
      int A = 0;int T = 0;int C = 0;int G = 0;
      for (int j = 0; j < s.length; j++)
      {
        if (DNA[j].charAt(s1[j]) == 'a') {
          A++;
        } else if (DNA[j].charAt(s1[j]) == 'c') {
          C++;
        } else if (DNA[j].charAt(s1[j]) == 'g') {
          G++;
        } else {
          T++;
        }
        s1[j] += 1;
      }
      if ((A >= T) && (A >= G) && (A >= C)) {
        motifString.append('a');
      } else if ((T >= A) && (T >= G) && (T >= C)) {
        motifString.append('t');
      } else if ((G >= A) && (G >= T) && (G >= C)) {
        motifString.append('g');
      } else {
        motifString.append('c');
      }
    }
    return motifString.toString();
  }
  
  public boolean isZeroArray(int[] s)
  {
    for (int i = 0; i < s.length; i++) {
      if (s[i] != 0) {
        return false;
      }
    }
    return true;
  }
  
  public Pair<int[], Integer> ByPass(int[] s, int i, int t, int k)
  {
    for (int j = i; j >= 0; j--) {
      if (s[j] < k)
      {
        s[j] += 1;
        return new Pair(s, Integer.valueOf(j));
      }
    }
    return new Pair(s, Integer.valueOf(-2));
  }
  
  public Pair<int[], Integer> NextVertex(int[] s, int i, int t, int k)
  {
    if (i < t - 1)
    {
      s[(i + 1)] = 0;
      return new Pair(s, Integer.valueOf(i + 1));
    }
    for (int j = t - 1; j >= 0; j--) {
      if (s[j] < k)
      {
        s[j] += 1;
        return new Pair(s, Integer.valueOf(j));
      }
    }
    return new Pair(s, Integer.valueOf(-2));
  }
  
  public int[] Nextleaf(int[] s, int t, int k)
  {
    for (int i = t - 1; i >= 0; i--)
    {
      if (s[i] < k)
      {
        s[i] += 1;
        return s;
      }
      s[i] = 0;
    }
    return s;
  }
  
  public int Score(int[] s, String[] DNA, int motifLength, int upto)
  {
    int totalscore = 0;
    int[] s1 = new int[s.length];
    System.arraycopy(s, 0, s1, 0, s.length);
    for (int i = 0; i < motifLength; i++)
    {
      int A = 0;int T = 0;int C = 0;int G = 0;
      for (int j = 0; j < upto; j++)
      {
        if (DNA[j].charAt(s1[j]) == 'a') {
          A++;
        } else if (DNA[j].charAt(s1[j]) == 'c') {
          C++;
        } else if (DNA[j].charAt(s1[j]) == 'g') {
          G++;
        } else {
          T++;
        }
        s1[j] += 1;
      }
      totalscore += Math.max(Math.max(C, G), Math.max(A, T));
    }
    return totalscore;
  }
}
