import java.io.*;
import java.util.*;

public class WordCount
{
	public static class MotifInSeq
	{
		public int score, seq_id, start_index;
		public String candidate;
	}

	public static void main(String[] args) throws IOException {
		char m_char[] = "acgt".toCharArray();
		int[] m = new int[8];

		List<String> lines = Files.readAllLines(Paths.get("c:/temp/promoters_data_clean.txt"), Charset.defaultCharset());

		for(m[0]=0; m[0]<4; m[0]++)
			for(m[1]=0; m[1]<4; m[1]++)
				for(m[2]=0; m[2]<4; m[2]++)
					for(m[3]=0; m[3]<4; m[3]++)
						for(m[4]=0; m[4]<4; m[4]++)
							for(m[5]=0; m[5]<4; m[5]++)
								for(m[6]=0; m[6]<4; m[6]++)
									for(m[7]=0; m[7]<4; m[7]++)
									{
										ArrayList<MotifInSeq> cache = new ArrayList<MotifInSeq>();
										
										int total = 0;
										String motif = "";
										for(int i = 0; i < 8; i++)
											motif += m_char[m[i]];

										for(int l=0; l<lines.size(); l++)
										{
											MotifInSeq ms = new MotifInSeq();
											String str = lines.get(l);
											int str_len = str.length();
												
											int min_dist = 10;
											int min_index = 0;

											for(int i = 0; i < str_len - 8; i++)
											{
												int dist = 0;
												for(int j = 0; j < 8; j++)
												{
													if(str.charAt(i + j) != motif.charAt(j))
														dist++;
												}
												if(dist < min_dist)
												{
													min_dist = dist;
													min_index = i;
												}
												if(min_dist == 0)
													break;
											}
											
											ms.score = min_dist;
											ms.seq_id = l;
											ms.start_index = min_index;
											ms.candidate = str.substring(min_index, min_index + 8);
											cache.add(ms);
											
											total += min_dist;
										}
										
										Collections.shuffle(cache);
										for (MotifInSeq val : cache)
										{
											System.out.println(motif + "\t" + val.candidate + "\t" + val.seq_id + "\t" + val.score + "\t" + val.start_index + "\t" + total);
										}
									}
	}
}
