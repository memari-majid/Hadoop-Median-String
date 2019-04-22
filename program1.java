import java.io.IOException;
import java.util.StringTokenizer;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;

public class program1
{
	public static class TokenizerMapper extends Mapper<Object, Text, Text, IntWritable>
	{
		private IntWritable value = new IntWritable(0);
		private Text        word  = new Text();

		public void map(Object key, Text input, Context context) throws IOException, InterruptedException
		{
			String str = input.toString();
			int str_len = str.length();
			char m_char[] = "acgt".toCharArray();
			int[] m = new int[8];

			for(m[0]=0; m[0]<4; m[0]++)
				for(m[1]=0; m[1]<4; m[1]++)
					for(m[2]=0; m[2]<4; m[2]++)
						for(m[3]=0; m[3]<4; m[3]++)
							for(m[4]=0; m[4]<4; m[4]++)
								for(m[5]=0; m[5]<4; m[5]++)
									for(m[6]=0; m[6]<4; m[6]++)
										for(m[7]=0; m[7]<4; m[7]++)
										{
											int min_dist = 10;
											String motif = "";
											for(int i = 0; i < 8; i++)
												motif += m_char[m[i]];

											for(int i = 0; i < str_len - 8; i++)
											{
												int dist = 0;
												for(int j = 0; j < 8; j++)
												{
													if(str.charAt(i + j) != motif.charAt(j))
														dist++;
												}
												if(dist < min_dist)
													min_dist = dist;
												if(min_dist == 0)
													break;
											}
											word.set(motif);
											value.set(min_dist);
											context.write(word, value);
										}
		}
	}

	public static class IntSumReducer extends Reducer<Text, IntWritable, Text, IntWritable> 
	{
		private IntWritable result = new IntWritable();

		public void reduce(Text key, Iterable<IntWritable> values, Context context) throws IOException, InterruptedException 
		{
			int sum = 0;
			for (IntWritable val : values)
			{
				sum += val.get();
			}
			result.set(sum);
			context.write(key, result);
		}
	}

	public static void main(String[] args) throws Exception
	{
		Configuration conf = new Configuration();
		Job job = Job.getInstance(conf, "word count");

		job.setJarByClass(WordCount.class);

		job.setMapperClass(TokenizerMapper.class);
		job.setCombinerClass(IntSumReducer.class);
		job.setReducerClass(IntSumReducer.class);

		job.setOutputKeyClass(Text.class);
		job.setOutputValueClass(IntWritable.class);

		FileInputFormat.addInputPath(job, new Path(args[0]));
		FileOutputFormat.setOutputPath(job, new Path(args[1]));
		System.exit(job.waitForCompletion(true) ? 0 : 1);
	}
}
