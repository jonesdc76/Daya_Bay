{
  SNF snf;
  snf.Configure("SNF.config");
  //table of reactor positions
  printf("\\begin{table}\n");
  printf("\\begin{center}\n");
  printf("\\caption{Surveyed position coordinates of detectors.}\n");
  printf("\\label{tab:detpos}\n");
  printf("\\begin{tabular}{|l|c|c|c|}\\hline");

  printf("Detector&x (m)&y (m)& z (m)\\\\\\hline\n");
  for(int i=0;i<snf.kDetector;++i)
    printf("%s&%0.3f&%0.3f&%0.3f\\\n",snf.GetDetectorName(i).Data(),snf.vDetPos[i][0],snf.vDetPos[i][1],snf.vDetPos[i][2]);
  printf("\\hline\n");
  printf("\\end{tabular}\n");
  printf("\\end{center}\n");
  printf("\\end{table}\n");
  

}
