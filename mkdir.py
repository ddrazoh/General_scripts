for f in `ls -F ${in_directory} | grep /`;
      do 
        cd "${f}"
        mkdir attle
        cd ..;
         
      done
      
for f in `ls -F ${in_directory} | grep /`;
      do 
      /Users/drake/RATTLE/rattle cluster -i G* -t 24 --iso -o rattle;
      done
      
#works right
for f in `ls -F ${in_directory} | grep /`;
      do
        cd "${f}"
        mkdir rattle

        /Users/drake/RATTLE/rattle cluster -i G* -t 24 --iso -o rattle
        /Users/drake/RATTLE/rattle correct -i *.fastq  -c ${f}/${f}/rattle/clusters.out  -t 24 -o rattle
        cd ..;
      done
