#!/bin/sh 

start=$(date +%s)
for dir in */; do
  if [ -d "$dir" ]; then
    cd "$dir"
    /Users/drake/opt/anaconda3/envs/strainline/strainline.sh -i . -o out -p ont --maxLD 0.01
    cd -
  fi
done
end=$(date +%s)
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Long reads Indeed. Time to relax with some music now!!"
cd /Users/drake/Music/Music/Unknown\ Artist/Unknown\ Album 
while true; do afplay "$(ls *.mp3 | shuf -n1)"; done

    