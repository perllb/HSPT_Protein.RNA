# Get all links for prot1

cut -f1 -d" " 9606.protein.links.v10.5.txt | sort | uniq -c | sed "s/[ \t]*//" | sed 's/ /,/g' | sed 's/9606.//g' > 9606.protein.links.v10.5.Freq.Prot1.txt

# Get all score > 400 links for prot1

awk ' { if ( $3 > 399 ) { print $0 } } '  9606.protein.links.v10.5.txt |  cut -f1 -d" " |  sort | uniq -c | sed "s/[ \t]*//" | sed 's/ /,/g' | sed 's/9606.//g' > 9606.protein.links.v10.5.Freq.Prot1.score400.txt


# Get only binding links ( for prot1 )
## All scores included
grep binding 9606.protein.actions.v10.5.txt | cut -f1  | sort | uniq -c | sed "s/[ \t]*//" | sed 's/ /,/g' | sed 's/9606.//g' > 9606.binding.links.v10.5.Freq.Prot1_allScores.txt

## Only score > 400 included
grep binding 9606.protein.actions.v10.5.txt |  awk ' { if ( $6 > 399 ) { print $0 } } ' | cut -f1  | sort | uniq -c | sed "s/[ \t]*//" | sed 's/ /,/g' | sed 's/9606.//g' > 9606.binding.links.v10.5.Freq.Prot1_score400.txt

## Get binding acting proteins
## All scores
grep binding 9606.protein.actions.v10.5.txt |  awk ' { if ( $5 == "t" ) { print $0 } } ' | cut -f1  | sort | uniq -c | sed "s/[ \t]*//" | sed 's/ /,/g' | sed 's/9606.//g' > 9606.binding.Acting.links.v10.5.Freq.Prot1_allScores.txt

## only score >400
grep binding 9606.protein.actions.v10.5.txt |  awk ' { if ( $5 == "t" ) { print $0 } } ' |  awk ' { if ( $6 > 399 ) { print $0 } } ' | cut -f1 | sort | uniq -c | sed "s/[ \t]*//" | sed 's/ /,/g' | sed 's/9606.//g' > 9606.binding.Acting.Links.v10.5.Freq.Prot1_score400.txt

## Get all acting interactions (binding, inhibition etc..)
