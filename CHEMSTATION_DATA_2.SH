
#Tuesday May 10
#extracting data from chemstation
#modified due to naming issue on sequence, usd fileinfo.txt to get file names insted
#usage locate this file in folder with al .D chemstation folders .D and run 
#bash CHEMSTATION DATA_2.SH


for sample_dir in */ ; do
    
	cd $sample_dir
    
   
    str_sid=$(sed '9!d' fileinfo.txt) 
	sid=$(awk '{ sub(/.*:/, ""); print }' <<< "$str_sid" | rev | cut -c2- | rev)

	str_date=$(sed '5!d' rteres.txt)
	date_time=$(echo $str_date | cut -c10-34)
	
	cut -d "," -f 2,3,9 RESULTS.CSV | sed '/^$/d' | sed '1,6d' | 
	awk -F, '{$1=(sample_id) FS $1;}1' sample_id="$sid" OFS=, | 
	awk -F, '{$(NF+1)=date;}1' date="$date_time" OFS=, | sed '1d' > sumtable.csv    
    
    cd ..
    
    cat $sample_dir/sumtable.csv  >> gcms_data_2.csv
    
done