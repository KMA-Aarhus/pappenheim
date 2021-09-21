import csv
import json

def align_csv(pangolearn_csv, usher_csv, outfile):
	with open(pangolearn_csv, newline='') as csvfile:
		pangolearn = list(csv.reader(csvfile, delimiter=','))
		headers = pangolearn[0]
		output = pangolearn[1]
		pangolearn_lineage = output[headers.index("lineage")]
	csvfile.close()	
	with open(usher_csv, newline='') as csvfile:
		usher = list(csv.reader(csvfile, delimiter=','))
		usher_lineage = usher[1][usher[0].index("lineage")]
	csvfile.close()	
	if pangolearn_lineage != usher_lineage:
		final_lineage = find_lineage(pangolearn_lineage, usher_lineage, snakemake.input[2])
		#final_lineage = find_lineage(pangolearn_lineage, usher_lineage, scriptdir+"alias.json")
		output[headers.index("lineage")] = final_lineage
	with open(outfile, 'w', newline='\n') as csvfile:
		lineagewriter = csv.writer(csvfile, delimiter=',')
		lineagewriter.writerow(headers)
		lineagewriter.writerow(output)
	csvfile.close()




def find_lineage(pangolearn_lineage, pango_usher_lineage, alias_key_json):
	with open(alias_key_json) as json_file:
		alias_keys = json.load(json_file)
	alias_keys_items = [v for k, v in alias_keys.items()]
	json_file.close()
	if pango_usher_lineage in pangolearn_lineage:
		return pango_usher_lineage
	if pango_usher_lineage.split(".")[0] in alias_keys.keys() and pango_usher_lineage not in alias_keys_items:
		print(alias_keys.items())
		pango_usher_parental = pango_usher_lineage.split(".")[0]
		pango_usher_lineage = alias_keys[pango_usher_parental]
		print(pango_usher_lineage)
		if isinstance(pango_usher_lineage, list):
			pango_usher_lineage = ",".join(pango_usher_lineage)
	if pangolearn_lineage.split(".")[0] in alias_keys.keys() and pangolearn_lineage not in alias_keys_items:
		pangolearn_parental = pangolearn_lineage.split(".")[0]
		pangolearn_lineage = alias_keys[pangolearn_parental]
		print(pangolearn_lineage)
		if isinstance(pangolearn_lineage, list):
			pangolearn_lineage = ",".join(pangolearn_lineage)
	if pangolearn_lineage == pango_usher_lineage:
		return pangolearn_lineage
	else:    
		return "UNKNOWN - pangoLEARN and USHER disagrees"

#scriptdir = "/Users/admin1/Documents/KMA_AUH/SARS-CoV-2/pappenheim/scripts/"
#align_csv(scriptdir+"20210920.0947_89333492_pangolearn.csv", scriptdir+"20210920.0947_89333492_usher.csv", scriptdir+"aligned.csv" )
align_csv(snakemake.input[0], snakemake.input[1], snakemake.output[0])
