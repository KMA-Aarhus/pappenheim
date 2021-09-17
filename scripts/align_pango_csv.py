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
		output[headers.index("lineage")] = final_lineage
	with open(outfile, 'w', newline='\n') as csvfile:
		lineagewriter = csv.writer(csvfile, delimiter=',')
		lineagewriter.writerow(headers)
		lineagewriter.writerow(output)
	csvfile.close()




def find_lineage(pangolearn_lineage, pango_usher_lineage, alias_key_json):
	with open(alias_key_json) as json_file:
		alias_keys = json.load(json_file)
	#inv_alias_keys = inv_map = {v: k for k, v in alias_keys.items()}
	json_file.close()
	print(alias_keys)
	if pango_usher_lineage.count(".") == 1:
		pango_usher_parental = pango_usher_lineage.split(".")[0]
		pango_usher_lineage = alias_keys[pango_usher_parental]
	if pangolearn_lineage.count(".") == 1:
		pangolearn_parental = pangolearn_lineage.split(".")[0]
		pangolearn_lineage = alias_keys[pangolearn_parental]
	if pangolearn_lineage == pango_usher_lineage:
		return pangolearn_lineage
	else:    
		return "UNKNOWN - pangoLEARN and USHER disagrees"

#scriptdir = "/Users/admin1/Documents/KMA_AUH/SARS-CoV-2/pappenheim/scripts/"
align_csv(snakemake.input[0], snakemake.input[1], snakemake.output[0])
