# USAGE: ruby common_peps.rb ../results/F001808_transitions.csv ../results/F001809_transitions.csv ../results/common_peps.csv

# find common peptides from the 2 lists and print all matches per peptide (from both experiments)

require 'rubygems'
require 'fastercsv'

exp1_csvfile = ARGV[0]
exp2_csvfile = ARGV[1]
outfile = ARGV[2]

exp1_peps = Hash.new { |h,k| h[k] = [] }
FasterCSV.foreach(exp1_csvfile) do |row|
	if !row[1].include? "PROT_ACCESSION"
		row_updated = row << 1 
		exp1_peps[[row[1],row[3]]] << row_updated
	end
end

common_peptides = Hash.new { |h,k| h[k] = [] }
FasterCSV.foreach(exp2_csvfile) do |row|
	if exp1_peps.has_key?([row[1],row[3]])
		if !common_peptides.has_key?([row[1],row[3]])
			exp1_peps[[row[1],row[3]]].each do |value|
				common_peptides[[row[1],row[3]]] << value
			end
		end
		row_updated = row << 2 
		common_peptides[[row[1],row[3]]] << row_updated
	end
end

FasterCSV.open(outfile,'w') do |csv|
	csv << %w{QUERY_NUM PROT_ACCESSION PROT_DESCRIPTION PEPTIDE MODIFICATION MOD_PEPTIDE SPECTRUM_ID CHARGE RET_TIME_MINS HIGHEST_SCORE CUTOFF PARENT_MZ PEP_CALC_MASS PEP_DELTA MS2_MZ MS2_INT y_INT_RANK y_INDEX y_ION_SEQ y-1 y y+2 EXPERIMENT}

	common_peptides.each_key do |key|
		common_peptides[key].each do |value|
			csv << value
		end
	end
end

