# NOTE: check has_both_modifications before running!
# USAGE: 
# Y: ruby common_peps_and_top_peaks.rb ../data/42-45_aspirin/elute_csvs/F003903_SCP_elute_Y_with_pipes.csv ../data/42-45_aspirin/elute_csvs/F003903.dat ../data/42-45_aspirin/elute_csvs/F003900_KAT_elute_Y_with_pipes.csv ../data/42-45_aspirin/elute_csvs/F003900.dat 100.0 ../results/F003903_SCP_elute_Y_transitions.csv ../results/F003900_KAT_elute_Y_transitions.csv
# S: ruby common_peps_and_top_peaks.rb ../data/42-45_aspirin/elute_csvs/F003894_SCP_elute_S_with_pipes.csv ../data/42-45_aspirin/elute_csvs/F003894.dat ../data/42-45_aspirin/elute_csvs/F003892_KAT_elute_S_with_pipes.csv ../data/42-45_aspirin/elute_csvs/F003892.dat 100.0 ../results/F003894_SCP_elute_S_transitions.csv ../results/F003892_KAT_elute_S_transitions.csv
# T: ruby common_peps_and_top_peaks.rb ../data/42-45_aspirin/elute_csvs/F003899_SCP_elute_T_with_pipes.csv ../data/42-45_aspirin/elute_csvs/F003899.dat ../data/42-45_aspirin/elute_csvs/F003896_KAT_elute_T_with_pipes.csv ../data/42-45_aspirin/elute_csvs/F003896.dat 100.0 ../results/F003899_SCP_elute_T_transitions.csv ../results/F003896_KAT_elute_T_transitions.csv

# Two subjects KAT and SCP.
# For each one of them, for each aa (Y,S,T) pick the peptides with both aspirin +42da and +45da (42 and 45 is not going to be found in the same peptide match, but in different matches).
# F003900_KAT_elute_Y.csv
# F003903_SCP_elute_Y.csv
# F003892_KAT_elute_S.csv
# F003894_SCP_elute_S.csv
# F003896_KAT_elute_T.csv
# F003899_SCP_elute_T.csv
# For each aa, get the common peptides from both subjects.
# For each common peptide get the top peaks for both subjects.

require 'rubygems'
require 'fastercsv'
require 'pep'
require 'mascot/dat'
require 'mascot_hits_csv_parser'

sample1_file = ARGV[0]
dat1_file = ARGV[1]
sample2_file = ARGV[2]
dat2_file = ARGV[3]
cutoff = ARGV[4]
transitions1_file = ARGV[5]
transitions2_file = ARGV[6]
modification = []
mod_positions = []
mod_positions_str = nil
spectrum = {}

# pick the peptides with both aspirin +42da and +45da modifications, from sample1 
sample1_peptides = {}
csvp = MascotHitsCSVParser.open(sample1_file, cutoff)
csvp.each_peptide do |peptide|
	if csvp.has_both_modifications(peptide)
		highest_scored_hit = csvp.highest_scored_hit_for_pep(peptide)
		sample1_peptides[peptide] = highest_scored_hit	
	end
end

# pick the peptides with both aspirin +42da and +45da modifications, from sample2 
sample2_peptides = {}
csvp2 = MascotHitsCSVParser.open(sample2_file, cutoff)
csvp2.each_peptide do |peptide|
	if csvp2.has_both_modifications(peptide)
		highest_scored_hit = csvp2.highest_scored_hit_for_pep(peptide)
		sample2_peptides[peptide] = highest_scored_hit
	end
end

# get the common peptides between the two created lists
common_peptides = {}
sample2_peptides.each do |peptide, hit|
	if sample1_peptides.has_key?(peptide)
		common_peptides[peptide] = hit
	end
end

# get the top peaks of the common peptides for sample1
FasterCSV.open(transitions1_file,'w') do |csv|
	csv << %w{QUERY_NUM PROT_ACCESSION PROT_DESCRIPTION PEPTIDE MODIFICATION MOD_PEPTIDE SPECTRUM_ID CHARGE RET_TIME_MINS HIGHEST_PEP_SCORE PEP_EXPECTANCY PEP_EXP_MR PEP_CALC_MR PEP_DELTA MS2_MZ MS2_INT y_INT_RANK y_INDEX y_ION_SEQ y-1 y y+2}

	sample1_peptides.each do |peptide, hit|
		if common_peptides.has_key?(peptide)
			peptide = hit.pep_seq
			prot_acc = hit.prot_acc
			prot_desc = hit.prot_desc
			modification.clear
			# Acetyl (Y), Acetyl:2H(3) (Y)
			hit.pep_var_mod.scan(/\d*Acetyl*.+\s\((\w)\)/).each do |i|
				modification << $1
			end
			#get modification positions from string
			mod_positions_str = hit.pep_var_mod_pos.split('.')[1].split('')
			mod_positions.clear
			mod_positions_str.each_index do |i|
				if !mod_positions_str[i].include?('0') && modification.include?(peptide.split('')[i])
					mod_positions << i + 1 #1-based
				end
			end
			pep_query = hit.pep_query
			pep_exp_mr = hit.pep_exp_mr
			pep_calc_mr = hit.pep_calc_mr
			pep_delta = hit.pep_delta
			pep_expect = hit.pep_expect
			pep_score = hit.pep_score	

			# take the ions table from dat file
			dat = Mascot::DAT.open(dat1_file, true)
			spectrum = dat.query(pep_query)
			title = spectrum.title
			charge = spectrum.charge
			rtinseconds = spectrum.rtinseconds
			ions1 = spectrum.peaks
			mzs = spectrum.mz
			intensities = spectrum.intensity

			pep = Pep.new(peptide, mzs, intensities, mod_positions)
			assigned_yions = pep.assigned_yions
			
			# Pick up the 5 most significant transitions from MS2 spectrum for each peptide & print transitions file
			ranked_idx = pep.ranked_yions_intensities_idx
			ranked_idx.each_with_index do |i,ii|
				if( assigned_yions[i] && !assigned_yions[i][0].nil? && 1
					assigned_yions[i][0] > 0 && ii < 5)
					yidx = assigned_yions[i][0] - 1
					csv << [
						pep_query,
						prot_acc,
						prot_desc,
						peptide,
						modification,
						pep.to_s,
						title,
						charge,
						rtinseconds.to_f / 60 ,
						pep_score,
						pep_expect,
						pep_exp_mr, #parent_mass
						pep_calc_mr, #calc_mass
						pep_delta,
						mzs[i],
						intensities[i],
						ii, # the mz intensity rank
						assigned_yions[i][0],
						pep.yions[yidx][0], # y ion sequence
						assigned_yions[i][1] - H ,
						assigned_yions[i][1],
						assigned_yions[i][1] + (H * 2)
					]
				end
			end
		end
	end
end

# get the top peaks of the common peptides for sample2
FasterCSV.open(transitions2_file,'w') do |csv|
	csv << %w{QUERY_NUM PROT_ACCESSION PROT_DESCRIPTION PEPTIDE MODIFICATION MOD_PEPTIDE SPECTRUM_ID CHARGE RET_TIME_MINS HIGHEST_PEP_SCORE PEP_EXPECTANCY PEP_EXP_MR PEP_CALC_MR PEP_DELTA MS2_MZ MS2_INT y_INT_RANK y_INDEX y_ION_SEQ y-1 y y+2}

	sample2_peptides.each do |peptide, hit|
		if common_peptides.has_key?(peptide)
			peptide = hit.pep_seq
			prot_acc = hit.prot_acc
			prot_desc = hit.prot_desc
			modification.clear
			# Acetyl (Y), Acetyl:2H(3) (Y)
			hit.pep_var_mod.scan(/Acetyl*.+\s\((\w)\)/).each do |i|
				modification << $1
			end
			#get modification positions from string
			mod_positions_str = hit.pep_var_mod_pos.split('.')[1].split('')
			mod_positions.clear
			mod_positions_str.each_index do |i|
				if !mod_positions_str[i].include?('0') && modification.include?(peptide.split('')[i])
					mod_positions << i + 1 #1-based
				end
			end
			pep_query = hit.pep_query
			pep_exp_mr = hit.pep_exp_mr
			pep_calc_mr = hit.pep_calc_mr
			pep_delta = hit.pep_delta
			pep_expect = hit.pep_expect
			pep_score = hit.pep_score	

			# take the ions table from dat file
			dat = Mascot::DAT.open(dat2_file, true)
			spectrum = dat.query(pep_query)
			title = spectrum.title
			charge = spectrum.charge
			rtinseconds = spectrum.rtinseconds
			ions1 = spectrum.peaks
			mzs = spectrum.mz
			intensities = spectrum.intensity

			pep = Pep.new(peptide, mzs, intensities, mod_positions)
			assigned_yions = pep.assigned_yions
			
			# Pick up the 5 most significant transitions from MS2 spectrum for each peptide & print transitions file
			ranked_idx = pep.ranked_yions_intensities_idx
			ranked_idx.each_with_index do |i,ii|
				if( assigned_yions[i] && !assigned_yions[i][0].nil? && 1
					assigned_yions[i][0] > 0 && ii < 5)
					yidx = assigned_yions[i][0] - 1
					csv << [
						pep_query,
						prot_acc,
						prot_desc,
						peptide,
						modification,
						pep.to_s,
						title,
						charge,
						rtinseconds.to_f / 60 ,
						pep_score,
						pep_expect,
						pep_exp_mr, #parent_mass
						pep_calc_mr, #calc_mass
						pep_delta,
						mzs[i],
						intensities[i],
						ii, # the mz intensity rank
						assigned_yions[i][0],
						pep.yions[yidx][0], # y ion sequence
						assigned_yions[i][1] - H ,
						assigned_yions[i][1],
						assigned_yions[i][1] + (H * 2)
					]
				end
			end
		end
	end
end
