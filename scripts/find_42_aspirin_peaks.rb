# USAGE: ruby find_42_aspirin_peaks.rb ../data/42-45_aspirin/F003886_with_pipes.csv ../data/42-45_aspirin/F003886.dat ../data/42-45_aspirin/lysine_acetylation_list_11142013_F003886.csv 100.0 ../results/F003886_transitions.csv

# ruby find_42_aspirin_peaks.rb ../data/42-45_aspirin/F003887_with_pipes.csv ../data/42-45_aspirin/F003887.dat ../data/42-45_aspirin/lysine_acetylation_list_11142013_F003887.csv 100.0 ../results/F003887_transitions.csv

require 'rubygems'
require 'fastercsv'
require 'pep'
require 'mascot/dat'
require 'mascot_hits_csv_parser'

csvfile = ARGV[0]
datfile = ARGV[1]
common_peps_file = ARGV[2]
cutoff = ARGV[3]
outfile = ARGV[4]
fieldnames = []
modification = []
mod_positions = []
mod_positions_str = nil
spectrum = {}

common_peptides = {}
FasterCSV.foreach(common_peps_file) do |row|
	common_peptides[row[4].to_s] = nil
end

FasterCSV.open(outfile,'w') do |csv|
	csv << %w{QUERY_NUM PROT_ACCESSION PROT_DESCRIPTION PEPTIDE MODIFICATION MOD_PEPTIDE SPECTRUM_ID CHARGE RET_TIME_MINS HIGHEST_SCORE CUTOFF PARENT_MZ PEP_CALC_MASS PEP_DELTA MS2_MZ MS2_INT y_INT_RANK y_INDEX y_ION_SEQ y-1 y y+2}

	csvp = MascotHitsCSVParser.open(csvfile, cutoff)
	csvp.each_peptide do |peptide|
		if common_peptides.has_key?(peptide)
			if petide.eq?('ICLDLQAPLYKK')
				highest_scored_hit = csvp.highest_scored_hit_for_pep(peptide)
				puts highest_scored_hit.pep_var_mod
				peptide = highest_scored_hit.pep_seq
				prot_acc = highest_scored_hit.prot_acc
				prot_desc = highest_scored_hit.prot_desc

				modification.clear
				highest_scored_hit.pep_var_mod.scan(/Acetyl\s\((\w)\)/).each do |i| # Acetyl (K)
					modification << $1
				end
				#get modification positions from string
				mod_positions_str = highest_scored_hit.pep_var_mod_pos.split('.')[1].split('')
				mod_positions.clear
				mod_positions_str.each_index do |i|
					if !mod_positions_str[i].include?('0') && modification.include?(peptide.split('')[i])
						mod_positions << i + 1 #1-based
					end
				end

				pep_query = highest_scored_hit.pep_query
				pep_exp_mr = highest_scored_hit.pep_exp_mr
				pep_calc_mr = highest_scored_hit.pep_calc_mr
				pep_delta = highest_scored_hit.pep_delta
				pep_expect = highest_scored_hit.pep_expect
				pep_score = highest_scored_hit.pep_score	

				# take the ions table from dat file
				dat = Mascot::DAT.open(datfile, true)
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
end

