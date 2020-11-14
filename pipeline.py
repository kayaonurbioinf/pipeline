import os
import sys
import platform
import subprocess
import multiprocessing
import subprocess
import tempfile
from greedyFAS import calcFAS
from greedyFAS.annoFAS import annoParserFAS
from pathlib import Path

sys.path.append(os.getcwd() + "/mhc_i/src")
sys.path.append(os.getcwd() + "/mhc_i/method/allele-info/allele_info")
sys.path.append(os.getcwd() + "/mhc_i/method/iedbtools-utilities/")

import argparse
import json
from predict_binding import Prediction
from seqpredictor import MHCBindingPredictions
from util import get_species, InputData


class pCommand:
	class info:
		info_text = "Use \"info predict\" for a clarification on how the command \"predict\" must be used.\n"\
                     "Use \"info ls\" to see a list of files in the current working directory.\n"\
                     "Use \"info cwd\" to see the current working directory."
		predict_text = "Takes 4 to 5 arguments: \"python pipeline.py predict(arg 1) [arg 2: file name] [arg 3: tool name] [arg 4: length of binding sites]\n\n"\
                        "Argument 2 - file name: The name of your input file (probably a fasta file). If the script cannot find your file,\n"\
                        "you may be excluding the file extension in your entry. Eg. \"my_input_file\" instead of \"my_input_file.fasta\".\n"\
                        "If you are having trouble getting the script to find your file, you may use the commands \"info ls\" to see a list of files in the current working directory\n"\
                        "or \"info cwd\" to see the path of the current working directory.\n\n"\
                        "Argument 3 - tool: IEDB offers epitope prediction tools for MHC I and MHC II.\n"\
                        "The only permitted entries for this argument are \"mhc_i\" for only MHC I epitope prediction,\n"\
                        "\"mhc_ii\" for only MHC II epitope prediction, or \"both\" for both of them. Both methods will use all alleles available in IEDB.\n\n"\
                        "Argument 4 - binding site length: Specification of the length of binding sites to look for.\n"\
                        "If you chose \"both\" tools to use, the length you pass for argument 4 will be used for MHC I.\n"\
                        "To look for binding sites of a single length, enter only the length. Eg: 15,\n"\
                        "To look for a range of binding sites, enter [lower boundary]-[upper boundary]. Eg: 12-15 -> 12,13,14,15.\n\n"\
                        "Argument 5 (If both MHC I and MHC II are to be used) - Lengths of binding sites to look for for MHC II\n\n"\
                        "Example 1: python pipeline.py predict my_input.fasta mhc_i 8-12\n"\
                        "Example 2: python pipeline.py predict my_input.fasta both 9-10 15-17\n"
		num_of_args = [1,2]
		arg2 = ["info", "predict", "ls", "cwd"]
	
		def print_info(): print(pCommand.info.info_text)
		def print_predict(): print(pCommand.info.predict_text)
		def print_ls():
			for f in [f for f in os.listdir(os.getcwd()) if os.path.isfile(os.path.join(os.getcwd(), f))]:
				print(f)
		def print_cwd():print(os.getcwd())
	
	
	class predict:
		num_of_args = [4,5]
		arg3 = ["mhc_i", "mhc_ii", "both"]
		length_dict={"MHC I": range(8,15),
					 "MHC II": range(11,31)}
	
	commands = ["info", "predict"]


class pFastaFile:
	def readFile(fname):
		slash = "\\" if platform.system() == 'Windows' else "/"
		filePath = os.getcwd() + slash + fname
		file = open (fname, 'r')
		sequences = {}
		seq_id = 0
		for line in file:
			if line[0] == '>':
				seq_id += 1
				sequences[seq_id] = {'seq': "", 'len': 0}
			else: 
				sequences[seq_id]['seq'] += line[:-1]
				sequences[seq_id]['len'] += len(line[:-1])
		return sequences


class pException:
	class commandNotFound(Exception):
		def __init__(self, command):
			self.command = command
			message = "Command {} does not exist."\
					  "Use \"info\" for a list of commands.".format(command)
			super().__init__(message)

	class tooFewArgs(Exception):
		def __init__(self, command, num_of_args):
			self.command = command
			self.num_of_args = num_of_args
			expected_num_of_args = pCommand.info.num_of_args if command == "info" else pCommand.predict.num_of_args
			message = "Arguments expected for \"{}\": {}, "\
				      "Arguments given: {}."\
				      "\nUse command \"info {}\" for more info.".format(
				      command, '-'.join([str(num) for num in expected_num_of_args]), num_of_args, command)
			super().__init__(message)

	class tooManyArgs(Exception):
		def __init__(self, command, num_of_args):
			self.command = command
			self.num_of_args = num_of_args
			message = "Arguments expected for \"{}\": {}, "\
				      "Arguments given: {}."\
				      "\nUse command \"info {}\" for more info.".format(
				      command, '-'.join([str(num) for num in pCommand.commands[command]["num_of_args"]]), num_of_args, command)
			super().__init__(message)

	class tooFewLengths(Exception):
		def __init__(self):
			message = "To use both MHC I and MHC II prediction tools, you need to pass 2 lengths, 1 length for each tool.\n"\
			          "Example: python pipeline.py predict my_input.fasta both 9-10 15-17              ||||9-10 -> MHC I | 15-17 -> MHC II"
			super().__init__(message)

	class tooManyLengths(Exception):
		def __init__(self, mhc):
			message = "2 lengths were passed while passing only 1 method: {}. Please pass only one length.".format(mhc)
			super().__init__(message)
	
	class lengthsDontMatch(Exception):
		def __init__(self, faulty_lengths):
			message = "\n"
			for length in faulty_lengths:
				message += "The range of lengths available for {} is {}-{}, length entered: {}.\n".format(
				length[0], min(pCommand.predict.length_dict[length[0]]), max(pCommand.predict.length_dict[length[0]]), length[1])
			super().__init__(message)
	
	class incorrectArgs(Exception):
		def __init__(self, incorrect_args):
			self.incorrect_args = incorrect_args
			message = "\n\n"
			for arg in incorrect_args:
				message_for_current_arg = "Faulty argument {}: {}:\n".format(arg["nr"], arg["arg"])
				message_for_current_arg += arg["note"] + '\n\n'
				message += message_for_current_arg
			super().__init__(message)



class pParse:
	def parseBSL(argBSL): #BSL = Binding Site Length
		"""
		Used for parsing args.
		Generates the list of binding site lengths to look for from the given argument.
		"15" -> [15]
		"12-16" -> [12,13,14,15,16]
		"""
		boundaries = [int(boundary) for boundary in argBSL.split('-')] # "12-16" -> [12,16]
		if len(boundaries) == 1: return boundaries
		binding_sites_length = [num for num in range(boundaries[0], boundaries[1]+1)] # [12,16] -> [12,13,14,15,16]
		return binding_sites_length


	def parseArgs():	
		cmd_args = sys.argv
		
		arg_dict = {}
		
		if cmd_args[1] not in pCommand.commands: raise pException.commandNotFound(cmd_args[1])
		elif cmd_args[1] == "predict":
			if len(cmd_args)-1 < min(pCommand.predict.num_of_args): raise pException.tooFewArgs(cmd_args[1], len(cmd_args)-1)
			elif len(cmd_args)-1 > max(pCommand.predict.num_of_args): raise pException.tooManyArgs(cmd_args[1], len(cmd_args)-1)
			incorrect_args = []
			if cmd_args[2] not in [f for f in os.listdir(os.getcwd()) if os.path.isfile(os.path.join(os.getcwd(), f))]:
				incorrect_args.append({"nr": 2, "arg": cmd_args[2]})
				incorrect_args[-1]["note"] = "File \"{}\" cannot be found in current working directory."\
				                             "Use command \"info ls\" to see a list of files in the current working directory.".format(cmd_args[2])
			if cmd_args[3] not in pCommand.predict.arg3:
				incorrect_args.append({"nr": 3, "arg": cmd_args[3]})
				incorrect_args[-1]["note"] = "The only accepted entries for third argument of \"{}\" are {}".format(
				cmd_args[1], ", ".join(pCommand.predict.arg3))
			try:
				pParse.parseBSL(cmd_args[4])
			except:
				incorrect_args.append({"nr": 4, "arg": cmd_args[4]})
				incorrect_args[-1]["note"] = "The fourth argument of \"{}\" (binding lengths) must be in form: \n"\
											 "[Lower boundary]-[Upper boundary] for a range of lengths or [Length] for a single length.\n"\
											 "(Make sure you aren't putting spaces between the numbers and the dash)". format(
											 cmd_args[1])
			if cmd_args[3] == "both" and len(cmd_args)-1 == max(pCommand.predict.num_of_args):
				try:
					pParse.parseBSL(cmd_args[5])
				except:
					incorrect_args.append({"nr": 5, "arg": cmd_args[5]})
					incorrect_args[-1]["note"] = "The fifthh argument of \"{}\" (binding lengths) must be in form: \n"\
												 "[Lower boundary]-[Upper boundary] for a range of lengths or [Length] for a single length.\n"\
												 "(Make sure you aren't putting spaces between the numbers and the dash)". format(
												 cmd_args[1])
			if len(incorrect_args) != 0: raise pException.incorrectArgs(incorrect_args)
			if cmd_args[3] == "both" and len(cmd_args)-1 != max(pCommand.predict.num_of_args):
				raise pException.tooFewLengths()
			elif cmd_args[3] != "both" and len(cmd_args)-1 != min(pCommand.predict.num_of_args):
				mhc = "MHC I" if cmd_args[3] == "mhc_i" else "MHC II"
				raise pException.tooManyLengths(mhc)
			if cmd_args[3] == "mhc_i":
				mhc_i_lengths = pParse.parseBSL(cmd_args[4])
				if mhc_i_lengths[0] < min(pCommand.predict.length_dict["MHC I"]) or mhc_i_lengths[-1] > max(pCommand.predict.length_dict["MHC I"]):
					raise pException.lengthsDontMatch([("MHC I", mhc_i_lengths)])
			elif cmd_args[3] == "mhc_ii":
				mhc_ii_lengths = pParse.parseBSL(cmd_args[4])
				if mhc_ii_lengths[0] < min(pCommand.predict.length_dict["MHC II"]) or mhc_ii_lengths[-1] > max(pCommand.predict.length_dict["MHC II"]):
					raise pException.lengthsDontMatch([("MHC II", mhc_ii_lengths)])
			elif cmd_args[3] == "both":
				faulty_lengths = []
				mhc_i_lengths = pParse.parseBSL(cmd_args[4])
				mhc_ii_lengths = pParse.parseBSL(cmd_args[5])
				if mhc_i_lengths[0] < min(pCommand.predict.length_dict["MHC I"]) or mhc_i_lengths[-1] > max(pCommand.predict.length_dict["MHC I"]):
					faulty_lengths.append(["MHC I", mhc_i_lengths])
				if mhc_ii_lengths[0] < min(pCommand.predict.length_dict["MHC II"]) or mhc_i_lengths[-1] > max(pCommand.predict.length_dict["MHC II"]):
					faulty_lengths.append(["MHC II", mhc_ii_lengths])
					raise pException.lengthsDontMatch(faulty_lengths)

			arg_dict["command"] = cmd_args[1]
			arg_dict["fname"] = cmd_args[2]
			arg_dict["method"] = cmd_args[3]
			arg_dict["bsl"] = pParse.parseBSL(cmd_args[4])
			if len(cmd_args)-1 == max(pCommand.predict.num_of_args): arg_dict["bsl2"] = pParse.parseBSL(cmd_args[5])
			arg_dict["json"] = True if input("Do you want to generate a json file?[y/n]: ") == 'y' else False
				
		elif cmd_args[1] == "info":
			if len(cmd_args)-1 > 2: raise pException.tooManyArgs(cmd_args[1], len(cmd_args)-1)
			elif len(cmd_args)-1 == 2 and cmd_args[2] not in pCommand.info.arg2:
				raise pException.commandNotFound(cmd_args[2])
			elif len(cmd_args)-1 == 2 and cmd_args[2] in pCommand.info.arg2:
				arg_dict["command"] = cmd_args[1]
				arg_dict["entry"] = cmd_args[2]
			elif len(cmd_args)-1 == 1:
				arg_dict["command"] = cmd_args[1]
				arg_dict["entry"] = None
		
		return arg_dict



class mhc_i:
	alleles = ["HLA-A*01:01", "HLA-A*02:01", "HLA-A*02:03", "HLA-A*02:06",
 		"HLA-A*03:01", "HLA-A*11:01", "HLA-A*23:01", "HLA-A*24:02",
		"HLA-A*26:01", "HLA-A*30:01", "HLA-A*30:02", "HLA-A*31:01",
		"HLA-A*32:01", "HLA-A*33:01", "HLA-A*68:01", "HLA-A*68:02",
		"HLA-B*07:02", "HLA-B*08:01", "HLA-B*15:01", "HLA-B*35:01",
		"HLA-B*40:01", "HLA-B*44:02", "HLA-B*44:03", "HLA-B*51:01",
		"HLA-B*53:01", "HLA-B*57:01", "HLA-B*58:01"]

	def predict(fname, binding_site_lengths):
		print("MHC I Prediction")
		combined_table_rows = []
		pPrediction = Prediction()
		proteins = pPrediction.read_protein(fname)
		species = [get_species(allele) for allele in mhc_i.alleles]
		progress = 1
		for length in binding_site_lengths:
			for allele in mhc_i.alleles:
				print("Allele: {}, Length: {} ({}/{})".format(length, allele, progress, len(mhc_i.alleles)*len(binding_site_lengths)))
				pInput = InputData(pPrediction.version, "netmhcpan_el", allele, '', length, proteins, species)
				predictor = MHCBindingPredictions(pInput)
				mhc_scores = predictor.predict(pInput.input_protein.as_amino_acid_text())
				use_cutoff = cutoff_value = None
				table_rows = pPrediction.format_binding(pInput, mhc_scores, "netmhcpan_el", use_cutoff, cutoff_value)
				combined_table_rows.extend(table_rows)
				progress += 1
		combined_table_rows.sort(key=lambda tup: tup[8])
		combined_table_rows.reverse()
		mhci_tsv = open(fname.split('.')[0]+"_mhci.tsv", 'w')
		write_file = '\n'.join(['\t'.join(map(str, row)) for row in combined_table_rows])
		header = '\t'.join(['allele','seq_num','start','end','length','peptide','core','icore',predictor.get_score_unit(),'rank'])+'\n'
		write_file = header + write_file
		mhci_tsv.write(write_file)
		mhci_tsv.close()
		return combined_table_rows
				

	
	def predict_process(my_allele, my_length, my_nr, lengths, fname, combined_table_rows, pPrediction, proteins, species):
		pInput = InputData(pPrediction.version, "netmhcpan_el", my_allele, '', my_length, proteins, species)
		predictor = MHCBindingPredictions(pInput)
		mhc_scores = predictor.predict(pInput.input_protein.as_amino_acid_text())
		print("Allele: {}, Length: {}    ({}/{})".format(my_allele, my_length, my_nr, len(mhc_i.alleles)*lengths))
		use_cutoff = cutoff_value = None
		table_rows = pPrediction.format_binding(pInput, mhc_scores, "netmhcpan_el", use_cutoff, cutoff_value)
		table_rows.sort(key=lambda tup: tup[6])
		combined_table_rows.extend(table_rows)
		return
		
	
	
	def predict_mp(fname, binding_site_lengths):
		combined_table_rows = []
		pPrediction = Prediction()
		proteins = pPrediction.read_protein(fname)
		species = [get_species(allele) for allele in mhc_i.alleles]
		
		processes = []
		for length in binding_site_lengths:
			for allele in mhc_i.alleles:
				process = multiprocessing.Process(target=mhc_i.predict_process, 
				args=(allele, length, mhc_i.alleles.index(allele)+1, len(binding_site_lengths), 
				fname, combined_table_rows, pPrediction, proteins, species))
				processes.append(process)
				process.start()

		combined_table_rows.sort(key=lambda tup: tup[2])
				
		return combined_table_rows


	def crop(result_table, fname):
		seqs = pFastaFile.readFile(fname)
		progress = 1
		cropped_table = []
		for row in result_table:
			#print("TO DO")
			if float(row[9]) < 1.0: cropped_table.append((row[1], seqs[int(row[1])]['len'], row[0], row[2], row[3], "mhc_i"))
		return cropped_table



class mhc_ii:
	alleles = ["HLA-DRB1*01:01", "HLA-DRB1*03:01", "HLA-DRB1*04:01", "HLA-DRB1*04:05",
		       "HLA-DRB1*07:01", "HLA-DRB1*08:02", "HLA-DRB1*09:01", "HLA-DRB1*11:01",
  		       "HLA-DRB1*12:01", "HLA-DRB1*13:02", "HLA-DRB1*15:01", "HLA-DRB3*01:01",
		       "HLA-DRB3*02:02", "HLA-DRB4*01:01", "HLA-DRB5*01:01",
		       "HLA-DQA1*05:01/DQB1*02:01", "HLA-DQA1*05:01/DQB1*03:01",
		       "HLA-DQA1*03:01/DQB1*03:02", "HLA-DQA1*04:01/DQB1*04:02",
		       "HLA-DQA1*01:01/DQB1*05:01", "HLA-DQA1*01:02/DQB1*06:02",
		       "HLA-DPA1*02:01/DPB1*01:01", "HLA-DPA1*01:03/DPB1*02:01", 
		       "HLA-DPA1*01:03/DPB1*04:01", "HLA-DPA1*03:01/DPB1*04:02", 
		       "HLA-DPA1*02:01/DPB1*05:01", "HLA-DPA1*02:01/DPB1*14:01"]

	
	def predict_process(my_nr, my_allele, my_length, total_progress, fname):
		result = ""
		with tempfile.TemporaryFile() as tempf:
			cmd = ['python', './mhc_ii/mhc_II_binding.py', 'IEDB_recommended', my_allele, fname, str(my_length)]
			proc = subprocess.Popen(cmd, stdout=tempf)
			proc.wait()
			tempf.seek(0)
			result = tempf.read().decode()
			f = open(str(my_nr)+".txt", 'w')
			f.close()
		print("Allele: {}, Length: {} ({}/{})".format(my_allele, my_length, my_nr, total_progress))


	def cartesian_product(lengths, fname):
		cartesian = []
		all_workers = len(lengths)*len(mhc_ii.alleles)
		worker_nr = 1
		for length in lengths:
			for allele in mhc_ii.alleles:
				cartesian.append((worker_nr,allele,length,all_workers,fname))
				worker_nr += 1
		return cartesian


	def predict_mp(fname, lengths):
		mp.freeze_support()
		manager = mp.Manager()
		processes = []
		mp.Pool(9).starmap(mhc_ii.predict_process, mhc_ii.cartesian_product(lengths, fname))
		
		print(results_dict)

		
	def predict(fname, lengths):
		print("MHC II Prediction")
		result = ""
		with tempfile.TemporaryFile() as tempf:
			progress = 1
			for allele in mhc_ii.alleles:
				for length in lengths:
					print("Allele: {}, Length: {} ({}/{})".format(allele, length, progress, len(mhc_ii.alleles)*len(lengths)))
					cmd = ['python', './mhc_ii/mhc_II_binding.py', 'IEDB_recommended', allele, fname, str(length)]
					proc = subprocess.Popen(cmd, stdout=tempf)
					proc.wait()
					tempf.seek(0)
					progress += 1
					result = tempf.read().decode()
		result = '\n'.join(list(dict.fromkeys(result.split('\n')))) # remove duplicate headers
		mhcii_tsv = open(fname.split('.')[0]+"_mhcii.tsv", 'w')
		mhcii_tsv.write(result)
		mhcii_tsv.close()
		return result


	def crop(result_table, fname):
		seqs = pFastaFile.readFile(fname)
		progress = 1
		cropped_table = []
		result_table = result_table.split('\n')
		result_table = [row.split('\t') for row in result_table]
		del result_table[1]
		for row in result_table:
			if len(row) > 1:
				if float(row[7]) < 1.0: 
					cropped_table.append((row[1], row[3], row[0], row[2], row[3], "mhc_ii"))
		return cropped_table



class FAS:
	def feature_dict_gen(feature_columns_table):
		feature_dict = {}
		for row in feature_columns_table:
			p_id, length, start, stop, tool = row[0], row[1], int(row[3]), int(row[4]), row[5]
			f_id = tool + '_'+ row[2]
			if p_id not in feature_dict:
				feature_dict[p_id] = {'length': str(length), tool: {}}
			if f_id not in feature_dict[p_id][tool]:
				feature_dict[p_id][tool][f_id] = {'evalue': 'NA', 'instance': []}
			feature_dict[p_id][tool][f_id]['instance'].append([start, stop, 'NA'])
		return feature_dict

	def json_gen(fname, feature_columns_table):
		feature_dict = FAS.feature_dict_gen(feature_columns_table)
		count = annoParserFAS.count_features(feature_dict)
		out_dict = {'feature': feature_dict, 'clan': {}, 'count': count}
		jsonOut = json.dumps(out_dict, ensure_ascii=False)
		out_file_path = os.getcwd()+"/my_annotations/"+fname.split('.')[0] + ".json"
		f = open(out_file_path, 'w')
		f.write(jsonOut)
		f.close()
	
	def feature_types_file_gen(tools):
		my_featuretypes_file = open(os.getcwd()+"/my_featuretypes_file", 'w')
		if tools == "mhc_i": my_featuretypes_file.write("#linearized\n#normal\nmhc_i")
		elif tools == "mhc_ii": my_featuretypes_file.write("#linearized\n#normal\nmhc_ii")
		elif tools == "both": my_featuretypes_file.write("#linearized\n#normal\nmhc_i\nmhc_ii")
		my_featuretypes_file.close()

	def calculateFAS(fname):
		args = argparse.Namespace(annotation_dir='my_annotations/', bidirectional=True, cpus=0,
			         	 domain=True, eFeature=0.001, eFlps=1e-07, eInstance=0.01,
			         	 extra_annotation=None, featuretypes='my_featuretypes_file', force=False,
			         	 max_cardinality=500, max_overlap=0, max_overlap_percentage=0.4,
			         	 org='euk', out_dir='fas_out/', out_name='pipe_output_file',
			         	 pairwise=None, phyloprofile=None, priority_mode=True,
			         	 priority_threshold=30, query=fname, query_id=None, raw=False,
			         	 ref_2=None, ref_proteome=None, score_weights=[0.7,0.0,0.3],
			         	 seed=fname, seed_id=None, silent=False, timelimit=3600,
			         	 toolPath = '', weight_constraints=None, weight_correction='loge')
		toolpath = calcFAS.__file__.replace('calcFAS.py', 'pathconfig.txt')
		print(args)
		print('\n\n', toolpath)
		calcFAS.fas(args, toolpath)



def main():
	arg_dict = pParse.parseArgs()
	
	feature_columns_table = []
	
	if arg_dict["command"] == "predict":
		
		if arg_dict["method"] == "mhc_i":
			feature_columns_table = mhc_i.crop(mhc_i.predict(arg_dict["fname"], arg_dict["bsl"]), arg_dict["fname"])
			FAS.json_gen(arg_dict['fname'], feature_columns_table)
		elif arg_dict["method"] == "mhc_ii":
			feature_columns_table = mhc_ii.crop(mhc_ii.predict(arg_dict["fname"], arg_dict["bsl"]), arg_dict["fname"])
			FAS.json_gen(arg_dict['fname'], feature_columns_table)
		elif arg_dict["method"] == "both":
			feature_columns_table_mhc_i = mhc_i.crop(mhc_i.predict(arg_dict["fname"], arg_dict["bsl"]))
			feature_columns_table_mhc_ii = mhc_ii.crop(mhc_ii.predict(arg_dict["fname"], arg_dict["bsl2"]))
			feature_columns_table = feature_columns_table_mhc_i + feature_columns_table_mhc_ii
			FAS.json_gen(arg_dict['fname'], feature_columns_table)
		pathconfigfile = os.path.realpath(__file__).replace('calcFAS.py', 'pathconfig.txt')
		with open(pathconfigfile) as f:
			toolpath = f.readline().strip()
		FAS.feature_types_file_gen(arg_dict["method"])
		FAS.calculateFAS(arg_dict['fname'])
		
	
	
	elif arg_dict["command"] == "info":
		if arg_dict["entry"] == "info" or arg_dict["entry"] == None: pCommand.info.print_info()
		elif arg_dict["entry"] == "predict": pCommand.info.print_predict()
		elif arg_dict["entry"] == "ls": pCommand.info.print_ls()
		elif arg_dict["entry"] == "cwd": pCommand.info.print_cwd()



main()
