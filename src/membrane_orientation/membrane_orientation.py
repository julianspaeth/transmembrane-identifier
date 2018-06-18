import argparse

from q_value_optimisation import process_results_file


def parse_arguments():
	"""
	Parses the input arguments from the terminal
	"""
	parser = argparse.ArgumentParser(description="Parse input for Transmembrane helix classifier.")
	parser.add_argument("-i", type=str, help="Path to results CSV of helix prediction", required=True)
	parser.add_argument("-pdb", type=str, help="Path to PDB file of the predicted protein", required=True)
	parser.add_argument("-o", type=str, help="Path to output JSON file", required=True)
	parser.add_argument("-plot", action='store_true', help="Optional. Will save figure of the Q-value distribution of the "
												"optimal angle.", required=False)
	parser.add_argument("-n_opti", type=int, help="Optional. Number of vectors to be tested in a 15° cone around "
												  "the preliminary orientation (1 vector ~ 10 seconds). Default: 20",
						required=False, default=20)

	return parser.parse_args()





if __name__ == "__main__":
	args = parse_arguments()
	results = process_results_file(args.i, args.pdb, args.o, args.plot, args.n_opti)
