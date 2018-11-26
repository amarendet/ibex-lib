//============================================================================
//                                  I B E X
//
//                               ************
//                                 IbexSolve
//                               ************
//
// Author      : Gilles Chabert
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Last Update : Oct 01, 2017
//============================================================================

#include "ibex.h"
#include "args.hxx"

#include <sstream>

#ifndef _IBEX_WITH_SOLVER_
#error "You need to install the IbexSolve plugin (--with-solver)."
#endif

using namespace std;
using namespace ibex;

struct VectorHelper {
	ibex::Vector value = ibex::Vector(1, 0.0);
};

struct IntervalVectorHelper {
	ibex::IntervalVector value = ibex::IntervalVector(1, 0.0);
};

namespace args {
template<>
struct ValueReader<VectorHelper> {
	bool operator()(const std::string& name, const std::string& value, VectorHelper& dest) {
		std::vector<double> numbers;
		if (value[0] != '(' || value[value.size() - 1] != ')') {
			throw args::ParseError(value + " is not a valid vector for " + name);
		}
		size_t pos = value.find(';', 1);
		size_t last_pos = 1;
		while (pos != string::npos) {
			try {
				numbers.emplace_back(std::stod(value.substr(last_pos, pos - last_pos)));
			} catch (std::invalid_argument& e) {
				throw args::ParseError(value + " is not a valid vector for " + name + ": " + e.what());
			}
			last_pos = pos + 1;
			pos = value.find_first_of(";)", pos + 1);
		}
		dest.value = ibex::Vector(numbers.size(), 0.0);
		for (size_t i = 0; i < numbers.size(); ++i) {
			dest.value[i] = numbers[i];
		}
		return true;
	}
};

template<>
struct ValueReader<IntervalVectorHelper> {
	bool operator()(const std::string& name, const std::string& value, IntervalVectorHelper& dest) {
		std::vector<ibex::Interval> numbers;
		if (value[0] != '(' || value[value.size() - 1] != ')') {
			throw args::ParseError(value + " is not a valid vector for " + name);
		}
		size_t pos = value.find(';', 1);
		size_t last_pos = 1;
		while (pos != string::npos) {
			try {
				if(value[last_pos] == '[') {
					last_pos += 1;
					size_t pos_virgule = value.find_first_of(",");
					double lb = std::stod(value.substr(last_pos, pos_virgule - last_pos));
					size_t pos_end = value.find_first_of("]");
					double ub = std::stod(value.substr(pos_virgule+1, pos_end - pos_virgule));
					numbers.emplace_back(ibex::Interval(lb, ub));
				} else {
					numbers.emplace_back(std::stod(value.substr(last_pos, pos - last_pos)));
				}
			} catch (std::invalid_argument& e) {
				throw args::ParseError(value + " is not a valid interval vector for " + name + ": " + e.what());
			}
			last_pos = pos + 1;
			pos = value.find_first_of(";)", pos + 1);
		}
		dest.value = ibex::IntervalVector(numbers.size(), 0.0);
		for (size_t i = 0; i < numbers.size(); ++i) {
			dest.value[i] = numbers[i];
		}
		return true;
	}
};
}

int main(int argc, char** argv) {

	stringstream _random_seed, _eps_x_min, _eps_x_max;
	_random_seed << "Random seed (useful for reproducibility). Default value is " << DefaultSolver::default_random_seed << ".";
	_eps_x_min << "Minimal width of output boxes. This is a criterion to _stop_ bisection: a "
			"non-validated box will not be larger than 'eps-min'. Default value is 1e" << round(::log10(DefaultSolver::default_eps_x_min)) << ".";
	_eps_x_max << "Maximal width of output boxes. This is a criterion to _force_ bisection: a "
			"validated box will not be larger than 'eps-max' (unless there is no equality and it is fully inside inequalities)."
			" Default value is +oo (none)";

	args::ArgumentParser parser("********* IbexSolve (defaultsolver) *********.", "Solve a Minibex file.");
	args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
	args::ValueFlag<double> eps_x_min(parser, "float", _eps_x_min.str(), {'e', "eps-min"});
	args::ValueFlag<double> eps_x_max(parser, "float", _eps_x_max.str(), {'E', "eps-max"});
	args::ValueFlag<double> timeout(parser, "float", "Timeout (time in seconds). Default value is +oo (none).", {'t', "timeout"});
	args::ValueFlag<string> input_file(parser, "filename", "Manifold input file. The file contains a "
			"(intermediate) description of the manifold with boxes in the MNF (binary) format.", {'i',"input"});
	args::ValueFlag<string> output_file(parser, "filename", "Manifold output file. The file will contain the "
			"description of the manifold with boxes in the MNF (binary) format.", {'o',"output"});
	args::Flag format(parser, "format", "Show the output text format", {"format"});
	args::Flag bfs(parser, "bfs", "Perform breadth-first search (instead of depth-first search, by default)", {"bfs"});
	args::Flag txt(parser, "txt", "Write the output manifold in a easy-to-parse text file. See --format", {"txt"});
	args::Flag trace(parser, "trace", "Activate trace. \"Solutions\" (output boxes) are displayed as and when they are found.", {"trace"});
	args::ValueFlag<string> boundary_test_arg(parser, "true|full-rank|half-ball|false", "Boundary test strength. Possible values are:\n"
			"\t\t* true:\talways satisfied. Set by default for under constrained problems (0<m<n).\n"
			"\t\t* full-rank:\tthe gradients of all constraints (equalities and potentially activated inequalities) must be linearly independent.\n"
			"\t\t* half-ball:\t(**not implemented yet**) the intersection of the box and the solution set is homeomorphic to a half-ball of R^n\n"
	        "\t\t* false: never satisfied. Set by default if m=0 or m=n (inequalities only/square systems)",
			{"boundary"});
	args::Flag sols(parser, "sols", "Display the \"solutions\" (output boxes) on the standard output.", {'s',"sols"});
	args::ValueFlag<double> random_seed(parser, "float", _random_seed.str(), {"random-seed"});
	args::Flag quiet(parser, "quiet", "Print no report on the standard output.",{'q',"quiet"});
	args::ValueFlag<string> forced_params(parser, "vars","Force some variables to be parameters in the parametric proofs.",{"forced-params"});

	args::Group group_path(parser, "This group specify the path for path finding:", args::Group::Validators::AllOrNone);
	args::ValueFlag<std::string> start_point(group_path, "vector", "Starting point for path finding (constant name)", { "start" });
	args::ValueFlag<std::string> goal_point(group_path, "vector", "End point for path finding (constant name)", { "goal" });
	args::Group group(parser, "This group specify the algorithm to use for path finding (default: Dijkstra):",
			args::Group::Validators::AtMostOne);
	args::Flag path_finding_dijkstra(group, "dijkstra", "Use Dijkstra's algorithm", { "dijkstra" });
	args::Flag path_finding_astar_distance(group, "astar_distance", "Use A* with distance heuristic", { "astar-dist" });

	args::Positional<std::string> filename(parser, "filename", "The name of the MINIBEX file.");
	
	try
	{
		parser.ParseCLI(argc, argv);
	}
	catch (args::Help&)
	{
		std::cout << parser;
		return 0;
	}
	catch (args::ParseError& e)
	{
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	}
	catch (args::ValidationError& e)
	{
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	}

	if (format) {
		cout << Manifold::format() << endl;
		exit(0);
	}

	if (filename.Get()=="") {
		ibex_error("no input file (try ibexsolve --help)");
		exit(1);
	}

	std::vector<std::string> accepted_options = { "--start", "--goal" };
	MinibexOptionsParser minibexParser(accepted_options);
	minibexParser.parse(filename.Get());
	vector<string> unsupported_options = minibexParser.unsupported_options();
	for (const string& s : unsupported_options) {
		ibex::ibex_warning("Unsupported option in minibex file: " + s);
	}

	vector<string> argv_vec = minibexParser.as_argv_list();
	int new_argc = argc + argv_vec.size();
	const char* new_argv[new_argc];
	new_argv[0] = argv[0];
	int arg_index = 1;
	for (const string& name : argv_vec) {
		new_argv[arg_index] = name.c_str();
		++arg_index;
	}
	for(int j = 1; j < argc; ++j) {
		new_argv[arg_index] = argv[j];
		++arg_index;
	}

	// second pass
	try {
		parser.ParseCLI(new_argc, new_argv);
	} catch (args::Help&) {
		std::cout << parser;
		return 0;
	} catch (args::ParseError& e) {
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	} catch (args::ValidationError& e) {
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	}

	try {

		// Load a system of equations
		System sys(filename.Get().c_str());

		string output_manifold_file; // manifold output file
		bool overwitten=false;       // is it overwritten?
		string manifold_copy;

		if (!quiet) {
			cout << endl << "***************************** setup *****************************" << endl;
			cout << "  file loaded:\t\t" << filename.Get() << endl;

			if (eps_x_min)
				cout << "  eps-x:\t\t" << eps_x_min.Get() << "\t(precision on variables domain)" << endl;

			// Fix the random seed for reproducibility.
			if (random_seed)
				cout << "  random seed:\t\t" << random_seed.Get() << endl;

			if (bfs)
				cout << "  bfs:\t\t\tON" << endl;
		}

		if (output_file) {
			output_manifold_file = output_file.Get();
		} else {
			// got from stackoverflow.com:
			string::size_type const p(filename.Get().find_last_of('.'));
			// filename without extension
			string filename_no_ext=filename.Get().substr(0, p);
			stringstream ss;
			ss << filename_no_ext << ".mnf";
			output_manifold_file=ss.str();

			ifstream file;
			file.open(output_manifold_file.c_str(), ios::in); // to check if it exists

			if (file.is_open()) {
				overwitten = true;
				stringstream ss;
				ss << output_manifold_file << "~";
				manifold_copy=ss.str();
				// got from stackoverflow.com:
				ofstream dest(manifold_copy, ios::binary);

			    istreambuf_iterator<char> begin_source(file);
			    istreambuf_iterator<char> end_source;
			    ostreambuf_iterator<char> begin_dest(dest);
			    copy(begin_source, end_source, begin_dest);
			}
			file.close();
		}

		if (!quiet) {
			cout << "  output file:\t\t" << output_manifold_file << "\n";
			if (txt)
				cout << "  output format:\tTXT" << endl;
		}

		Solver* solver = nullptr;

		ibex::CellBuffer* buffer = nullptr;
		bool pathFinding = false;
		if (start_point && goal_point) {
			if(!quiet) {
				cout << "  path finding: ON\n";
			}
			ibex::IntervalVector start_vector = sys.csts.find(start_point.Get())->second.get().get_vector_value();
			ibex::IntervalVector goal_vector = sys.csts.find(goal_point.Get())->second.get().get_vector_value();
			/*double c1[] = {-3,-3,5,5};
			double c2[] = {3,-3,5,5};
			ibex::Vector start_vector(4,c1);
			ibex::Vector goal_vector(4,c2);*/
			ibex::CellBufferNeighborhood::Heuristic heuristic = ibex::CellBufferNeighborhood::Heuristic::DIJKSTRA;
			if (path_finding_dijkstra) {
				heuristic = ibex::CellBufferNeighborhood::Heuristic::DIJKSTRA;
			} else if (path_finding_astar_distance) {
				heuristic = ibex::CellBufferNeighborhood::Heuristic::A_STAR_DISTANCE;
			}
			buffer = new ibex::CellBufferNeighborhood(start_vector, goal_vector, heuristic);
			pathFinding = true;
			solver = new DefaultSolver(sys,
				eps_x_min ? eps_x_min.Get() : DefaultSolver::default_eps_x_min,
				eps_x_max ? eps_x_max.Get() : DefaultSolver::default_eps_x_max,
				buffer,
				random_seed? random_seed.Get() : DefaultSolver::default_random_seed);
		} else {
			solver = new DefaultSolver(sys,
				eps_x_min ? eps_x_min.Get() : DefaultSolver::default_eps_x_min,
				eps_x_max ? eps_x_max.Get() : DefaultSolver::default_eps_x_max,
				!bfs,
				random_seed? random_seed.Get() : DefaultSolver::default_random_seed);
		}

		Solver& s = *solver;
		// Build the default solver
		/*DefaultSolver s(sys,
				eps_x_min ? eps_x_min.Get() : DefaultSolver::default_eps_x_min,
				eps_x_max ? eps_x_max.Get() : DefaultSolver::default_eps_x_max,
				!bfs,
				random_seed? random_seed.Get() : DefaultSolver::default_random_seed);*/

		if (boundary_test_arg) {

			if (boundary_test_arg.Get()=="true")
				s.boundary_test = Solver::ALL_TRUE;
			else if (boundary_test_arg.Get()=="full-rank")
				s.boundary_test = Solver::FULL_RANK;
			else if (boundary_test_arg.Get()=="half-ball")
				s.boundary_test = Solver::HALF_BALL;
			else if (boundary_test_arg.Get()=="false")
				s.boundary_test = Solver::ALL_FALSE;
			else {
				cerr << "\nError: \"" << boundary_test_arg.Get() << "\" is not a valid option (try --help)\n";
				exit(0);
			}

			if (!quiet)
				cout << "  boundary test:\t\t" << boundary_test_arg.Get() << endl;
		}

		if (forced_params) {
			SymbolMap<const ExprSymbol*> symbols;
			for (int i=0; i<sys.args.size(); i++)
				symbols.insert_new(sys.args[i].name, &sys.args[i]);

			string vars=args::get(forced_params);

			vector<const ExprNode*> params;
			int j;
			do {
				j=vars.find("+");
				if (j!=-1) {
					params.push_back(&parse_indexed_symbol(symbols,vars.substr(0,j)));
					vars=vars.substr(j+1,vars.size()-j-1);
 				} else {
 					params.push_back(&parse_indexed_symbol(symbols,vars));
 				}
			} while (j!=-1);

			if (!params.empty()) {
				s.set_params(VarSet(sys.f_ctrs,params,false)); //Array<const ExprNode>(params)));
				for (vector<const ExprNode*>::iterator it=params.begin(); it!=params.end(); it++) {
					cleanup(**it,false);
				}
			}
		}

		// This option limits the search time
		if (timeout) {
			if (!quiet)
				cout << "  timeout:\t\t" << timeout.Get() << "s" << endl;
			s.time_limit=timeout.Get();
		}

		// This option prints each better feasible point when it is found
		if (trace) {
			if (!quiet)
				cout << "  trace:\t\tON" << endl;
			s.trace=trace.Get();
		}

		if (!quiet) {
			cout << "*****************************************************************" << endl << endl;
		}

		if (!quiet)
			cout << "running............" << endl << endl;

		// Get the solutions
		if (input_file)
			s.solve(input_file.Get().c_str());
		else
			s.solve(sys.box);

		if (trace) cout << endl;

		if (!quiet) s.report();

		if (sols) cout << s.get_manifold() << endl;

		if (txt)
			s.get_manifold().write_txt(output_manifold_file.c_str());
		else
			s.get_manifold().write(output_manifold_file.c_str());

		if (!quiet) {
			cout << " results written in " << output_manifold_file << "\n";
			if (overwitten)
				cout << " (old file saved in " << manifold_copy << ")\n";
		}
		//		if (!quiet && !sols) {
//			cout << " (note: use --sols to display solutions)" << endl;
//		}

		if(pathFinding) {
			CellBufferNeighborhood* path_buffer = static_cast<CellBufferNeighborhood*>(buffer);
			if(path_buffer->isPathFound) {
				auto path = path_buffer->pathFound;
				cout << "Begin path:" << endl;
				cout << print_mma_path(path) << endl;
				cout << "End path." << endl;
			} else {
				cout << "No path found." << endl;
			}
		}
		delete solver;
		if(buffer) delete buffer;
	}
	catch(ibex::UnknownFileException& e) {
		cerr << "Error: cannot read file '" << filename.Get() << "'" << endl;
	}
	catch(ibex::SyntaxError& e) {
		cout << e << endl;
	}
	catch(ibex::DimException& e) {
		cout << e << endl;
	}
}
