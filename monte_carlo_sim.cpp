#include <iostream>
#include <random>
#include <vector>
#include <fstream>
#include <string>
#include <omp.h>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <cstring>


std::vector<size_t> rankSort(const std::vector<double>& v_temp) {
	std::vector<std::pair<double, size_t> > v_sort(v_temp.size());

	for (size_t i = 0U; i < v_sort.size(); ++i) {
		v_sort[i] = std::make_pair(v_temp[i], i);
	}

	std::sort(v_sort.begin(), v_sort.end());

	std::pair<double, size_t> rank;
	std::vector<size_t> result(v_temp.size());

	for (size_t i = 0U; i < v_sort.size(); ++i) {
		if (v_sort[i].first != rank.first) {
			rank = std::make_pair(v_sort[i].first, i);
		}
		result[v_sort[i].second] = rank.second + 1;
	}
	return result;
}



int main(int argc, char **argv)
{
	int thread_n, sim_n;
	std::string eta_file, N_file, LRT_file, omega_file, date, file_n;

	std::istringstream argv_thread(argv[1]);
	std::istringstream argv_sim(argv[2]);
	std::istringstream argv_eta(argv[3]);
	std::istringstream argv_N(argv[4]);
	std::istringstream argv_LRT(argv[5]);
	std::istringstream argv_omega(argv[6]);
	std::istringstream argv_date(argv[7]);
	std::istringstream argv_file(argv[8]);

	if (!(argv_thread >> thread_n))
	{
		std::cerr << "Invalid number for multicore: " << argv[1] << std::endl;
		exit(1);
	}
	if (!(argv_sim >> sim_n))
	{
		std::cerr << "Invalid number for simulation time: " << argv[2] << std::endl;
		exit(1);
	}
	if (!(argv_eta >> eta_file))
	{
		std::cerr << "Invalid file name for eta: " << argv[3] << std::endl;
		exit(1);
	}
	if (!(argv_N >> N_file))
	{
		std::cerr << "Invalid file name for N: " << argv[4] << std::endl;
		exit(1);
	}
	if (!(argv_LRT >> LRT_file))
	{
		std::cerr << "Invalid file name for LRT: " << argv[5] << std::endl;
		exit(1);
	}
	if (!(argv_omega >> omega_file))
	{
		std::cerr << "Invalid file name for omega: " << argv[6] << std::endl;
		exit(1);
	}
	if (!(argv_date >> date))
	{
		std::cerr << "Invalid file name for date: " << argv[7] << std::endl;
		exit(1);
	}
	if (!(argv_file >> file_n))
	{
		std::cerr << "Invalid file number: " << argv[8] << std::endl;
		exit(1);
	}

	std::vector<std::vector<double>> eta;
	std::vector<std::vector<int>> N;
	std::vector<std::vector<double>> omega;
	std::vector<std::vector<double>> LRT;
	

	std::ifstream IN_ETA(eta_file);
	if (IN_ETA.is_open())
	{
		std::string str;
		while (std::getline(IN_ETA, str))
		{
			std::vector<double> vec_tmp;
			const char *loc = str.c_str();
			vec_tmp.push_back(std::atof(loc));
			loc = std::strstr(loc, "\t");
			while (loc != NULL)
			{
				vec_tmp.push_back(atof(loc + 1));
				loc = std::strstr(loc + 1, "\t");
			}
			eta.push_back(vec_tmp);
		}
		IN_ETA.close();
	}
	else
	{
		fprintf(stderr, "There was an error opening the eta file.\n");
		exit(1);
	}

	std::ifstream IN_N(N_file);
	if (IN_N.is_open())
	{
		std::string str;
		while (std::getline(IN_N, str))
		{
			std::vector<int> vec_tmp;
			const char *loc = str.c_str();
			vec_tmp.push_back(std::atoi(loc));
			loc = std::strstr(loc, "\t");
			while (loc != NULL)
			{
				vec_tmp.push_back(atoi(loc + 1));
				loc = std::strstr(loc + 1, "\t");
			}
			N.push_back(vec_tmp);
		}
		IN_N.close();
	}
	else
	{
		fprintf(stderr, "There was an error opening the N file.\n");
		exit(1);
	}

	std::ifstream IN_OMEGA(omega_file);
	if (IN_OMEGA.is_open())
	{
		std::string str;
		while (std::getline(IN_OMEGA, str))
		{
			std::vector<double> vec_tmp;
			const char *loc = str.c_str();
			vec_tmp.push_back(std::atof(loc));
			loc = std::strstr(loc, "\t");
			while (loc != NULL)
			{
				vec_tmp.push_back(atof(loc + 1));
				loc = std::strstr(loc + 1, "\t");
			}
			omega.push_back(vec_tmp);
		}
		IN_OMEGA.close();
	}
	else
	{
		fprintf(stderr, "There was an error opening the omega file.\n");
		exit(1);
	}

	std::ifstream IN_LRT(LRT_file);
	if (IN_LRT.is_open())
	{
		std::string str;
		while (std::getline(IN_LRT, str))
		{
			std::vector<double> vec_tmp;
			const char *loc = str.c_str();
			vec_tmp.push_back(std::atof(loc));
			loc = std::strstr(loc, "\t");
			while (loc != NULL)
			{
				vec_tmp.push_back(atof(loc + 1));
				loc = std::strstr(loc + 1, "\t");
			}
			LRT.push_back(vec_tmp);
		}
		IN_LRT.close();
	}
	else
	{
		fprintf(stderr, "There was an error opening the LRT file.\n");
		exit(1);
	}

	if (!(N.size() == omega.size() && N.size() == LRT.size() && N[0].size() == eta.size() && eta.size() == omega[0].size()))
	{
		std::cerr << "dimensions were not compatible: " << N.size() << " " << omega.size() << " " << LRT.size() << " " << N[0].size() << " " << eta.size() << " " << omega[0].size() << std::endl;
	}
	int gene_n = N.size();
	int sample_n = eta[0].size();
	int type_n = eta.size();
	std::vector<double> abs_LRT(gene_n);

	for (size_t k = 0; k < gene_n; k++)
	{
		abs_LRT[k] = std::abs(LRT[k][0]);
	}

	std::vector<double> p2sides_noN(gene_n);

	std::vector<double> p2sides_noN_adj(gene_n);

	omp_set_num_threads(thread_n);
	
	#pragma omp declare reduction(vec_plus : std::vector<double> : \
											 std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                    initializer(omp_priv = omp_orig)

	#pragma omp parallel for reduction(vec_plus:p2sides_noN)
	for (int k = 0; k < gene_n; k++)
	{
		std::default_random_engine generator(k);
		std::poisson_distribution<> rpois;
		for (int t = 0; t < sim_n; t++)
		{
			
			double numerator = 0;
			
			double denominator_noN = 0;
			for (int j = 0; j < type_n; j++)
			{
				double numerator_help = 0;
				
				double denominator_noN_help = 0;
				for (int i = 0; i < sample_n; i++)
				{
					rpois = std::poisson_distribution<>(N[k][j] * eta[j][i]);
					numerator_help += rpois(generator) / eta[j][i];
					
					denominator_noN_help += 1.0 / eta[j][i];
				}
				numerator_help -= sample_n * N[k][j];
				numerator += omega[k][j] * numerator_help;
				
				denominator_noN += omega[k][j] * omega[k][j] * denominator_noN_help;
			}
			
			
			double monte_carlo_T_noN = numerator / std::sqrt(denominator_noN);
			
			double abs_monte_carlo_T_noN = std::abs(monte_carlo_T_noN);
			
			for (int pk = 0; pk < gene_n; pk++)
			{
				
				if (abs_monte_carlo_T_noN >= abs_LRT[pk])
				{
					p2sides_noN[pk]++;
				}
			}
		}
		
		
	}


	std::transform(p2sides_noN.begin(), p2sides_noN.end(), p2sides_noN.begin(), [&sim_n, &gene_n](double value) {return value / (sim_n*gene_n); });
	std::vector<size_t> p2sides_noN_rank = rankSort(p2sides_noN);
	std::transform(p2sides_noN.cbegin(), p2sides_noN.cend(), p2sides_noN_rank.cbegin(), p2sides_noN_adj.begin(), [&gene_n](const double &p, const size_t &rank) {return p * gene_n / rank; });
	
	std::ofstream ofile_p(date + "_" + file_n + "_p.tmp");

	if (ofile_p.is_open())
	{
		for (size_t k = 0; k < gene_n; k++)
		{
			ofile_p << p2sides_noN[k] << "\t" << p2sides_noN_adj[k] << std::endl;
		}
		
		ofile_p.close();
	}
	else
	{
		std::cout << "opening error for the output file for p" << std::endl;
	}

    return 0;
}

