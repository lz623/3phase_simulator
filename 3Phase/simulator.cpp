#include "DiscreteProblem.hpp"
#include "Intel_Pardiso.hpp"
#include "NewtonSolver.hpp"
#include "R_Precond_Solver.hpp"
#include "Intel_Prec_GMRES.hpp"
#include "Intel_ILU0.hpp"
#include "fastl/containers/pod_vector_unbounded.hpp"
#include <fstream>
#include <dlib/matrix.h>

typedef GENSOL::R_Precond_Solver< GENSOL::Intel_Prec_GMRES< GENSOL::Intel_ILU0 > >  LINEARSOLVER;
//typedef GENSOL::Intel_Prec_GMRES< GENSOL::Intel_ILU0 >  LINEARSOLVER;
typedef GENSOL::NewtonSolver< DiscreteProblem, LINEARSOLVER > STDN;

class simulator
{
public:

	void
	read_from_file(const char * filename, std::vector<double> &v)
	{
		double tmp;
		std::ifstream strm(filename);
		if (strm.good())
		{
			for (std::size_t i = 0; i < v.size(); ++i)
			if (strm >> tmp) v[i] = tmp;
		}
	}

	void
		dump_solution(const char * filename, const DiscreteProblem::StateVector & v)
	{
			std::ofstream strm(filename);
			for (std::size_t i = 0; i < v.size(); ++i)
				strm << v[i].Po.value() << "\t"
				<< v[i].Sw.value() << "\t"
				<< v[i].Sg.value() << "\t"
				<< 1.0 - v[i].Sw.value() - v[i].Sg.value() << "\t"
				<< v[i].Rso.value() << std::endl;
			strm.close();
		}

	void
		dump_field(const char * filename,
		const std::vector<double> &phi,
		const std::vector<double> &kx,
		const std::vector<double> &ky,
		const std::vector<double> &kz)
	{
			std::ofstream strm(filename);
			for (std::size_t i = 0; i < phi.size(); ++i)
				strm << phi[i] << "\t"
				<< kx[i] << "\t"
				<< ky[i] << "\t"
				<< kz[i] << std::endl;
			strm.close();
		}

	void read_reservoir()
	{
		std::ifstream strm("./input/reservoir.dat");
		strm >> NX;
		strm >> NY;
		strm >> NZ;
		strm >> LX;
		strm >> LY;
		strm >> LZ;
	}
	void read_schedule(std::vector<double> &times)
	{
		double tmp;
		std::ifstream strm("./input/schedule.dat");
		while (strm.good())
		{
			strm >> tmp;
			times.push_back(tmp);

		}
	}
	void m_to_vector(dlib::matrix<double> &si,
		std::vector<double> &phi,
		std::vector<double> &kx,
		std::vector<double> &ky,
		std::vector<double> &kz)
	{
		size_t n = si.size();
		size_t nc = n / 4;

		for (unsigned i = 0; i < nc; i++)
		{
			kx[i] = exp(si(i, 0));
			ky[i] = exp(si(i + nc, 0));
			kz[i] = exp(si(i + 2 * nc, 0));
			phi[i] = si(i + 3 * nc, 0);
		}
	}
	simulator()
	{
		read_reservoir();
		read_schedule(times);
		vPORO.resize(NX*NY*NZ);
		vKX.resize(NX*NY*NZ);
		vKY.resize(NX*NY*NZ);
		vKZ.resize(NX*NY*NZ);
	}

	void run(dlib::matrix<double> &m, dlib::matrix<double> &v, 
		dlib::matrix<double> &d,double &NPV)
	{
		m_to_vector(m, vPORO, vKX, vKY, vKZ);
		std::ofstream out("./output/timestep.out");
		dump_field("./output/field.out", vPORO, vKX, vKY, vKZ);
		const double OPT_NLNITER = (3 < MAX_NLNITER ? MAX_NLNITER : 3);
		DiscreteProblem model(NX, NY, NZ, LX, LY, LZ, vKX, vKY, vKZ, vPORO);
		DiscreteProblem::StateVector uOld, uNew;
		model.initialize_state(uOld);
		LINEARSOLVER lnsolver(model.max_num_eqns(),
			model.max_num_nnz());
		STDN newton(model, MAX_NLNITER, 1);
		double DT = DT_INIT;
		double time = 0.0;
		for (std::size_t i = 0; i < times.size(); i++)
		{
			while (time < times[i])
			{
				uNew = uOld;
				STDN::report_t stdsmry = newton.solve_timestep(uNew, uOld, DT, model, lnsolver);
				if (stdsmry.is_converged)
				{
					uOld = uNew;
					time += DT;
					if (stdsmry.niter < OPT_NLNITER) DT *= DT_GROW;
					DT = std::min(DT, (times[i] - time));
					DT == 0 ? DT = DT_INIT : DT = DT;
					std::cout << "CONVERGED t = " << time << " days" << std::endl;
					out << DT << std::endl;
				}
				else
				{
					DT *= DT_CUT;
					std::cout << "FAILED " << std::endl;
				}
			}
		}
	}

		private:
			std::vector<double> times;
			const std::size_t MAX_NLNITER = 12;
			const double      DT_INIT = 5.0;
			const double      DT_CUT = 0.5;
			const double      DT_GROW = 2.0;
			const double      T_FINAL = 365.0;
			std::size_t NX, NY, NZ;
			double      LX, LY, LZ;
			std::vector<double> vPORO;
			std::vector<double> vKX;
			std::vector<double> vKY;
			std::vector<double> vKZ;
};