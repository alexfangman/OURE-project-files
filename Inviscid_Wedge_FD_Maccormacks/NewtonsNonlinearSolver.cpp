#include "stdafx.h"
#include "NewtonsNonlinearSolver.h"

void Calc_Primitive_Vars(const vector< vector<double> > &F, vector< vector<double> > &prim_vars, const double xi)
{
	double A, B, C, inner, rho, u, v, p, T;

	for (int i = 0; i <= (F.size() - 1); i++)
	{
		A = (pow(F[i][2] , 2) / (2 * F[i][0])) - F[i][3];
		B = (gamma / (gamma - 1)) * F[i][0] * F[i][1];
		C = -(( gamma + 1 ) / ( 2 * ( gamma - 1 ) )) * pow(F[i][0] , 3);
		inner = pow(B , 2) - (4 * A * C);

		//density
		rho = ( -B + pow(inner , 0.5) ) / (2 * A);
		//u
		u = F[i][0] / rho;

		if (xi < wedge_location)
		{
			v = 0;
		}
		else
		{
			//v
			v = F[i][2] / F[i][0];
		}

		
		
		//p
		p = F[i][1] - ( F[i][0] * u );
		//T
		T = p / (rho*R);

		//Place in prim_vars vector
		prim_vars[i][0] = u;
		prim_vars[i][1] = v;
		prim_vars[i][2] = p;
		prim_vars[i][3] = T;
		prim_vars[i][4] = rho;
	}

	return;
}