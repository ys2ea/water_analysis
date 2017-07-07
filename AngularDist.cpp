#include "AngularDist.h"


int main(int argc, char * argv[]) {
	int N;              //bug here if N is not const, don't know why yet!!!  ---fixed
	double step;
	double L;
	double blength;
	double mlength, OOlength, HOOangle;
	int startiter, maxiter, Np;
	
	double Mx=0., My=0., Mz=0., Msq=0.;
	
	std::ifstream input;
	input.open(argv[1]);
	if(input==NULL) { 
		std::cout << "FEED ME THE TRAJECTORY FILE!!\n";
		exit(1); }
		
	std::ifstream paraminput;
	paraminput.open(argv[2]);
	if(paraminput==NULL) { 
		std::cout << "FEED ME THE PARAM FILE!!\n";
		exit(1); }
	
	paraminput >> N;
	paraminput >> step;
	paraminput >> L;
	paraminput >> blength;
	paraminput >> mlength;
	paraminput >> OOlength;
	paraminput >> HOOangle;
	paraminput >> startiter;
	paraminput >> maxiter;
	
	paraminput.close();
	
	MD_system a(N, step, L);
	int iter = 0;
	std::string line;
	double Hstat[8];
	for(int j = 0; j < 8; j++) Hstat[j] = 0.;
	
	FILE* doutput;
	doutput = fopen("dipole.dat","w");
	
	int count = 0;
	while(iter < maxiter) {

		getline(input, line);
		getline(input, line);
		//be sure to change latter
		for(int i = 0; i< N; i ++) {
			getline(input, line);
			if(iter == 0) a.init_read(line, i);
			else  a.pos_read(line,i);

		}
		
		for(int i = 0; i< N*4/3; i ++) {
			getline(input, line);
			if(iter == 0) a.init_read(line, i);
			else  a.pos_read(line,i);

		}

		if(iter == 0) a.decide_molecule(blength);
		
		if(iter >= startiter && iter%20 == 1) {
			//a.fill_Ahist(1.3, mlength);
			if(count%1000 == 1) {
			std::cout << "iter :  " << iter << "\n";
			}
			
			double* Hc = a.H_bond_stat(OOlength, HOOangle, "ad");
			
			for(int k = 0; k < 8; k ++)
				Hstat[k] = Hstat[k]*count/(count+1) + Hc[k]/(count+1);
			

			double *dp = a.mol_dipole();
			fprintf(doutput, "%lf\t%lf\t%lf\n", dp[0], dp[1], dp[2]);
			Mx = Mx*count/(count+1) + dp[0]/(count+1);
			My = My*count/(count+1) + dp[1]/(count+1);
			Mz = Mz*count/(count+1) + dp[2]/(count+1);
			Msq = Msq*count/(count+1) + (dp[0]*dp[0]+dp[1]*dp[1]+dp[2]*dp[2])/(count+1);
			count ++;
		}  

		iter ++;
	}
	input.close();
	//a.Show_hist("angular_dist.txt");
	
	std::cout << " Average Dipole: (" << Mx << "," << My << "," << Mz << ")"  << "Variance of dipole: " << Msq - Mx*Mx - My*My - Mz*Mz << "\n";	
	for(int i = 0; i < 8; i++)
		std::cout << Hstat[i] << ",";
	std::cout << "\n";
	std::cout << Hstat[0]+Hstat[1]+Hstat[2]+Hstat[3] << "\n";
	std::cout << Hstat[4]+Hstat[5]+Hstat[6]+Hstat[7] << "\n";
	
	fclose(doutput);
	return 1;
}
