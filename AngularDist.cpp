#include "AngularDist.h"

//a1 should be t=0. a1 and a2 should be ordered
int HB_overlap(std::vector<std::pair<int,int> > a1, std::vector<std::pair<int,int> > a2) {
	size_t m= 0, n = 0; 
	int count = 0;

	while(m<a1.size() && n < a2.size()) {
		if(a1[m].first < a2[n].first || (a1[m].first==a2[n].first && a1[m].second < a2[n].second)) m++;
		else if(a1[m].first==a2[n].first && a1[m].second==a2[n].second) {m++; n++; count ++; }
		//std::cout << "(" << a1[m].first << "," a1[m].second << ") ; (" << a2[n].first << "," << a2[n].second << ")\n";}
		else {n ++;}
	}
	
	return count;
}

int HB_overlap(std::vector<std::pair<int,int> > a1, std::vector<std::pair<int,int> > a2, int* (&HB_count), int* (&hist), int N) {
	size_t m= 0, n = 0; 
	int count = 0;

	while(m<a1.size() && n < a2.size()) {
		int x = a1[m].first + a1[m].second*N/3;
		int y = a2[n].first + a2[n].second*N/3;
		if(a1[m].first < a2[n].first || (a1[m].first==a2[n].first && a1[m].second < a2[n].second)) 
			{hist[HB_count[x]] ++; HB_count[x] = 0; m ++;}
		else if(a1[m].first==a2[n].first && a1[m].second==a2[n].second) {HB_count[x] ++; m++; n++; count ++; }
		//std::cout << "(" << a1[m].first << "," a1[m].second << ") ; (" << a2[n].first << "," << a2[n].second << ")\n";}
		else {n ++; HB_count[y] = 1;}
	}
	
	return count;
}


int main(int argc, char * argv[]) {
	int N;              //bug here if N is not const, don't know why yet!!!  ---fixed
	double step;
	double L;
	double blength;
	double mlength, OOlength, HOOangle;
	int startiter, maxiter;
	std::vector<std::vector<std::pair<int,int> > > HB_map;;
	double Mx=0., My=0., Mz=0., Msq=0.;
	double timestep = 0.5;   //MD step in fs;
	
	bool do_angular = 1;
	bool do_HBstat = 0;
	bool do_HBlife = 0;
	bool do_dipole = 0;
	
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
		
	/*	for(int i = 0; i< N*4/3; i ++) {
			getline(input, line);
			if(iter == 0) a.init_read(line, i);
			else  a.pos_read(line,i);

		}  */                 //have to include when doing dipole. Pay attention to this, improve later.

		if(iter == 0) a.decide_molecule(blength);
		
		if(iter >= startiter) {
			
			if(count%1000 == 1) {
			std::cout << "iter :  " << iter << "\n";
			}
			
			if(do_angular) a.fill_Ahist(2, mlength);
			
			if(do_HBstat) {
				double* Hc = a.H_bond_stat(OOlength, HOOangle, "ad");
				for(int k = 0; k < 8; k ++)
					Hstat[k] = Hstat[k]*count/(count+1) + Hc[k]/(count+1);
			}
			
			if(do_HBlife) HB_map.push_back(a.H_bond_map(OOlength, HOOangle, "sf"));
			
			
			if(do_dipole) {

				double *dp = a.mol_dipole();
				//fprintf(doutput, "%lf\t%lf\t%lf\n", dp[0], dp[1], dp[2]);
				Mx = Mx*count/(count+1) + dp[0]/(count+1);
				My = My*count/(count+1) + dp[1]/(count+1);
				Mz = Mz*count/(count+1) + dp[2]/(count+1);
				Msq = Msq*count/(count+1) + (dp[0]*dp[0]+dp[1]*dp[1]+dp[2]*dp[2])/(count+1);
			}
			
			count ++;
		}  //if iter  

		iter ++;
	}
	input.close();
	
	if(do_angular)	a.Show_hist("angular_dist.txt");
	
	
	//post read-in calculation, mostly for HB dynamics;  
	
	
	//This doesn't take account of HB reforming
/*	int tau_max = 10000;						//range of the correlation function, in unit of steps.
	std::vector<double> corr;	//correlation function;
	for(int jj = 0;  jj < tau_max; ++ jj ) {
		int c0=0, ct=0;
		for(size_t kk = 0; kk < HB_map.size() - jj; ++ kk) {
		

			c0 += HB_map[kk].size();
			ct += HB_overlap(HB_map[kk],HB_map[kk+jj]);
			
		if(jj == 3000 && kk == 1000) {
			for(size_t ii = 0; ii < HB_map[kk].size(); ++ ii) std::cout << "(" << HB_map[kk][ii].first << "," << HB_map[kk][ii].second << ")"; 
			std::cout << "\n";
			for(size_t ii = 0; ii < HB_map[kk+jj].size(); ++ ii) std::cout << "(" << HB_map[kk+jj][ii].first << "," << HB_map[kk+jj][ii].second << ")";
			
			std::cout << "------(" << HB_map[kk].size() << "," << HB_overlap(HB_map[kk],HB_map[kk+jj]) << ")\n";
		}
			
			//std::cout << "(" << c0 << "," << ct << ")\n";
		}
		corr.push_back(double(ct)/double(c0));
		
	}  */
	
	
	if(do_HBlife) {
		int* HB_count;
		int* hist;
		int Nstep = 20000;
		HB_count = (int*)malloc(sizeof(int)*N*N/9);
		hist = (int*) malloc(sizeof(int)*Nstep);

		for(int ii = 0; ii < N*N/9; ++ii)
			HB_count[ii] = 0;
			
		for(int ii = 0; ii < Nstep; ++ii)
			hist[ii] = 0;
		
		for(size_t kk = 0; kk < HB_map[1].size(); ++kk) {
			int x = HB_map[1][kk].first + N/3*HB_map[1][kk].second;
			HB_count[x] = 1;
		}	
		for(size_t jj = 1; jj < HB_map.size()-1; ++jj) {
			HB_overlap(HB_map[jj], HB_map[jj+1], HB_count, hist, N);
		}
		FILE* corrout;
		corrout = fopen("corr.dat", "w");

		int* corr;
		corr = (int*) malloc(sizeof(int)*Nstep);
		corr[0] = 0;
		for(int ii = 0; ii < Nstep; ++ii)
			corr[0] += hist[ii];
		
		for(int ii = 1; ii < Nstep; ++ii)
			corr[ii] = corr[ii-1] - hist[ii-1];
		for(int i = 0; i < Nstep; i ++)
			fprintf(corrout, "%lf\t%lf\n", i*timestep/1000. , double(corr[i])/double(corr[0]));
	}	
	
	if(do_dipole) {
		std::cout << " Average Dipole: (" << Mx << "," << My << "," << Mz << ")"  << "Variance of dipole: " << Msq - Mx*Mx - My*My - Mz*Mz << "\n";	
	}
	
	if(do_HBstat) {
		for(int i = 0; i < 8; i++)
			std::cout << Hstat[i] << ",";
		std::cout << "\n";
		std::cout << Hstat[0]+Hstat[1]+Hstat[2]+Hstat[3] << "\n";
		std::cout << Hstat[4]+Hstat[5]+Hstat[6]+Hstat[7] << "\n";
	}
	
	return 1;
}


