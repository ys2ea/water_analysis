/* A code to read cp2k output and 
** calculate angular distribution function
** and maybe other quantities...
** This is code is written for periodic bc
** 
** Created by Yifei Shi 2017.2.1
** Modified 2017. 4. 25, adding H-bond statistics
**
** This is not only a code that do angules, but also analyze many other properties.
** Added dielectric constant 2017/6/29
*/ 

#ifndef ANGULARDIST_H
#define ANGULARDIST_H

const double PI = 3.1415926536;

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

class atom {
private:
	char type_;
	int idx_;
	double x_, y_, z_;
	//int charge_;
	
public:
	atom(char, int, double, double, double);
	std::vector<int> bond_list;
	void set_coord(const double, const double, const double);
	double x();
	double y();
	double z();
	char type();
	//set_atom_bond();
};

class molecule {
private:
	int Natoms_;
	std::vector<atom*> a_mol_;
	std::vector<molecule*> nbmolecule_;
public:
	molecule(atom);
	void add_atom(atom *);
	std::vector<atom*> show_atom();
	int num_atom();
};


class MD_system {
private: 
	int nmolecule_;
	int natom_;
	std::vector<atom> alist_;   //List of atoms in the system.
	std::vector<atom> wcenterlist_;    //List of wannier centers, they are treated as an atom, but not listed in the alist.
	std::vector<int> mlist_;    //List of molecules in the system. Right now it's just the index of the O atom, should make a molecule class....
	std::vector<molecule> mmlist_;
	std::vector<int> Ahist_;
	double L_;	//assume cubic symmetry
	double step_;
public:
	MD_system(const int, const double, const double);
	void init_read(const std::string &, const int);	//initial read, which create all the atoms
	void pos_read(const std::string &, const int &);		//update read, only change the positions of the atoms
										//Note that this code doesn't consider velocities at all
	void decide_molecule(double);		//find covalent bonds using geometry criterion
	void fill_Ahist(double, double); 		//fill hist using geometry criterion
	void Show_hist(const char*);
	double* H_bond_stat(const double, const double, const char*);  //Decide how many H-bond a molecule forms.
													//param: distance(max), angle(max), outputfile name 
	double dipole();	//wrong
	double* mol_dipole();    //this one attatches wannier centers to molecules										
};


molecule::molecule( atom a) {
	Natoms_ = 1;
	a_mol_.push_back(&a);	
}

int molecule::num_atom() {
	return Natoms_;
}

std::vector<atom*> molecule::show_atom() {
	return a_mol_;
}

atom::atom(char type, int idx, double x, double y, double z) {
	type_ = type;
	idx_ = idx;
	x_ = x;
	y_ = y;
	z_ = z;
}


void atom::set_coord(const double x, const double y, const double z) {
	x_ = x; y_ = y; z_ = z;
}

double atom::x() {
	return x_;
}

double atom::y() {
	return y_;
}

double atom::z() {
	return z_;
}

char atom::type() {
	return type_;
}

MD_system::MD_system(const int na, const double step, const double L) {
	natom_ = na;
	nmolecule_ = na/3;
	step_ = step;
	Ahist_.resize(int(90./step));
	std::fill(Ahist_.begin(), Ahist_.end(), 0);
	L_ = L;
}

void MD_system::init_read(const std::string &line, const int idx) {
	char type[10];
	double x, y, z;
	sscanf(line.c_str(), "%s %lf %lf %lf", type, &x, &y, &z);
	if(type[0] != 'X') 
		alist_.push_back(atom(type[0], idx, x, y, z));
	else
		wcenterlist_.push_back(atom(type[0], idx, x, y, z));
}

void MD_system::pos_read(const std::string &line, const int & idx) {
	char type[10];
	double x, y, z;
	sscanf(line.c_str(), "%s\t%lf\t%lf\t%lf", type, &x, &y, &z);
	if(type[0] != 'X')
		alist_[idx].set_coord(x, y, z);
	else
		wcenterlist_[idx].set_coord(x, y, z);
}

void MD_system::decide_molecule(double bl) {
	double dist;
	for(int i = 0; i < natom_; ++ i) {
		//store the index of O atoms
		if(alist_[i].type()=='O') {
			mlist_.push_back(i);
			//mmlist_.push_back(molecule(alist_[i]));
		}
			
		for(int j = i + 1; j < natom_; ++ j) {

				double dx = alist_[i].x()-alist_[j].x();
				double dy = alist_[i].y()-alist_[j].y();
				double dz = alist_[i].z()-alist_[j].z();
			
				if(dx > L_/2.) dx = L_-dx;
				if(dx < -L_/2.) dx = L_+dx;
				if(dy > L_/2.) dy = L_-dy;
				if(dy < -L_/2.) dy = L_+dy;
				if(dz > L_/2.) dz = L_-dz;
				if(dz < -L_/2.) dz = L_+dz;
									 
				dist = dx*dx + dy*dy + dz*dz;
				  
				if(dist < bl) {
					alist_[i].bond_list.push_back(j);
					alist_[j].bond_list.push_back(i);
				}

		}
	}
	

	//A check specifically for water molecule
	for(int i = 0; i < natom_; ++ i) {
		if(alist_[i].type()=='O' && alist_[i].bond_list.size() != 2)
			std::cout << "Oxegen (" << alist_[i].x() << "," << alist_[i].y() << "," << alist_[i].z() << ") have wrong bonds!\n";
		if(alist_[i].type()=='H' && alist_[i].bond_list.size() != 1)
			std::cout << "Hydrogen (" << alist_[i].x() << "," << alist_[i].y() << "," << alist_[i].z() << ") have wrong bonds!\n";
		//std::cout << "(" << alist_[i].type() << "," << alist_[i].bond_list.size() << ")\n";
		
	}
	
	
	//attatch wannier centers to the molecules
}

//poorly written, just want it to work.
//Fill the hist of HB angle distribution
void MD_system::fill_Ahist(double bsq1, double bsq2) {
	double dist, dx, dy, dz;
	for(int i = 0; i < natom_; ++ i) {
		for(int j = i+1; j < natom_; ++ j) {
			dx = alist_[i].x()-alist_[j].x();
			dy = alist_[i].y()-alist_[j].y();
			dz = alist_[i].z()-alist_[j].z();
			
			if(dx > L_/2.) dx -= L_;
			if(dx < -L_/2.) dx += L_;
			if(dy > L_/2.) dy -= L_;
			if(dy < -L_/2.) dy += L_;
			if(dz > L_/2.) dz -= L_;
			if(dz < -L_/2.) dz += L_;
			
			dist = dx*dx + dy*dy + dz*dz;
			  
			if(dist > bsq1 && dist < bsq2) {
				if(alist_[i].type() == 'O' && alist_[j].type() == 'H') {
				
					//the usual definition of H bond
					double vx1 = alist_[alist_[j].bond_list[0]].x()-alist_[i].x();
					double vy1 = alist_[alist_[j].bond_list[0]].y()-alist_[i].y();
					double vz1 = alist_[alist_[j].bond_list[0]].z()-alist_[i].z();
					
					double vx2 = alist_[alist_[j].bond_list[0]].x()-alist_[j].x();
					double vy2 = alist_[alist_[j].bond_list[0]].y()-alist_[j].y();
					double vz2 = alist_[alist_[j].bond_list[0]].z()-alist_[j].z(); 
					
					//Check if the O and O atoms are far away, 
					//also change the range of angle to 0-180 if this is used !!!!!!!!!!!!!!!!!!!!!!
				/*	double vx1 = alist_[j].x()-alist_[i].x();
					double vy1 = alist_[j].y()-alist_[i].y();
					double vz1 = alist_[j].z()-alist_[i].z();
					
					double vx2 = -alist_[alist_[j].bond_list[0]].x()+alist_[j].x();
					double vy2 = -alist_[alist_[j].bond_list[0]].y()+alist_[j].y();
					double vz2 = -alist_[alist_[j].bond_list[0]].z()+alist_[j].z(); */
					
					if(vx1 > L_/2.) vx1 = vx1-L_;
					if(vx1 < -L_/2.) vx1 = vx1+L_;
					if(vy1 > L_/2.) vy1 = vy1-L_;
					if(vy1 < -L_/2.) vy1 = L_+vy1;
					if(vz1 > L_/2.) vz1 = vz1-L_;
					if(vz1 < -L_/2.) vz1 = L_+vz1;
					if(vx2 > L_/2.) vx2 = vx2-L_;
					if(vx2 < -L_/2.) vx2 = L_+vx2;
					if(vy2 > L_/2.) vy2 = vy2-L_;
					if(vy2 < -L_/2.) vy2 = L_+vy2;
					if(vz2 > L_/2.) vz2 = vz2-L_;
					if(vz2 < -L_/2.) vz2 = L_+vz2; 
					
					double angle = 180.*acos((vx1*vx2+vy1*vy2+vz1*vz2)/sqrt(vx1*vx1+vy1*vy1+vz1*vz1)/sqrt(vx2*vx2+vy2*vy2+vz2*vz2))/PI;
					
					if(angle>0. && angle<90.)
					Ahist_[int(angle/step_)] ++;
				}
				
				else if(alist_[i].type() == 'H' && alist_[j].type() == 'O') {
					//the usual definition of H bond
					double vx1 = alist_[alist_[i].bond_list[0]].x()-alist_[i].x();
					double vy1 = alist_[alist_[i].bond_list[0]].y()-alist_[i].y();
					double vz1 = alist_[alist_[i].bond_list[0]].z()-alist_[i].z();
					
					double vx2 = alist_[alist_[i].bond_list[0]].x()-alist_[j].x();
					double vy2 = alist_[alist_[i].bond_list[0]].y()-alist_[j].y();
					double vz2 = alist_[alist_[i].bond_list[0]].z()-alist_[j].z(); 
					
					//Check if the O and O atoms are far away
				/*	double vx1 = alist_[alist_[i].bond_list[0]].x()-alist_[i].x();
					double vy1 = alist_[alist_[i].bond_list[0]].y()-alist_[i].y();
					double vz1 = alist_[alist_[i].bond_list[0]].z()-alist_[i].z();
					
					double vx2 = alist_[j].x()-alist_[i].x();
					double vy2 = alist_[j].y()-alist_[i].y();
					double vz2 = alist_[j].z()-alist_[i].z(); */
					
					if(vx1 > L_/2.) vx1 = vx1-L_;
					if(vx1 < -L_/2.) vx1 = vx1+L_;
					if(vy1 > L_/2.) vy1 = vy1-L_;
					if(vy1 < -L_/2.) vy1 = L_+vy1;
					if(vz1 > L_/2.) vz1 = vz1-L_;
					if(vz1 < -L_/2.) vz1 = L_+vz1;
					if(vx2 > L_/2.) vx2 = vx2-L_;
					if(vx2 < -L_/2.) vx2 = L_+vx2;
					if(vy2 > L_/2.) vy2 = vy2-L_;
					if(vy2 < -L_/2.) vy2 = L_+vy2;
					if(vz2 > L_/2.) vz2 = vz2-L_;
					if(vz2 < -L_/2.) vz2 = L_+vz2; 
					
					
					double angle = 180.*acos((vx1*vx2+vy1*vy2+vz1*vz2)/sqrt(vx1*vx1+vy1*vy1+vz1*vz1)/sqrt(vx2*vx2+vy2*vy2+vz2*vz2))/PI;
					
					if(angle>0. && angle<90.)
					Ahist_[int(angle/step_)] ++;				
				}
			}
		}
	}
}

//this doesn't work
double MD_system::dipole() {
	double dx = 0., dy = 0., dz = 0.;
	double px, py, pz;
	int ncharge;
	int c = 0;
	for(size_t i = 0; i < alist_.size(); ++ i ) {
		if(alist_[i].type()=='H') {
			ncharge = 1;
		}
		else if(alist_[i].type()=='O') {
			ncharge = 6;
		}
		else {
			std::cout << "Not water molecule?\n";
		}
		px = alist_[i].x() > L_/2 ? alist_[i].x()-L_ : alist_[i].x();
		py = alist_[i].y() > L_/2 ? alist_[i].y()-L_ : alist_[i].y();
		pz = alist_[i].z() > L_/2 ? alist_[i].z()-L_ : alist_[i].z();
		
		dx += px*ncharge;
		dy += py*ncharge;
		dz += pz*ncharge;
		
		c += ncharge;
	}
	
	for(size_t i = 0; i < wcenterlist_.size(); ++ i ) {
		px = wcenterlist_[i].x() > L_/2 ? wcenterlist_[i].x()-L_ : wcenterlist_[i].x();
		py = wcenterlist_[i].y() > L_/2 ? wcenterlist_[i].y()-L_ : wcenterlist_[i].y();
		pz = wcenterlist_[i].z() > L_/2 ? wcenterlist_[i].z()-L_ : wcenterlist_[i].z();
		
		dx += -2.*px;
		dy += -2.*py;
		dz += -2.*pz;
		
		c -= 2;
	}
	
	if(c!=0) std::cout << "Total charge nonezero!!!!!!!!!!!!!!!!!\n";
	
	return sqrt(dx*dx+dy*dy+dz*dz);
	
}

//dipole for each molecule. Assign wannier centers to molecules
double* MD_system::mol_dipole() {

	std::vector<std::vector<int> > mol_wc;
	
	double dx, dy, dz, distsq;
	double Mx=0., My=0., Mz=0.;
	double daverage = 0.;
	
	for(size_t i = 0; i < mlist_.size(); ++i ) {

		//std::vector<int> wcs;
		std::vector<std::vector<double> > ccarray;  //array that contains <charge, x, y, z> of atoms in a molecule. the coord are relative to O.
//		std::cout << "M " << i << " Coord: " << "(" << alist_[mlist_[i]].x() << "," << alist_[mlist_[i]].y() << "," << alist_[mlist_[i]].z() << ")\n";
		if(alist_[mlist_[i]].bond_list.size() != 2) std::cout << "Bond number wrong!\n";
		//The hydrogen atoms
		for(int j = 0; j < 2; ++ j) {
			dx = alist_[mlist_[i]].x() - alist_[alist_[mlist_[i]].bond_list[j]].x();
			dy = alist_[mlist_[i]].y() - alist_[alist_[mlist_[i]].bond_list[j]].y();
			dz = alist_[mlist_[i]].z() - alist_[alist_[mlist_[i]].bond_list[j]].z();
			
			if(dx > L_/2.) dx -= L_;
			else if(dx < -L_/2.) dx += L_;
			if(dy > L_/2.) dy -= L_;
			else if(dy < -L_/2.) dy += L_;
			if(dz > L_/2.) dz -= L_;
			else if(dz < -L_/2.) dz += L_;
//			std::cout << "H atom " << i << " (" << dx << "," << dy << "," << dz << "," << sqrt(dx*dx + dy*dy + dz*dz) << ")\n";
			std::vector<double> atm;
			atm.push_back(1.0);
			atm.push_back(dx);
			atm.push_back(dy);
			atm.push_back(dz);
			
			ccarray.push_back(atm);			
		} 
		
		//find the wcenters
		for(size_t j = 0; j < wcenterlist_.size(); ++j ) {
			dx = alist_[mlist_[i]].x() - wcenterlist_[j].x();
			dy = alist_[mlist_[i]].y() - wcenterlist_[j].y();
			dz = alist_[mlist_[i]].z() - wcenterlist_[j].z();
			
			if(dx > L_/2.) dx -= L_;
			else if(dx < -L_/2.) dx += L_;
			if(dy > L_/2.) dy -= L_;
			else if(dy < -L_/2.) dy += L_;
			if(dz > L_/2.) dz -= L_;
			else if(dz < -L_/2.) dz += L_;

			distsq = dx*dx + dy*dy + dz*dz;
			
			//a distance criterion to decide which molecule wannier center belongs to 
			if(distsq < 1.1) {
				std::vector<double> atm;
				atm.push_back(-2.0);
				atm.push_back(dx);
				atm.push_back(dy);
				atm.push_back(dz);
//std::cout << "Molecule: " << i << " Wcenter: " << "(" << dx << "," << dy << "," << dz << "," << sqrt(dx*dx + dy*dy + dz*dz) << ")\n";
				ccarray.push_back(atm);	
				//wcs.push_back(j);
			}
			

		}

		//if(wcs.size()!=4) std::cout << "molecule charge nonezero!!!!!!!!!!!!!!!::" << wcs.size() << "\n";

		//mol_wc.push_back(wcs);
		double mmx=0., mmy=0., mmz=0.;  //dipoles for a molecule
		
		if(ccarray.size()!=6) std::cout << "Molecule charge none-zero!\n";
		for(size_t k = 0; k < ccarray.size(); ++ k) {
			mmx += ccarray[k][0]*ccarray[k][1];
			mmy += ccarray[k][0]*ccarray[k][2];
			mmz += ccarray[k][0]*ccarray[k][3];
			
		} 
//		std::cout << "Molecule: " << i << "Dipole "<< "(" << mmx << "," << mmy << "," << mmz << ") " << sqrt(mmx*mmx+mmy*mmy+mmz*mmz) << "\n" ;
		Mx += mmx;
		My += mmy;
		Mz += mmz;	
		
		daverage += sqrt(mmx*mmx+mmy*mmy+mmz*mmz);
		//std::cout << "(" << mmx << "," << mmy << "," << mmz << ")" << sqrt(mmx*mmx+mmy*mmy+mmz*mmz) << "\n";
	}
	static double result[3];
	result[0] = Mx;
	result[1] = My;
	result[2] = Mz;
	
//	std::cout << "Average dipole for molecule: " << daverage/mlist_.size() << "\n";
	return result;
}


//Find the HB stat for a step. Returns # of molecules that donate and accept 0-3 HBs
double*  MD_system::H_bond_stat(const double maxdis, const double maxangle, const char* filename) {
	double distsq, dx, dy, dz;
	//for all the O atoms.
	//Hbond donator is the one with O-H group
	std::vector<int> Ndonate;
	std::vector<int> Naccept;
	static double result[8];     // having 0 1 2 3 donors and acceptors.
	for(int i = 0; i < 8; i ++) result[i] = 0.;
	
	for(int i = 0; i < nmolecule_; ++i) {
		Ndonate.push_back(0);
		Naccept.push_back(0);
	} 
	for(size_t i = 0; i < mlist_.size(); ++ i) {

		for(size_t j = i+1; j < mlist_.size(); ++ j) {
			dx = alist_[mlist_[i]].x()-alist_[mlist_[j]].x();
			dy = alist_[mlist_[i]].y()-alist_[mlist_[j]].y();
			dz = alist_[mlist_[i]].z()-alist_[mlist_[j]].z();
									
			if(dx > L_/2.) dx -= L_;
			if(dx < -L_/2.) dx += L_;
			if(dy > L_/2.) dy -= L_;
			if(dy < -L_/2.) dy += L_;
			if(dz > L_/2.) dz -= L_;
			if(dz < -L_/2.) dz += L_;

			distsq = dx*dx + dy*dy + dz*dz;
		
			if(distsq < maxdis && distsq > 6) {
			
				//donor is i;
				
				//try both H atoms of i molecule.
				
				//first H atom.
				
				if(alist_[mlist_[i]].bond_list.size() != 2) std::cout << "error!\n";
				double vx1 = alist_[mlist_[i]].x()-alist_[alist_[mlist_[i]].bond_list[0]].x();
				double vy1 = alist_[mlist_[i]].y()-alist_[alist_[mlist_[i]].bond_list[0]].y();
				double vz1 = alist_[mlist_[i]].z()-alist_[alist_[mlist_[i]].bond_list[0]].z();
				
				double vx2 = alist_[mlist_[i]].x()-alist_[mlist_[j]].x();
				double vy2 = alist_[mlist_[i]].y()-alist_[mlist_[j]].y();
				double vz2 = alist_[mlist_[i]].z()-alist_[mlist_[j]].z();
			
				if(vx1 > L_/2.) vx1 = vx1-L_;
				if(vx1 < -L_/2.) vx1 = vx1+L_;
				if(vy1 > L_/2.) vy1 = vy1-L_;
				if(vy1 < -L_/2.) vy1 = L_+vy1;
				if(vz1 > L_/2.) vz1 = vz1-L_;
				if(vz1 < -L_/2.) vz1 = L_+vz1;
				if(vx2 > L_/2.) vx2 = vx2-L_;
				if(vx2 < -L_/2.) vx2 = L_+vx2;
				if(vy2 > L_/2.) vy2 = vy2-L_;
				if(vy2 < -L_/2.) vy2 = L_+vy2;
				if(vz2 > L_/2.) vz2 = vz2-L_;
				if(vz2 < -L_/2.) vz2 = L_+vz2; 
			
			     double angle1 = 180.*acos((vx1*vx2+vy1*vy2+vz1*vz2)/sqrt(vx1*vx1+vy1*vy1+vz1*vz1)/sqrt(vx2*vx2+vy2*vy2+vz2*vz2))/PI;
			     
				if(angle1 < maxangle) {
					Ndonate[i] ++;
					Naccept[j] ++;
				}
				
				//second H atom;
				vx1 = alist_[mlist_[i]].x()-alist_[alist_[mlist_[i]].bond_list[1]].x();
				vy1 = alist_[mlist_[i]].y()-alist_[alist_[mlist_[i]].bond_list[1]].y();
				vz1 = alist_[mlist_[i]].z()-alist_[alist_[mlist_[i]].bond_list[1]].z();
				
				if(vx1 > L_/2.) vx1 = vx1-L_;
				if(vx1 < -L_/2.) vx1 = vx1+L_;
				if(vy1 > L_/2.) vy1 = vy1-L_;
				if(vy1 < -L_/2.) vy1 = L_+vy1;
				if(vz1 > L_/2.) vz1 = vz1-L_;
				if(vz1 < -L_/2.) vz1 = L_+vz1;								
				
				double angle2 = 180.*acos((vx1*vx2+vy1*vy2+vz1*vz2)/sqrt(vx1*vx1+vy1*vy1+vz1*vz1)/sqrt(vx2*vx2+vy2*vy2+vz2*vz2))/PI;
				
				if(angle2 < maxangle) {
					Ndonate[i] ++;
					Naccept[j] ++;
				}
				
				//switch i j and again. j is the donor.
								
				//first H atom.
				vx1 = alist_[mlist_[j]].x()-alist_[alist_[mlist_[j]].bond_list[0]].x();
				vy1 = alist_[mlist_[j]].y()-alist_[alist_[mlist_[j]].bond_list[0]].y();
				vz1 = alist_[mlist_[j]].z()-alist_[alist_[mlist_[j]].bond_list[0]].z();
				
				vx2 = alist_[mlist_[j]].x()-alist_[mlist_[i]].x();
				vy2 = alist_[mlist_[j]].y()-alist_[mlist_[i]].y();
				vz2 = alist_[mlist_[j]].z()-alist_[mlist_[i]].z();
			
				if(vx1 > L_/2.) vx1 = vx1-L_;
				if(vx1 < -L_/2.) vx1 = vx1+L_;
				if(vy1 > L_/2.) vy1 = vy1-L_;
				if(vy1 < -L_/2.) vy1 = L_+vy1;
				if(vz1 > L_/2.) vz1 = vz1-L_;
				if(vz1 < -L_/2.) vz1 = L_+vz1;
				if(vx2 > L_/2.) vx2 = vx2-L_;
				if(vx2 < -L_/2.) vx2 = L_+vx2;
				if(vy2 > L_/2.) vy2 = vy2-L_;
				if(vy2 < -L_/2.) vy2 = L_+vy2;
				if(vz2 > L_/2.) vz2 = vz2-L_;
				if(vz2 < -L_/2.) vz2 = L_+vz2; 
			
			     angle1 = 180.*acos((vx1*vx2+vy1*vy2+vz1*vz2)/sqrt(vx1*vx1+vy1*vy1+vz1*vz1)/sqrt(vx2*vx2+vy2*vy2+vz2*vz2))/PI;
			     
				if(angle1 < maxangle) {
					Ndonate[j] ++;
					Naccept[i] ++;
				}
				
				//second H atom;
				vx1 = alist_[mlist_[j]].x()-alist_[alist_[mlist_[j]].bond_list[1]].x();
				vy1 = alist_[mlist_[j]].y()-alist_[alist_[mlist_[j]].bond_list[1]].y();
				vz1 = alist_[mlist_[j]].z()-alist_[alist_[mlist_[j]].bond_list[1]].z();
				
				if(vx1 > L_/2.) vx1 = vx1-L_;
				if(vx1 < -L_/2.) vx1 = vx1+L_;
				if(vy1 > L_/2.) vy1 = vy1-L_;
				if(vy1 < -L_/2.) vy1 = L_+vy1;
				if(vz1 > L_/2.) vz1 = vz1-L_;
				if(vz1 < -L_/2.) vz1 = L_+vz1;								
				
				angle2 = 180.*acos((vx1*vx2+vy1*vy2+vz1*vz2)/sqrt(vx1*vx1+vy1*vy1+vz1*vz1)/sqrt(vx2*vx2+vy2*vy2+vz2*vz2))/PI;
				
				if(angle2 < maxangle) {
					Ndonate[j] ++;
					Naccept[i] ++;
				}
				
			} //end if(distsq is good)
		}//end for j
	}//end for i
	
	for(int i = 0; i < nmolecule_; i ++) {
	//	std::cout << "(" << i << "," << Ndonate[i] << "," << Naccept[i] << ")" << "\n";
		if(Ndonate[i]<4) result[Ndonate[i]] += 1.;
		if(Naccept[i]<4) result[4+Naccept[i]] += 1.;
	}
	
	for(int j = 0; j < 8; j ++)
		result[j] = result[j]/nmolecule_;
		
	return result;
}
				
				
				
void MD_system::Show_hist(const char* outname) {
	FILE* output;
	output = fopen(outname, "w");
	for(size_t i = 0; i < Ahist_.size(); i ++) 
		fprintf(output, "%lf\t%d\n", i*step_, Ahist_[i]);
}
#endif
