/*
 * Copyright (c) 2005-2016
 *     Andrew Chen, Alberto Pellizzoni, Alessio Trois (IASF-Milano),
 *     Andrea Bulgarelli, Andrea Zoli (IASF-Bologna),
 *     Tomaso Contessi (Nuove Idee sas)
 *
 * Any information contained in this software
 * is property of the AGILE TEAM and is strictly
 * private and confidential.
*/

#include <iostream>
#include <string.h>

//#define DEBUG 1

#include <FitsUtils.h>
#include <Selection.h>
#include <Eval.h>
#include <PilParams.h>
#include <CalibUtils.h>

using std::cout;
using std::endl;
using std::vector;



int main(int argc, char *argv[])
{
	char* edpfilename = argv[1];
	cout << edpfilename << endl;
	EdpGrid edp;
	edp.Read(edpfilename);
	VecF  m_edptrueenergy = edp.TrueEnergies();
	cout << "true energy " << m_edptrueenergy.Dim(0) << endl;
	for(int i=0; i<m_edptrueenergy.Dim(0); i++)
		cout << i+1 << " " << m_edptrueenergy[i] <<" check index " << m_edptrueenergy.GeomIndex(m_edptrueenergy[i]) <<  endl;
	VecF  m_edpobsenergy = edp.ObsEnergies();
	cout << "m_edpobsenergy " << m_edpobsenergy.Dim(0) << endl;
	for(int i=0; i<m_edpobsenergy.Dim(0); i++)
		cout << i << " " << m_edpobsenergy[i] << " check index " << m_edpobsenergy.GeomIndex(m_edpobsenergy[i]) << endl;
	
	VecF  m_edptheta = edp.Thetas();
	cout << "theta POL (SLICE)" << m_edptheta.Dim(0) << endl;
	for(int i=0; i<m_edptheta.Dim(0); i++)
		cout << i+1 << " " << m_edptheta[i] << " check index: " << m_edptheta.LinearIndex(m_edptheta[i]) << " " << i << endl;
	
	VecF  m_edpphi = edp.Phis();
	cout << "phi AZIMUT (DATACUBE)" << m_edpphi.Dim(0) << endl;
	for(int i=0; i<m_edpphi.Dim(0); i++)
		cout << i+1 << " " << m_edpphi[i] << " check index: " << m_edpphi.LinearIndex(m_edpphi[i]) << " " << i << endl;
	
	Mat4F m_edpgrid = edp.Values();
	//DIM0 -> PHI
	//DIM1 -> THETA
	//DIM2 -> OBS
	//DIM3 -> TRUE
	cout << "Dimension of m_edpgrid: " <<m_edpgrid.Dim(0) << " " << m_edpgrid.Dim(1) << " " << m_edpgrid.Dim(2) << " " << m_edpgrid.Dim(3)  << endl;
	
	cout << "m_edptrueenergy " << m_edptrueenergy.Dim(0) << endl;
	for(int i=0; i<m_edptrueenergy.Dim(0); i++)
		cout << i << " " << m_edptrueenergy[i] << endl;
	
	//for(int i_true=0; i_true<m_edptrueenergy.Dim(0); i_true++)
	int i_true=2-1; //true energy
	double sum = 0;
	double val, val2, val3;
	for(int j=0; j<m_edpobsenergy.Dim(0); j++) {

		// CURRENT m_edp.Val(m_energy[etrue], m_energy[eobs], m_theta[thetaind], m_phi[phiind]);
		
		val = edp.Val(m_edptrueenergy[i_true], m_edpobsenergy[j], 20, 270); //CORRETTO
		int thetain = 5 - 1;
		int phiin = 7-1;
		//int l = m_edptrueenergy.GeomIndex(m_edptrueenergy[i_true]);//16
		int l = 1;
		//cout << "start#######" << endl;
		int m = m_edpobsenergy.GeomIndex(m_edpobsenergy[j]);//16
		//cout << "end####### " << m << endl;
		int n = m_edptheta.LinearIndex(m_edptheta[thetain]);//19
		int p = m_edpphi.LinearIndex(m_edpphi[phiin]);//8
		cout << i_true << " " << l << " / " << j << " " << m << endl;
		/*
		 NAXIS1  =                   16 /
		 NAXIS2  =                   16 /
		 NAXIS3  =                   19 /
		 NAXIS4  =                    8 /
		 CTYPE1  = 'TRUE-TAB'
		 CUNIT1  = 'MeV     '
		 CTYPE2  = 'OBSS-TAB '
		 CUNIT2  = 'chan    '
		 CTYPE3  = 'POL--TAB'
		 CUNIT3  = 'deg     '
		 CTYPE4  = 'AZIM-TAB'
		 CUNIT4  = 'deg
		 -> indici al contrario rispetto alle matrici in memoria
		 */
		//phi-azimut, theta-pol, obs, true
		val3 = m_edpgrid(phiin, thetain, j, i_true);//CORRETTO
		val2 = m_edpgrid(p, n, m, l); //CORRETTO
		
		cout << m_edpobsenergy[j] << " " << m_edptrueenergy[i_true] << " " << val << " " << val2 << endl;
		sum += val;
	}
	cout << sum << endl;
    return 0;
	/*
	 NAXIS1  =                   16 /
	 NAXIS2  =                   19 /
	 NAXIS3  =                    8 /
	 CTYPE1  = 'ENRG-TAB  '
	 CUNIT1  = 'MeV     '
	 CTYPE2  = 'POL--TAB'
	 CUNIT2  = 'deg     '
	 CTYPE3  = 'AZIM-TAB'
	 CUNIT3  = 'deg     '
	 m_aeffgrid(phiind, thetaind, etrue);
	 */
}
