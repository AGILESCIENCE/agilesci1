////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       AGILE Science Tools
//       AG testedp
//		 2018
//		 Author: Andrea Bulgarelli (INAF/IASF Bologna)
//
// NOTICE
//       Any information contained in this software
//       is property of the AGILE TEAM and is strictly
//       private and confidential.
//       Copyright (C) 2005-2019 AGILE Team. All rights reserved.
/*
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
////////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <string.h>

//#define DEBUG 1

#include <FitsUtils.h>
#include <Selection.h>
#include <Eval.h>
#include <PilParams.h>
#include <CalibUtils.h>

#include "TNamed.h"
#include "TVirtualFitter.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMinuit.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"
#include "Math/GaussLegendreIntegrator.h"

using std::cout;
using std::endl;
using std::vector;

double UpdateNormPL(double eMin, double eMax, double index)
{

	//0 - PL -> (k E^-{\index})
	//Analytical expression
	index = 1.0-index;
	double m_normFactor = (pow(eMin, index)-pow(eMax, index));
	return m_normFactor;
	/*
	 for (int i=0; i<eneChanCount-1; i++)
	 	specwt[i] = pow(double(m_energy[i]), 1.0-m_index) - pow(double(m_energy[i+1]), 1.0-m_index);
	 specwt[eneChanCount-1] = pow(double(m_energy[eneChanCount-1]), 1.0-m_index);
	 */
	/*
	 // cout << "NUMINT PL" << endl;
	 TF1 f("PowerLaw", "x^(-[0])", m_eInf, m_eSup);
	 f.SetParameter(0, index);
	 ROOT::Math::WrappedTF1 wf1(f);
	 ROOT::Math::GaussIntegrator ig;
	 ig.SetFunction(wf1);
	 ig.SetRelTolerance(0.001);
	 m_normFactor = ig.Integral(eMin, eMax) / ig.Integral(m_eInf, m_eSup);
	 */

}

double UpdateNormPLSuperExpCutOff(double eMin, double eMax, double index, double m_par2, double m_par3, double m_eInf, double m_eSup)
{
	//3 - PLSuperExpCutoff k E^-{\index} e^ ( - pow(E / E_c, gamma2) ) -> par2 = E_c, par3 = gamma2, index=gamma1
	TF1 f("PLSuperExpCutoff", "x^(-[0]) * exp(- pow(x / [1], [2]))", m_eInf, m_eSup);
	f.SetParameter(0, index);
	f.SetParameter(1, m_par2);
	f.SetParameter(2, m_par3);
	ROOT::Math::WrappedTF1 wf1(f);
	ROOT::Math::GaussIntegrator ig;
	ig.SetFunction(wf1);
	ig.SetRelTolerance(0.001);
	double m_normFactor = ig.Integral(eMin, eMax);
	return m_normFactor;
}

double UpdateNormPLExpCutOff(double eMin, double eMax, double index, double m_par2, double m_eInf, double m_eSup)
{
	//1 - PLExpCutoff -> k E^-{\index} e^ ( - E / E_c ) -> par2 = E_c
	TF1 f("PLExpCutoff", "x^(-[0]) * exp(- x / [1])", m_eInf, m_eSup);
	f.SetParameter(0, index);
	f.SetParameter(1, m_par2);
	ROOT::Math::WrappedTF1 wf1(f);
	ROOT::Math::GaussIntegrator ig;
	ig.SetFunction(wf1);
	ig.SetRelTolerance(0.001);
	double m_normFactor = ig.Integral(eMin, eMax);
	return m_normFactor;
}

double UpdateNormLogParabola(double eMin, double eMax, double index, double m_par2, double m_par3, double m_eInf, double m_eSup)
{
	TF1 f("LogParabola", "( x / [1] ) ^ ( -( [0] + [2] * log ( x / [1] ) ) )", m_eInf, m_eSup);
	f.SetParameter(0, index);
	f.SetParameter(1, m_par2);
	f.SetParameter(2, m_par3);
	ROOT::Math::WrappedTF1 wf1(f);
	ROOT::Math::GaussIntegrator ig;
	ig.SetFunction(wf1);
	ig.SetRelTolerance(0.001);
	double m_normFactor = ig.Integral(eMin, eMax);
	return m_normFactor;
}

double detCorrectionSpectraFactor(EdpGrid &edp, int iMin, int iMax, double index, double par1, double par2, double par3, int typefun) {
	cout << "--------------" << endl;
	VecF  m_edptrueenergy = edp.TrueEnergies();
	VecF  m_edpobsenergy = edp.ObsEnergies();
	VecF  m_edptheta = edp.Thetas();
	VecF  m_edpphi = edp.Phis();
	int eneChanCount = m_edptrueenergy.Dim(0);

	cout << m_edptrueenergy[iMin] << " " << m_edptrueenergy[iMax] << endl;
	double normsumpl=0, normsumple = 0;
	for (int i=iMin; i<=iMax; i++) {
		double lastenergy = m_edptrueenergy[i+1];
		cout << i << endl;
		double udp1 = 0;
		if(typefun == 0)
			udp1 = UpdateNormPL(m_edptrueenergy[i], lastenergy, par1);
		if(typefun == 1)
			udp1 = UpdateNormPLExpCutOff(m_edptrueenergy[i], lastenergy, par1, par2, m_edptrueenergy[0], m_edptrueenergy[eneChanCount-1]);
		if(typefun == 2)
			udp1 = UpdateNormPLSuperExpCutOff(m_edptrueenergy[i], lastenergy, par1, par2, par3, m_edptrueenergy[0], m_edptrueenergy[eneChanCount-1]);
		normsumple += udp1;

		udp1 = UpdateNormPL(m_edptrueenergy[i], lastenergy, index);
		normsumpl += udp1;
	}
	cout << "A " << normsumpl << " " << normsumple << endl;
	VecF edpArr(eneChanCount);
	edpArr = 0.0f;

	for(int thetaind=0; thetaind<m_edptheta.Dim(0); thetaind++) {
		for(int phiind=0; phiind<m_edpphi.Dim(0); phiind++){
			int phiindcor = phiind%2?phiind-1:phiind;
			float avgValuePL = 0.0f;
			float avgValuePLE = 0.0f;
			for (int etrue = 0; etrue < eneChanCount-1; etrue++)  {
				double lastenergy = m_edptrueenergy[etrue+1];
				//if(etrue == iMax && GetEmax() == 50000)
				//	lastenergy = 50000;
				
				for (int eobs = iMin;  eobs <= iMax; eobs++) {
					//cout << "edp: " << m_edptrueenergy[etrue] << " " << m_edpobsenergy[eobs] << " " << m_edptheta[thetaind] << " " << m_edpphi[phiindcor] << " " <<  edp.Val(m_edptrueenergy[etrue], m_edpobsenergy[eobs], m_edptheta[thetaind], m_edpphi[phiindcor]) << endl;
					edpArr[etrue] += edp.Val(m_edptrueenergy[etrue], m_edpobsenergy[eobs], m_edptheta[thetaind], m_edpphi[phiindcor]); //CORRETTO
				}
				avgValuePL += edpArr[etrue] * UpdateNormPL(m_edptrueenergy[etrue], lastenergy, index) * 1;
				//cout << UpdateNormPL(m_edptrueenergy[etrue], m_edptrueenergy[etrue+1], 2.1) << " " << edpArr[etrue] << endl;
				if(typefun == 0)
					avgValuePLE += edpArr[etrue] * UpdateNormPL(m_edptrueenergy[etrue], lastenergy, par1) * 1;
				if(typefun == 1)
					avgValuePLE += edpArr[etrue] * UpdateNormPLExpCutOff(m_edptrueenergy[etrue], lastenergy, par1, par2, m_edptrueenergy[0], m_edptrueenergy[eneChanCount-1]) * 1;
				if(typefun == 2)
					avgValuePLE += edpArr[etrue] * UpdateNormPLSuperExpCutOff(m_edptrueenergy[etrue], lastenergy, par1, par2, par3,  m_edptrueenergy[0], m_edptrueenergy[eneChanCount-1]) * 1;
			}
			double avgpl = 0, avgple = 0;
			//avgalue dipende da theta e phi
			avgpl = avgValuePL/normsumpl;
			avgple = avgValuePLE/normsumple;
			double corr = avgpl/avgple;
			cout << "A " << m_edptheta[thetaind] << " " << m_edpphi[phiindcor] << " " << avgpl << " " << avgple << " " << avgpl - avgple << " PL/" << avgpl/avgple << endl;
			return corr;
		}
	}

}

double detCorrectionSpectraFactorSimple(EdpGrid &edp, int iMin, int iMax, double index, double par1) {
	cout << "--------------" << endl;
	VecF  m_edptrueenergy = edp.TrueEnergies();
	VecF  m_edpobsenergy = edp.ObsEnergies();
	VecF  m_edptheta = edp.Thetas();
	VecF  m_edpphi = edp.Phis();
	int eneChanCount = m_edptrueenergy.Dim(0);

	cout << m_edptrueenergy[iMin] << " " << m_edptrueenergy[iMax] << endl;
	double normsumpl=0, normsumple = 0;
	for (int i=iMin; i<=iMax; i++) {
		double lastenergy = m_edptrueenergy[i+1];
		cout << i << endl;
		double udp1 = 0;

		udp1 = UpdateNormPL(m_edptrueenergy[i], lastenergy, par1);
		normsumple += udp1;

		udp1 = UpdateNormPL(m_edptrueenergy[i], lastenergy, index);
		normsumpl += udp1;
	}
	cout << "A " << normsumpl << " " << normsumple << endl;
	VecF edpArr(eneChanCount);
	edpArr = 0.0f;

	
	int thetaind=0;
	int phiind=0;
	

	
	float avgValuePL = 0.0f;
	float avgValuePLE = 0.0f;
	for (int etrue = 0; etrue < eneChanCount-1; etrue++)  {
		double lastenergy = m_edptrueenergy[etrue+1];
		//if(etrue == iMax && GetEmax() == 50000)
		//	lastenergy = 50000;
		cout << "edp: thetaind  phiind  etrue  eobs m_edptrueenergy[etrue] m_edpobsenergy[eobs] m_edptheta[thetaind] m_edpphi[phiind] edp.Val(m_edptrueenergy[etrue], m_edpobsenergy[eobs], m_edptheta[thetaind], m_edpphi[phiind])" << endl;
		for (int eobs = iMin;  eobs <= iMax; eobs++) {
			cout << "edp: " << thetaind << " " << phiind << " " << etrue << " " << eobs << " " << m_edptrueenergy[etrue] << " " << m_edpobsenergy[eobs] << " " << m_edptheta[thetaind] << " " << m_edpphi[phiind] << " " <<  edp.Val(m_edptrueenergy[etrue], m_edpobsenergy[eobs], m_edptheta[thetaind], m_edpphi[phiind]) << endl;
			edpArr[etrue] += edp.Val(m_edptrueenergy[etrue], m_edpobsenergy[eobs], m_edptheta[thetaind], m_edpphi[phiind]); //CORRETTO
		}
		avgValuePL  += edpArr[etrue] * UpdateNormPL(m_edptrueenergy[etrue], lastenergy, index) * 1;
		avgValuePLE += edpArr[etrue] * UpdateNormPL(m_edptrueenergy[etrue], lastenergy, par1) * 1;
	}
	double avgpl = 0, avgple = 0;
	//avgalue dipende da theta e phi
	avgpl = avgValuePL/normsumpl;
	avgple = avgValuePLE/normsumple;
	double corr = avgpl/avgple;
	cout << "A " << m_edptheta[thetaind] << " " << m_edpphi[phiind] << " " << avgpl << " " << avgple << " " << avgpl - avgple << " PL/" << avgpl/avgple << endl;
	return corr;


}


int mainPSF(int argc, char *argv[])
{
	char* psffilename = argv[1];
	cout << psffilename << endl;
	PsfGrid psf;
	psf.Read(psffilename);
	VecF  m_edptrueenergy = psf.Energies();
	cout << "psf energy " << m_edptrueenergy.Dim(0) << endl;
	for(int i=0; i<m_edptrueenergy.Dim(0); i++)
		cout << i+1 << " " << m_edptrueenergy[i] <<" check index " << m_edptrueenergy.GeomIndex(m_edptrueenergy[i]) <<  endl;
	return 0;
}

int mainCheckBoundaries(int argc, char *argv[]) {

	char* edpfilename = argv[3];
	cout << edpfilename << endl;
	EdpGrid edp;
	edp.Read(edpfilename);
	VecF  m_energy = edp.TrueEnergies();
	cout << "true energy " << m_energy.Dim(0) << endl;
	for(int i=0; i<m_energy.Dim(0); i++)
		cout << i+1 << " " << m_energy[i] <<" check index " << m_energy.GeomIndex(m_energy[i]) <<  endl;

	int resultMask = 0;

	double m_emin = atoi(argv[1]);
	double m_emax = atoi(argv[2]);

	int iMin = m_energy.GeomIndex(m_emin);
	if (m_emin!=m_energy[iMin])
		resultMask = resultMask | 1;	/// Using different energy lower bound

	int eneChanCount = m_energy.Size();
	int iMax = eneChanCount-1;
	if (m_emax<=m_energy[iMax]) {
		iMax = m_energy.GeomIndex(m_emax);
		if (m_emax!=m_energy[iMax])
			resultMask = resultMask | 2;	/// Using different energy upper bound
		if (iMax>iMin)
			--iMax;
		//iMax;
	}
	else
		resultMask = resultMask | 4;	/// Upper bound treated as infinity

	cout << "Boundaries: " << m_energy[iMin] <<  " " << m_energy[iMax] << " resultMask " << resultMask << endl;
	return 0;
}

int main2(int argc, char *argv[]) {

	char* edpfilename = argv[1];
	cout << edpfilename << endl;
	EdpGrid edp;
	edp.Read(edpfilename);
	VecF  m_energy = edp.TrueEnergies();
	cout << "true energy " << m_energy.Dim(0) << endl;
	for(int i=0; i<m_energy.Dim(0); i++)
		cout << i+1 << " " << m_energy[i] <<" check index " << m_energy.GeomIndex(m_energy[i]) <<  endl;
	Mat4F m_edpgrid = edp.Values();
	VecF  m_edptheta = edp.Thetas();
	//Correction factor
	int numtheta = m_edpgrid.Dim(1);
	int numphi = m_edpgrid.Dim(0);
	int eneChanCount = m_energy.Dim(0);
	TCanvas* c1 =  new  TCanvas();
	gPad->SetLogy();
	//for (int thetaind = 0; thetaind < numtheta; thetaind++) {
		//for (int phiind = 0; phiind < numphi; phiind++) {
	int thetaind = 8;
			int phiind = 0;
			double val3 = 0;

			double sum=0;
			for (int etrue = 0; etrue < eneChanCount; etrue++) {
				TString title; title += "True energy (low boundary): "; title += m_energy[etrue]; title += " MeV";
				TH1D* h1 = new TH1D(title, title, eneChanCount, 0, eneChanCount);
				TH1D* h2 = 0;
				TString name;
				name += m_edptheta[thetaind]; name += "_";
				name += phiind; name += "_";
				name += m_energy[etrue] ; name += "";
				name += ".png";
				for (int eobs = 0;  eobs < eneChanCount; eobs++) {
					val3 = m_edpgrid(phiind, thetaind, eobs, etrue);//CORRETTO
					cout << eneChanCount << " " << eobs << " " << val3 << endl;
					sum += val3;
					h1->SetBinContent(eobs+1, val3);
					if(eobs == eneChanCount-2) {
						h2 = (TH1D*) h1->Clone("edp2");
						h2->SetBinContent(eobs+1, h1->GetBinContent(eneChanCount-3+1) / 2.0);
						cout << "*" << eneChanCount << " " << eobs << " " << h2->GetBinContent(eobs+1) << endl;
					}
					if(eobs == eneChanCount-1) {

						h2->SetBinContent(eobs+1, h1->GetBinContent(eneChanCount-3+1) / 4.0);
						cout << "*" << eneChanCount << " " << eobs << " " << h2->GetBinContent(eobs+1) << endl;
						double scalefactor = h2->Integral();
						h2->Scale(1.0/scalefactor);
					}


				}
				cout << "SUM " << sum << endl;
				sum = 0;
				h1->Draw();
				h2->SetLineColor(kRed);
				h2->Draw("HIST SAME");
				c1->SaveAs(name);


			}
		//}
	//}

}

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
	int eneChanCount = m_edptrueenergy.Dim(0);
	for (int i=0; i<eneChanCount-1; i++) {
		double udp1 = UpdateNormPLExpCutOff(m_edptrueenergy[i], m_edptrueenergy[i+1], 1.57, 1678, m_edptrueenergy[0], m_edptrueenergy[m_edptrueenergy.Dim(0)-1]);
		cout << m_edptrueenergy[i] << " " << m_edptrueenergy[i+1] << " " << UpdateNormPL(m_edptrueenergy[i], m_edptrueenergy[i+1], 2.1) << " " << udp1 << " " << UpdateNormPL(m_edptrueenergy[i], m_edptrueenergy[i+1], 2.1)  -  udp1 << endl;
	}

	int i1=4, i2=8;
	double udp1 = UpdateNormPLExpCutOff(m_edptrueenergy[i1], m_edptrueenergy[i2], 1.57, 1678, m_edptrueenergy[0], m_edptrueenergy[m_edptrueenergy.Dim(0)-1]);
	cout << m_edptrueenergy[i1] << " " << m_edptrueenergy[i2] << " " << UpdateNormPL(m_edptrueenergy[i1], m_edptrueenergy[i2], 2.1) << " " << udp1 << " " << UpdateNormPL(m_edptrueenergy[i1], m_edptrueenergy[i2], 2.1)  -  udp1 << endl;
	i1=8; i2=10;
	udp1 = UpdateNormPLExpCutOff(m_edptrueenergy[i1], m_edptrueenergy[i2], 1.57, 1678, m_edptrueenergy[0], m_edptrueenergy[m_edptrueenergy.Dim(0)-1]);
	cout << m_edptrueenergy[i1] << " " << m_edptrueenergy[i2] << " " << UpdateNormPL(m_edptrueenergy[i1], m_edptrueenergy[i2], 2.1) << " " << udp1 << " " << UpdateNormPL(m_edptrueenergy[i1], m_edptrueenergy[i2], 2.1)  -  udp1 << endl;
	i1=10; i2=12;
	udp1 = UpdateNormPLExpCutOff(m_edptrueenergy[i1], m_edptrueenergy[i2], 1.57, 1678, m_edptrueenergy[0], m_edptrueenergy[m_edptrueenergy.Dim(0)-1]);
	cout << m_edptrueenergy[i1] << " " << m_edptrueenergy[i2] << " " << UpdateNormPL(m_edptrueenergy[i1], m_edptrueenergy[i2], 2.1) << " " << udp1 << " " << UpdateNormPL(m_edptrueenergy[i1], m_edptrueenergy[i2], 2.1)  -  udp1 << endl;

	//PROCEDURE FROM HERE
	//input iMin, iMax, indexPL, index1, par2, par3, funtype

	double par1 = 1.6692;
	double par2 = 3403; //ec
	double par3 = 0;
	double index = 2.1; //index, gamma1
	int typefun = 0;
	//int iMin = 10; int iMax=12;
	detCorrectionSpectraFactor      (edp, 4, 12, index, par1, par2, par3, typefun);
	detCorrectionSpectraFactorSimple(edp, 4, 12, index, par1);
	/*
	detCorrectionSpectraFactor(edp, 4, 6-1, index, par1, par2, par3, typefun);
	detCorrectionSpectraFactor(edp, 6, 8-1, index, par1, par2, par3, typefun);
	detCorrectionSpectraFactor(edp, 8, 10-1, index, par1, par2, par3, typefun);
	detCorrectionSpectraFactor(edp, 10, 12-1, index, par1, par2, par3, typefun);
	detCorrectionSpectraFactor(edp, 12, 13, index, par1, par2, par3, typefun);
	detCorrectionSpectraFactor(edp, 3, 12, index, par1, par2, par3, typefun);
	 */
	return 0;
	//for(int i_true=0; i_true<m_edptrueenergy.Dim(0); i_true++)
	int i_true=2-1; //true energy
	double sum = 0;
	double val, val2, val3;
	for(int j=0; j<m_edpobsenergy.Dim(0); j++) {

		// CURRENT m_edp.Val(m_energy[etrue], m_energy[eobs], m_theta[thetaind], m_phi[phiind]);

		val = edp.Val(m_edptrueenergy[i_true], m_edpobsenergy[j], 20, 270); //CORRETTO
		int thetain = 0;
		int phiin = 0;
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
	sum=0;
	int numtheta = m_edpgrid.Dim(1);
	int numphi = m_edpgrid.Dim(0);
	//int eneChanCount = m_edpgrid.Dim(2);
	for (int thetaind = 0; thetaind < numtheta/2; thetaind++) {
		//for (int phiind = 0; phiind < numphi; phiind++) {
		int phiind = 0;
			for (int etrue = 0; etrue < eneChanCount; etrue++) {
				for (int eobs = 0;  eobs < eneChanCount; eobs++) {
					val = edp.Val(m_edptrueenergy[etrue], m_edpobsenergy[eobs], m_edptheta[thetaind], m_edpphi[phiind]);
					cout << "EDP VALUE: " << etrue << " " << m_edptrueenergy[etrue] << " " << eobs << " " << m_edpobsenergy[eobs] << " " << thetaind << " " << m_edptheta[thetaind] << " " << phiind << " " << m_edpphi[phiind] << " " << val << endl;
					sum += val;
					//avgValue += edpArr[etrue] * specwt[etrue] * m_aeffgrid(phiind, thetaind, etrue);
				}
				cout << "SUM" << sum << endl;
				sum = 0;
			}
		//m_avgValues(phiind, thetaind) = avgValue/normsum;
		//}
	}


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
