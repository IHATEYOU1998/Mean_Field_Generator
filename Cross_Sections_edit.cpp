#include <iostream>
#include <cmath>
#include <cstdlib>
#include <math.h>

#include "TVector3.h"

#include "constants.h"
#include "Cross_Sections.h"

// Given to me by Andrew Denniston on April 5th 2019 (friday)
// Changes made since then: deleted helpers.h (with funtion sq). Changed all sq to pow(,2) to sqaure it

Cross_Sections::Cross_Sections()
{
  // Set defaults
  myModel=kelly;
  myMethod=cc1;
}

Cross_Sections::Cross_Sections(csMethod thisMeth, ffModel thisMod)
{
  std::cerr << "Cross_Sections: you have selected configuration: " << thisMeth << " " << thisMod <<"\n";
  myModel=thisMod;
  myMethod=thisMeth;
}

Cross_Sections::~Cross_Sections(){}

double Cross_Sections::sigma_eN(double Ebeam,TVector3 k, TVector3 p, bool isProton)
{
  switch (myMethod)
    {
    case onshell:
      return sigma_onShell_by_Etheta(Ebeam,k,isProton);
    case cc1:
      return sigmaCC1(Ebeam,k,p,isProton);
    case cc2:
      return sigmaCC2(Ebeam,k,p,isProton);
    default:
      {
	std::cerr << "Invalid cross section method! Double check and fix!\n";
        exit(-1);
      }
    }
  return 0;
}

// Because DeForest uses the opposite Lorentz convention
double dot4(double x0, TVector3 x, double y0, TVector3 y)
{
  return ((x0*y0)-(x*y));
}

double Cross_Sections::sigmaCCn(double Ebeam, TVector3 k, TVector3 p, bool isProton, int n)

{
  TVector3 q = TVector3(0.,0.,Ebeam) - k;
  TVector3 pM = p-q;
  
  double omega = Ebeam - k.Mag();
  double QSq = q.Mag2() - pow(omega,2);
  double E = sqrt(p.Mag2() + pow(mN,2));
  double Ebar = sqrt(pM.Mag2() + pow(mN,2));
  double omegabar = E-Ebar;
  double QSqbar = q.Mag2() - pow(omegabar,2);

  // Calculate form factors
  double GE = (isProton)? GEp(QSq) : GEn(QSq);
  double GM = (isProton)? GMp(QSq) : GMn(QSq);

  double F1 = (GE + GM * QSq/(4.*pow(mN,2)))/(1. + QSq/(4.*pow(mN,2)));
  double kF2 = (GM - GE)/(1. + QSq/(4.*pow(mN,2)));

  double wC;
  double wT;
  double wS;
  double wI;
  
  if (n==1)
    {
      wC = (pow(E+Ebar,2)*(pow(F1,2) + QSqbar/(4.*mN*mN) * pow(kF2,2)) - q.Mag2()*pow(F1 + kF2,2))/(4.*E*Ebar);
      wT = QSqbar*pow(F1 + kF2,2)/(2.*Ebar*E);
      wS = p.Mag2() * pow(sin(p.Angle(q)),2) * (pow(F1,2) + QSqbar/(4.*mN*mN) * pow(kF2,2))/(E*Ebar);
      wI = -p.Mag()*sin(p.Angle(q))*(Ebar + E)*(pow(F1,2) + QSqbar/(4.*mN*mN) * pow(kF2,2))/(E*Ebar);
    }
  else if (n==2)
    {  
      double pbarp = dot4(Ebar,pM,E,p);
      double pbarq = dot4(Ebar,pM,omega,q);
      double pq = dot4(E,p,omega,q);
      double qbarq = dot4(omegabar,q,omega,q);
      double sumq = dot4((Ebar+E),(pM+p),omega,q);

      wC = (E*Ebar
	    + 0.5 * (pbarp + pow(mN,2)) * pow(F1,2)
		   - 0.5 * q.Mag2() * F1 * kF2
		   - ((pbarq*E + pq*Ebar)*omega
		      - Ebar * E * QSq
		      + pbarq * pq
		      - 0.5 * (pbarp - pow(mN,2)) * q.Mag2())
	    * pow(kF2,2)/(4*pow(mN,2))
		   )/(E*Ebar);
      wT = (-(pbarp + pow(mN,2)) * pow(F1,2)
		   + qbarq * F1 * kF2
		   + (2*pbarq*pq
		      - (pbarp - pow(mN,2))*QSq)
	    * pow(kF2,2)/(4*pow(mN,2))
		   )/(Ebar*E);
      wS = p.Mag2() * pow(sin(p.Angle(q)),2) * (pow(F1,2)
						+ QSq/(4.*mN*mN) * pow(kF2,2))/(E*Ebar);
      wI = p.Mag()*sin(p.Angle(q))*(-(Ebar + E) * pow(F1,2)
					   + (sumq * omega
					      - (Ebar + E) * QSq)
				    * pow(kF2,2)/(4*pow(mN,2))
					   )/(E*Ebar);
    }
  else
    {
      std::cerr << "Invalid cross section designation. Check and fix. Exiting\n\n\n";
    }
      
  double sigmaMott = cmSqGeVSq * 4. * pow(alpha,2) * k.Mag2() * pow(cos(k.Theta()/2.),2) / pow(QSq,2);

  double phi = q.Cross(k).Angle( q.Cross(p) );
  return sigmaMott * ( pow(QSq,2)/q.Mag2() * wC +
                       (QSq/(2.*q.Mag2()) + pow(tan(k.Theta()/2.),2)) * wT +
                       QSq/q.Mag2() * sqrt(QSq/q.Mag2() + pow(tan(k.Theta()/2.),2)) * wI * cos(phi) +
                       (QSq/q.Mag2() * pow(cos(phi),2) + pow(tan(k.Theta()/2.),2)) * wS
                       );
}

double Cross_Sections::sigmaCC1(double Ebeam, TVector3 k, TVector3 p, bool isProton)
{
  return sigmaCCn(Ebeam, k, p, isProton, 1);
}

double Cross_Sections::sigmaCC2(double Ebeam, TVector3 k, TVector3 p, bool isProton)
{
  return sigmaCCn(Ebeam, k, p, isProton, 2);
}

double Cross_Sections::sigma_onShell_by_Etheta(double Ebeam, TVector3 k, bool isProton)
{
  double theta=k.Theta();
  double E3 = Ebeam * mN/ (mN + Ebeam*(1.-k.CosTheta()));
  double QSq = 2. * Ebeam * E3 * (1.-k.CosTheta());
  double tau = QSq/(4.*mN*mN);
  double GE = isProton ? GEp(QSq) : GEn(QSq);
  double GM = isProton ? GMp(QSq) : GMn(QSq);
  double epsilon = epsilon = 1./(1.+2.*(1.+tau)*pow(tan(theta/2.),2));

  double sigmaMott = cmSqGeVSq * pow(2.*alpha*E3 * cos(theta/2.)/QSq,2) * (E3/Ebeam);

  double sigmaRosenbluth = sigmaMott * (pow(GE,2) + tau/epsilon * pow(GM,2))/(1. + tau);
  return sigmaRosenbluth * Ebeam / (E3 * (2.*tau + 1.));
}

double Cross_Sections::GEp(double QSq)
{
  switch (myModel)
    {
    case dipole:
      return Gdipole(QSq);
    case kelly:
      return Gkelly(QSq,-0.24,10.98,12.82,21.97);
    default:
      std::cerr << "Error in GEp: invalid form factor model!\n";
      exit(-1);
    }
  return 0.;
}

double Cross_Sections::GEn(double QSq) // This will use the Galster parameterization
{
  double tau = QSq/(4.*mN*mN);
  return 1.70 * tau / (1. + 3.3 * tau) * Gdipole(QSq); // params from Kelly paper
  //return mu_n * tau / (1. + 5.6 * tau) * Gdipole(QSq); // the original Galster numbers
}

double Cross_Sections::GMp(double QSq)
{
  switch (myModel)
    {
    case dipole:
      return mu_p * Gdipole(QSq);
    case kelly:
      return mu_p * Gkelly(QSq,0.12,10.97,18.86,6.55);
    default:
      std::cerr << "Error in GMp: invalid form factor model!\n";
      exit(-1);
    }
  return 0.;
}

double Cross_Sections::GMn(double QSq)
{
  switch (myModel)
    {
    case dipole:
      return mu_n * Gdipole(QSq);
    case kelly:
      return mu_n * Gkelly(QSq,2.33,14.72,24.20,84.1);
    default:
      std::cerr << "Error in GMn: invalid form factor model!\n";
      exit(-1);
    }
  return 0.;
}

double Cross_Sections::Gdipole(double QSq){ return 1. / pow(1 + QSq/0.71,2); }

double Cross_Sections::Gkelly(double QSq,double a1, double b1, double b2, double b3)
{
  double tau = QSq/(4.*mN*mN);
  double denom = 1. + b1*tau + b2*tau*tau + b3*tau*tau*tau;
  double numer = 1. + a1*tau;
  return numer/denom;
}