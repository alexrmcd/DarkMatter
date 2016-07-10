////////////////////////////////////////////////////
////		FluxDensity.cpp 			5-19-16  ///
////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <cmath>
#include <ctime>
#include <vector>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <iomanip> 


#include "Constants.h"

/* to compile and run, use command 

g++ -o NSD_IC_flux NSD_IC_flux.cpp  -I/home/alex/research/darksusy-5.1.2/include -L/home/alex/research/darksusy-5.1.2/lib -lgsl -lgslcblas -ldarksusy -lFH -lHB -lgfortran
*/

//////////////////////////// routines for calling fortran/darksusy stuff /////////////////////////

extern"C" {	 								//interface with fortran code to initialize darksusy

void dsinit_();

}	

extern"C" {									// dshayield gives injection spectrum
double dshayield_(double *mwimp, double *emuthr, int *ch,  int *yieldk, int *istat);
}


//////////////////////////////////////////////////////////////////////////////////////////////////


class Cluster{  
	
	public:
	std::string name;
	double z ;
	double rh; //halo radius

	double B0 ; 				//microGauss
	double Bmu;
	double rcore ; 		//0.291; //core radius in Mpc for coma
	double root_dv;

	double alpha;		// D(E) ~ E^alpha

	Cluster() :	name(""),
				z(0.0232), 						//redshift
				rh(1),						//halo radius Mpc
				B0(4.7), rcore(0.291), Bmu(1.2)		//Bfield Params
			//DIffusion parameter, not really a cluster thing but easy access is good

	{	//everything labelled or Coma, sgould work in user options
		std::cout << "creating cluster... " << std::endl;
	}

	double bfield_model (double r) {

		double beta = 0.75;
		double eta = 0.5;
		double rc = rcore*mpc2cm;
		double B_field = B0 * pow(( 1 + r*r/(rc*rc)),(-1.5*beta*eta));		// Storm et al 2013 

		return B_field;

	};	

	double bloss(double E , double r){

		double bloss = 0.0253*pow(bfield_model(r), 2)*E*E 					//bsyn bfield_model(r)
						+ 0.265 * pow(1 + z, 4 )*E*E  					//bIC
						+ 1.51*n*(0.36 + log(E/me/n) )			//+ 1.51*n*(0.36 + log(E/me) )*E						//bbrem Note this is most likely incorrect, no factor of E, and should be E/me/nin log
						+ 6.13*n*( 1 + log(E/me/n)/75); 					//bcoul

		bloss = bloss *1e-16;	

		return bloss;
	};


	double bloss(double E ){ 												//overload so that we can just have bloss(E), could also have same as b(E, r) but set r=0 ??  kinda sloppy
		double ne = 1.3e-3;
		//double Bmu = 10;													//should us average or something later

		double bloss = 0.0254*pow(Bmu, 2.0)*E*E 								//bsyn bfield_model(r)
						+ 0.25 /** pow(1 + c.z, 4 )*/*E*E  					//bIC
						+ 1.51*ne*(0.36 + log(E/me/ne) )						//brem , in Emma's code  has + 1.51*n*(0.36 + log(E/me) )*E
						+ 6.13*ne*( 1 + log(E/me/ne)/75); 					//bcoul

		bloss = bloss *1e-16;	
		//std::cout <<Bmu<<std::endl;
		return bloss;
	};


	double DM_profile(double r){
		double rhos = 0.039974 ;//DM char density in Gev /cm^3 Storm13 (Mathematica)
		double rs = 0.404*mpc2cm;  //DM scale radius in Mpc Storm13  (Mathematica)

		double rho_N04 =0.08296;		//use Colafrancesco 2006 Eq 5 and Eq 33 
		double rs_N04 = 0.28*mpc2cm; 

		double rc = rcore * mpc2cm;
		//NFW
		//double rho = rhos / ( r/rs * pow(1 + r/rs , 2 )); //x? + 1e-100


		// N04
		double rho = (rho_N04) * exp(-2.0/0.17 * ( pow(r/rs_N04, 0.17) - 1 ));


		return rho;

	};




};

class Particle {
	public:

	int ch;
	double mx;
	double sv;

	Particle() : ch(0) , mx(0), sv(4.7e-25){}

};

Cluster c;
Particle p;




////////////////////////////////////////////////

double ddist(double z  , void * params){

	double distint =  mpc2cm * clight / ( H0 * sqrt( OmegaM * pow(1 + z , 3)  + OmegaL ) ); 

	return distint;

}

//distance as a function of redshift
double Dist(){ 

	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;

	gsl_function F;
	F.function = &ddist;

	gsl_integration_qags (&F, 0, c.z, 0, 1e-3, 1000,
	                w, &result, &error); 

	gsl_integration_workspace_free (w);

	return result;
 		
}


double darksusy (double Ep){

	int yieldk = 151;
	int istat;

	double ds = dshayield_(&p.mx, &Ep, &p.ch, &yieldk, &istat);

	return ds;
}

double ddiffusion(double Ep, void * params){

	double ddiffusion = darksusy(Ep) ;

	return ddiffusion;

}

double gslInt_diffusion( double E){			// int over Ep


	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;


	gsl_function F;
	F.function = &ddiffusion;
								//pass Ep to rootdv(), pass r from dndeeq as well, 
	gsl_set_error_handler_off();
	gsl_integration_qags (&F, E, p.mx, 0, 1e-2, 1000, 
	                    w, &result, &error); 

	gsl_integration_workspace_free (w);
		//std::cout << "result: " <<result << std::endl;
	return result;

}

double dndeeq(double E, double r ){

	double dndeeq =(1 / c.bloss(E))* gslInt_diffusion(E);

	return dndeeq;

}

double G(double eps, double E_gamma, double boost){

	//double boost = E/(me * pow( clight , 2))
	double GAMMA = 4 * eps * boost/me;
	double q = E_gamma / (GAMMA * (boost * me  - E_gamma) );

	double G = 2 * q* log(q) + (1+2*q)*( 1 - q ) + pow(GAMMA*q , 2)*(1-q) / ( 2*(1+GAMMA*q ) );//
	return G;
}




double IC_cross_Section( double eps, double E_gamma, double E){


	double boost = E/me;

	double sigma_thompson = 6.6524e-35; //in km^2

	double sigma = 3 * sigma_thompson/(4* eps * pow(boost ,2))* G(eps, E_gamma, boost);

//if (sigma<0)
	//std::cout << "sigma( eps = " << eps <<", "<< "E_gamma = " << E_gamma << " , E = "<<E<<" ) = "<< sigma<<std::endl;
	return sigma;
}

double CMB_bbSpectrum(double eps){
	//double hplanckGeV = J2Gev*hplanck;
	//kb *= J2Gev;
	double T = 2.73;

	double nu = eps/(hplanck*J2Gev);

	double CMB_bbSpectrum = 8*pi* pow(nu, 2.0)/pow(clight, 3.0)/(	exp(hplanck * nu /(kb * T )) - 1	);

	
	//if (CMB_bbSpectrum == 0)
	//std::cout << CMB_bbSpectrum << "  , nu("<< eps <<")  = " << nu <<  " denom = " << 	exp(hplanck * nu /(kb * T )) -1<< std::endl;;//CMB_bbSpectrum  << std::endl;

	return CMB_bbSpectrum;
}

double dpIC(double eps, void * params ){

	//double nu = 0.02;	//Ghz
	std::vector<double> pICParams = *(std::vector<double> *)params;
	double E_gamma = pICParams[0] ;
	double E = pICParams[1] ;
	
	double q = pow(me,2)*E_gamma / ( 4 * eps * E * (E  - E_gamma) );
	
	double dpIC;
	if( q>1 or q < 1/(4* pow(E/me , 2) ) or eps/hplanck/J2Gev > 1e13){
		dpIC = 0;
		//std::cout << eps/hplanck/J2Gev  <<std::endl;
	}
	else{
		dpIC = CMB_bbSpectrum(eps) *  IC_cross_Section(eps, E_gamma, E); 
	};
	return dpIC;
}


double gslInt_pIC(double nu, double E){			//int over eps

		
	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (5000);

	double result, error;
	size_t size;

	std::vector<double> pICParams (2);

	double E_gamma = hplanck * nu * J2Gev;
	double boost = E/me;
	pICParams[0] = E_gamma;
	pICParams[1] = E;

	//if(E_gamma >= E )
	//	result = 0;
	//else{
	double eps_max = E_gamma*E/(E - E_gamma );//1.24e-12;//
	double eps_min = E_gamma/(4 * boost/me*( E - E_gamma)); //1.24e-15;//
	if(E<E_gamma)
		std::cout << eps_min/hplanck/J2Gev << "  " << eps_max/hplanck/J2Gev <<"--> " <<E << " " <<E_gamma<<std::endl;
	
	if (eps_max > 1e13*hplanck*J2Gev)
		eps_max = 1e13*hplanck*J2Gev;


	gsl_function F;
	F.function = &dpIC;
	F.params = &pICParams;
	gsl_set_error_handler_off();
	gsl_integration_qng (&F,  eps_min, eps_max,  0, 1e-3, &result, &error, &size); 

	gsl_integration_workspace_free (w);

	result *= clight* E_gamma; 
	
	
	//if(result != 0)
	//std::cout << "pIC( nu = " << nu <<", "<< "E_gamma = " << E_gamma << " , E = "<<E<<" boost = "<<boost << " ) => " << eps_min/hplanck/J2Gev << " " << eps_max/hplanck/J2Gev <<std::endl;
	//std::cout << "pIC( "<<E_gamma<<") = "<< result<<std::endl;
	return result;

}


double djIC(double E , void * params){

	std::vector<double> jICParams = *(std::vector<double> *)params;
	double nu = jICParams[0];
	double r = jICParams[1];

	double djIC = 2* gslInt_pIC(nu, E)* dndeeq(E , r);
	//if(djIC == 0)
	//std::cout << djIC<<"  pIC ("<< nu <<")= " << gslInt_pIC(nu, E) << "  dndeeq = " << dndeeq(E,r) <<std::endl;

	return djIC;
}


double gslInt_jIC(double nu, double r){ 				// int over E

			///////////
	std::clock_t start;
	double duration;
	start = std::clock();
	int a ; 
		///////////


	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;
	std::vector<double> jICParams (2);

	jICParams[0] = nu;
	jICParams[1] = r;

	gsl_function F;
	F.function = &djIC;
	F.params = &jICParams;

	gsl_integration_qags (&F, nu*hplanck*J2Gev , p.mx, 0, 1e-2, 1000,
	                    w, &result, &error); 


	gsl_integration_workspace_free (w);
	
	//std::cout << " . ";
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	//std::cout << me << " " <<  p.mx << std::endl;
	//std::cout << "jIC( r = " << r/mpc2cm*1000 <<" ) "<< result <<" -> duration: " << duration <<std::endl;  //      ~30s-60s	

	return result;

}


double dsIC( double r, void * params ){

	double nu = *(double *)params;


	double dist_z = Dist() / (1+c.z);

	double sICIntegrand = 4 *pi /pow(dist_z , 2) *pow(r,2)  *  pow(c.DM_profile(r) , 2)*gslInt_jIC(nu, r);
//if ( sICIntegrand < 0)
	//	std::cout << "sICIntegrand(r = "<< r/mpc2cm*1000<< ") = " <<sICIntegrand<<std::endl;
	
	return sICIntegrand;
}






double gslInt_sIC( double nu, double r ){				// int over r

	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;



	gsl_function F;
	F.function = &dsIC;
	F.params = &nu;

	gsl_integration_qags (&F, 1e-16, r, 0, 1e-1, 1000,
	                w, &result, &error); 






	//std::cout <<"result = " << result <<", ";
	result *=  0.00160218 *p.sv/(8* pi*pow( p.mx , 2.0 ));
		//std::cout <<"new result  = " << result <<", nu = "<<nu<<std::endl;
	return result;

}







void runFlux(double mx, double r){ //runs through values of root_dv
	
	p.mx = mx;
	int n_nu = 35;


	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();
	int a ; 
	///////before algorithm

	std::ostringstream makefilename;
	makefilename << "IC_SED_NSD" << p.mx << "Gev_coma.txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());

	for (int i = 0 ; i < n_nu + 1; ++i  ){
		//total time timer start
		std::clock_t start;
		double duration;
		start = std::clock();
		///////before algorithm

		double nu_min = 1e12;
		double nu_max = 1e24;
	
		double nu = nu_min * ( exp(    (log(nu_max) - log(nu_min))/ n_nu * i));
		double data = gslInt_sIC( nu , r)*nu ;
		file << log10(nu) << "\t" <<  data <<std::endl;

		duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
		std::cout << i<< "/" << n_nu << " " << log10(nu) << "\t" << data <<  " --> " << duration << std::endl;

	};
	
	////////after algorithm
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "Total time:  " << duration <<std::endl;
}






//NOTE that we need to set a cross section for this, use eq 1 in STorm for Snu.

main(){
	dsinit_(); //initialixe DarkSUSY
	

	double r = c.rh*mpc2cm;
	p.ch = 25;				
	p.mx = 40;			//darksusy channel
	//std::cout << gslInt_sIC( pow(10, 17.5) , r)*pow(10, 17.5) << std::endl;
	//std::cout << c.bfield_model(r) << std::endl;
 	runFlux( 40, r );

	p.ch = 13;	
	p.sv = 8.8e-26;	
	c.Bmu = 8;					//darksusy channel
	//std::cout << c.bfield_model(r) << std::endl;
	runFlux( 81, r ) ;



}




