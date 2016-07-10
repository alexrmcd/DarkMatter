/////NSD_Flux_synch.cpp///////

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

g++ -o NSD_Flux_synch NSD_Flux_synch.cpp  -I/home/alex/research/darksusy-5.1.2/include -L/home/alex/research/darksusy-5.1.2/lib -lgsl -lgslcblas -ldarksusy -lFH -lHB -lgfortran
*/

//////////////////////////// routines for calling fortran/darksusy stuff /////////////////////////

extern"C" {	 								//interface with fortran code to initialize darksusy

void dsinit_();

}	

extern"C" {									// dshayield gives injection spectrum
double dshayield_(double *mwimp, double *emuthr, int *ch,  int *yieldk, int *istat);
}


//////////////////////////////////////////////////////////////////////////////////////////////////


// consider putting this in its own header?? problem is that we need to access and change member variables
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
		//double Bmu = 1;													//should us average or something later

		double bloss = 0.0254*pow(Bmu, 2.0)*E*E 								//bsyn bfield_model(r)
						+ 0.25 /** pow(1 + c.z, 4 )*/*E*E  					//bIC
						+ 1.51*ne*(0.36 + log(E/me/ne) )						//brem , in Emma's code  has + 1.51*n*(0.36 + log(E/me) )*E
						+ 6.13*ne*( 1 + log(E/me/ne)/75); 					//bcoul

		bloss = bloss *1e-16;	

		return bloss;
	};


	double DM_profile(double r){
		double rho = 0.039974 ;//DM char density in Gev /cm^3 Storm13 (Mathematica)
		double rs = 0.404*mpc2cm;  //DM scale radius in Mpc Storm13  (Mathematica)

		double rho_N04 =0.08296;		//use Colafrancesco 2006 Eq 5 and Eq 33 
		double rs_N04 = 0.28*mpc2cm; 

		double rc = rcore * mpc2cm;
		//NFW
		//double rho = rhos / ( r/rs * pow(1 + r/rs , 2 )); //x? + 1e-100


		// N04
		double rho = (rho_N04) * exp(-2.0/0.17 * ( pow(r/rs_n04, 0.17) - 1 ));


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

//SYnchrotron emmission spectral function form Cola2006.
double fff(double x){

	double fff = 1.25 * pow( x , 1.0/3.0) * exp(-x) * pow((648 + x*x) , 1.0/12.0);

	return fff;
}

double dpsyn(double theta, void * params ){

	//double nu = 0.02;	//Ghz
	std::vector<double> psynParams = *(std::vector<double> *)params;
	double E = psynParams[0];
	double r = psynParams[1];
	double nu = psynParams[2];


	double psyn0 = 1.46323e-25 ; // Gev/s/Hz
	double x0 = 62.1881 ;			// dimensionless constant
	double nu_em =  nu/1000.0; //( 1 + c.z )* (observing freq)*(1+z)/1000 convert from Mhz to GHz 


	//std::cout << E << " , " << r/mpc2cm*1000 << std::endl;

	double x = x0 *nu_em / ( c.bfield_model( r ) * pow( E , 2) );
	double xc = 2.0/3.0 * 1e9 * nu_em / (2.8 * c.Bmu *pow( E , 2) / pow(me , 2) );/* 
			pow( 1 + pow(E/me * 8980.0*sqrt(1.3e-3)/(nu_em*1e9) ,2),3.0/2.0);
	double xp = pow( 1 + pow(E/me * 8980.0*sqrt(1.3e-3)/(nu_em*1e9) ,2),3.0/2.0); */
	//std::cout << pow(E/me * 8980*sqrt(1.3e-3)/(nu*1e6) ,2) << std::endl;
	//if(xp > 30)
	//std::cout << "x0 = " << xp <<" E = " << E<< std::endl;;//  << " B = " << c.bfield_model( r ) << std::endl;
	double dpsyn = psyn0 * c.Bmu * 0.5 * pow(  sin(theta) , 2) * fff( x  /sin(theta) ); 
	
	
	//std::cout << " dpsyn = "<< dpsyn << ", B = " << c.bfield_model(r) << std::endl;
	return dpsyn;

}

double gslInt_psyn(  double nu, double E, double r){			//int over theta

		
	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;

	std::vector<double> psynParams (3);

	psynParams[0] = E;
	psynParams[1] = r;
	psynParams[2] = nu;

	gsl_function F;
	F.function = &dpsyn;
	F.params = &psynParams;

	gsl_integration_qags (&F, 1e-16, pi, 0, 1e-2, 1000, //x?
	                    w, &result, &error); 

	gsl_integration_workspace_free (w);

	return result;

}


double djsyn(double E , void * params){

	std::vector<double> jsynParams = *(std::vector<double> *)params;
	double nu = jsynParams[0];
	double r = jsynParams[1];

	double djsyn = 2* gslInt_psyn(nu, E, r)* dndeeq(E , r);


	return djsyn;
}


double gslInt_jsyn(double nu, double r){ 				// int over E

			///////////
	std::clock_t start;
	double duration;
	start = std::clock();
	int a ; 
		///////////


	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;
	std::vector<double> jsynParams (2);

	jsynParams[0] = nu;
	jsynParams[1] = r;

	gsl_function F;
	F.function = &djsyn;
	F.params = &jsynParams;

	gsl_integration_qags (&F, me, p.mx, 0, 1e-2, 1000,
	                    w, &result, &error); 


	gsl_integration_workspace_free (w);

	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	//std::cout << "alpha = "<< c.alpha << ", jsyn( r = " << r/mpc2cm*1000 <<" ) duration: " << duration <<std::endl;  //      ~30s-60s
	//std::cout << " . ";


	return result;

}


double dssyn( double r, void * params ){

	double nu = *(double *)params;


	double dist_z = Dist() / (1+c.z);

	double ssynIntegrand = 4 *pi /pow(dist_z , 2) *pow(r,2)  *  pow(c.DM_profile(r) , 2)*gslInt_jsyn(nu, r);	
	
	return ssynIntegrand;
}


double gslInt_ssyn( double nu, double r ){				// int over r

	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;



	gsl_function F;
	F.function = &dssyn;
	F.params = &nu;

	gsl_integration_qags (&F, 1e-16, r, 0, 1e-2, 1000,
	                w, &result, &error); 
	 result *= GeVJy * p.sv/(8* pi*pow( p.mx , 2.0 ));
	return result;

}

void runFlux(double mx, double r){ //runs through values of root_dv
	
	p.mx = mx;
	int n = 50;


	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();
	int a ; 
	///////before algorithm

	std::ostringstream makefilename;
	makefilename << "NSD_3_Flux_" << p.mx << "Gev_coma.txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());

	for (int i = 0 ; i < n + 1; ++i  ){


		double nu_min = 10;
		double nu_max = 1e5;
	
		double nu = nu_min * ( exp(    (log(nu_max) - log(nu_min))/ n * i));
		double data = gslInt_ssyn( nu , r);
		file << nu << "\t" <<  data <<std::endl;
		std::cout << "\n" << nu << "\t" << data << std::endl;

	};
	
	////////after algorithm
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "Total time:  " << duration <<std::endl;
}
//NOTE that we need to set a cross section for this, use eq 1 in STorm for Snu.

main(){
	dsinit_(); //initialixe DarkSUSY
	

	double r = c.rh*mpc2cm;


	p.ch = 13;	
	p.sv = 8.8e-26;	
	c.Bmu = 8;					//darksusy channel
	//std::cout << c.bfield_model(r) << std::endl;
	runFlux( 81, r ) ;

	p.ch = 25;		
	p.sv = 4.7e-25;	
	c.Bmu = 1.2;							//darksusy channel
	//std::cout << c.bfield_model(r) << std::endl;
 	runFlux( 40, r );





}



