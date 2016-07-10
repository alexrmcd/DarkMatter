// 4-19-16
//
// Diffusion Equation with spatial diffusion, diffusion const independent of spatial coordinate -> D(E). Uses GSL and Darksusy
//
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

g++ -o FluxDensity_Synch FluxDensity_Synch.cpp  -I/home/alex/research/darksusy-5.1.2/include -L/home/alex/research/darksusy-5.1.2/lib -lgsl -lgslcblas -ldarksusy -lFH -lHB -lgfortran
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

	//Cluster Parameters
	double z ;					//redshift
	double rh; 					//halo radius in kpc
	double rcore ; 				//core radius in kpc

	//Bfield parameters
	double B0 ; 				//microGauss
	double Bmu;
	double beta;
	double eta;

	//DM density parameters
	int DM;						//selects DM profile, 0 -> NFW, 1-> Einasto
	double rhos_NFW;			//characteristic density for NFW in GeV/cm^3
	double rs_NFW;				//scale radius for NFW in kpc
	double rs_NFWa;
	double rhos_Ein;			//characteristic density for Einasto in GeV/cm^3
	double rs_Ein;				//scale radius for Einasto in kpc
	double alpha;

	//energy loss coeff. in 1e-16 Gev/s 
	double bsynch;				
	double bIC;
	double bbrem;
	double bcoul;
	
	// Diffusion Parameters
	int SD;						//switch to turn on diffusion, 0 -> NSD, 1-> SD
	double gamma;				// D(E) ~ E^gamma
	double db;					// minimum scale of uniformity for mag field, from Colafrancesco. 
	double D0;					//in cm^2/s

	// v lookup table parameters
	int vsize;
	double vscale;
	std::vector<double> vlookup;

	// Greens lookupt table parameters 
	int n_r;
	int n_rootdv;
	int rootdv_max;
	std::vector< std::vector<double> > GLUT;

	Cluster() :	name("Coma_default"),
				z(0.0232), 						
				rh(1000.0),						
				rcore(291.0),

				B0(5.0),		
				Bmu(1.0),
				beta(0.75),
				eta(0.5),

				DM(0),
				rhos_NFW(0.039974),			
				rs_NFW(25),	
				rs_NFWa(404),			
				rhos_Ein(0.08296),			
				rs_Ein(280.0),				
				alpha(0.17),
				
				bsynch(0.0253),
				bIC(0.265),
				bbrem(1.51),
				bcoul(6.13),

				SD(1),
				gamma(0.3), 
				db(20.0),
				D0(3.1e28),				
			
				vsize(1000000),					
				vlookup(vsize),
				vscale(1),					//calculated in create_vLUT()
				
				n_r((int)(rh) + 1),
				n_rootdv(1000 + 1),
				rootdv_max(100),			//calculated in createGLUT()
				GLUT( n_r , std::vector<double>(n_rootdv) )

	{	//everything labelled or Coma, should work in user options
		std::cout << "creating cluster... " << std::endl;
	}


	double bfield_model (double r) {

		double rc = rcore * kpc2cm;

		double B_field = B0 * pow(( 1 + r*r/(rc*rc)),(-1.5*beta*eta));		// Storm et al 2013 

		return B_field;

	};	

	double bloss(double E , double r){

		double bloss = bsynch*pow(bfield_model(r), 2)*E*E 					//bsyn bfield_model(r)
						+ bIC * pow(1 + z, 4 )*E*E  					//bIC
						+ bbrem*nele*(0.36 + log(E/me/nele) )			//+ 1.51*n*(0.36 + log(E/me) )*E						//bbrem Note this is most likely incorrect, no factor of E, and should be E/me/nin log
						+ bcoul*nele*( 1 + log(E/me/nele)/75); 					//bcoul

		bloss *=1e-16;	

		return bloss;
	};


	double bloss(double E ){ 												//overload so that we can just have bloss(E), could also have same as b(E, r) but set r=0 ??  kinda sloppy

		double bloss = bsynch*pow(Bmu, 2.0)*E*E 								//bsyn bfield_model(r)
						+ bIC /** pow(1 + c.z, 4 )*/*E*E  					//bIC
						+ bbrem*nele*(0.36 + log(E/me/nele) )						//brem , in Emma's code  has + 1.51*n*(0.36 + log(E/me) )*E
						+ bcoul*nele*( 1 + log(E/me/nele)/75); 					//bcoul

		//bloss *=1e-16;	

		return bloss;
	};

	double DM_profile(double r){
		double rho;
		if (DM == 0){
			//NFW
			double rs = rs_NFW * kpc2cm;
			rho = rhos_NFW / ( r/rs * pow(1 + r/rs , 2 )) ; 
		
		}
		else if (DM == 1){
			// Einasto
			double rs = rs_Ein * kpc2cm;
			rho = (rhos_Ein) * exp(-2.0/alpha * ( pow(r/rs, alpha) - 1 ));
		}

		return rho;

	};

	double D(double E){

		double D = D0 *pow( db , 2.0/3.0 )* pow(E, gamma)/pow(Bmu, 1.0/3.0);

		return D;
	};



};

class Particle {
	public:

	int ch;
	double mx;
	double sv;

	Particle() : ch(25) , mx(1000), sv(0){}

};

Cluster c;
Particle p;


///////////////////        root_dv        ////////////////////////////////


double dv(double E , void * params){

	double dv = c.D(E)/c.bloss(E);

	return dv;
}



double v( double E ){

	gsl_integration_workspace * w 
	= gsl_integration_workspace_alloc (1000);

	double result, error;

	gsl_function F;
	F.function = &dv;


	gsl_integration_qags (&F, E, p.mx, 0, 1e-8, 1000, 
	                w, &result, &error); 

	gsl_integration_workspace_free (w);

	result *= 1e16;
	return result;
}

double root_dv(double Ep,  double vE){

	double root_dv  =   sqrt( ( vE ) - v(Ep)   ) ;

	return root_dv;
}

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


double rconst(double rcm){

	double thetaB = 25; // beam size in arcsec
	double dist_z = Dist() / (1 + c.z);


	double rb = 0.5 * dist_z *thetaB/(3600) * (pi/180); //beam size in cm


	double rconst;

	if( rb < rcm){
		rconst = rcm;
	}

	else if(rb > rcm){
		rconst = rb;
	};
	
	return rconst;

}



double dGreens(double rp, void * params ){ 

	std::vector<double> greenParam = *(std::vector<double> *)params;
	double ri = greenParam[0];
	double r = greenParam[1];
	double root_dv = greenParam[2]; // 

	double dGreens = pow(root_dv , -1) * rp/ri * (exp( - pow( (rp-ri)/(2*root_dv) , 2)) 
		- exp( - pow( ( rp + ri)/(2*root_dv) , 2)) ) * pow( c.DM_profile(rp),2)/pow( c.DM_profile(r),2);
	//if(rp/kpc2cm >= 200)
	//	std::cout << r << " "<< r/kpc2cm <<", "<< rp << " "<< rp/kpc2cm<< ", "<< pow( c.DM_profile(rp),2)<<"  "<<pow( c.DM_profile(r),2) << std::endl;
	return dGreens;

}



double gslInt_Greens(double ri , double r, double root_dv, double rh){


	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;

	std::vector<double> greenParam (3);


	greenParam[0] = ri;
	greenParam[1] = r;
	greenParam[2] = root_dv;

	gsl_function F;
	F.function = &dGreens;
	F.params = &greenParam;
	gsl_set_error_handler_off();
	gsl_integration_qags (&F, 1e-16, rh, 0, 1e-3, 1000, //x?
	                w, &result, &error); 

	gsl_integration_workspace_free (w);

	return result;

}


double Greens (double r, double root_dv) {  //called by ddsyn

	double rh = c.rh * kpc2cm ;
	int imNum = 7; //number of image charges = 2*imNum + 1
	double Gsum = 0 ;
	for (int i = - imNum; i < imNum + 1; ++i ){


		double ri;
		
		if (i == 0)
			ri = r;
		else
			ri = (pow(-1 , i)*r + 2*i*rh);

		Gsum += pow(-1, i) * gslInt_Greens(ri, r, root_dv, rh);

	}

	double Greens = pow(4*pi , -1.0/2.0)*Gsum ;

		//std::cout << "rdv " << root_dv/kpc2cm << std::endl;
	return Greens;

};




double darksusy (double Ep){

	int yieldk = 151;
	int istat;

	double ds = dshayield_(&p.mx, &Ep, &p.ch, &yieldk, &istat);

	return ds;
}

double ddiffusion(double Ep, void * params){
	std::vector<double> diffusionParams = *(std::vector<double> *)params;

	double E = diffusionParams[0];
	double vE = diffusionParams[1];
	double r = diffusionParams[2];

	double ddiffusion;


	if(c.SD == 1){
		double Ep_scaled = (int)(Ep/c.vscale) ;
		double rootdv = sqrt( std::abs(vE - c.vlookup[Ep_scaled]) ); // 0.035*mpc2cm ; //  	//	
		//std::cout << "E = " << E << " vE = "<< vE <<" Ep = " << Ep << " vEp = "<< c.vlookup[Ep_scaled]<< " rootdv = "<<rootdv << "\n";
		
		double rootdv_max = c.rootdv_max*kpc2cm;

		double r_scale = c.rh/c.n_r;
		//std::cout << c.rh <<" "<< r_scale <<"\n";
		double rootdv_scale = rootdv_max/c.n_rootdv;

		
		double r_int = (int)(( 1 + r/kpc2cm )/r_scale);
		double rootdv_int = (int)(rootdv/rootdv_scale);
	//std::cout <<" rdv in kpc = "<<rootdv*1000/mpc2cm << "\n";
	//std::cout <<" rdv scaled = "<<rootdv*1000/mpc2cm/rootdv_scale << "\n";
	//std::cout <<" rdv int = "<<rootdv_int << "\n";

		if(r_int >= c.n_r)
			r_int = c.n_r - 1;

		ddiffusion = c.GLUT[r_int][rootdv_int] * darksusy(Ep);
		if(isnan(c.GLUT[r_int][rootdv_int]) == 1)
		std::cout << "GLUT " << c.GLUT[r_int][rootdv_int] << ", r = "<<r/kpc2cm << " r_int = "<< r_int <<", rootdv = " << rootdv/kpc2cm<<std::endl;
	}
	else{
		ddiffusion = darksusy(Ep);
	}

	return ddiffusion;

}

double gslInt_diffusion( double E, double r){			// int over Ep
	double E_scaled = (int)(E/c.vscale);

	double vE = c.vlookup[E_scaled];

	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;

	std::vector<double> diffusionParams (3);

	diffusionParams[0] = E;
	diffusionParams[1] = vE;
	diffusionParams[2] = r;

	gsl_function F;
	F.function = &ddiffusion;
	F.params = &diffusionParams; 								//pass Ep to rootdv(), pass r from dndeeq as well, 
	gsl_set_error_handler_off();
	gsl_integration_qags (&F, E, p.mx, 0, 1e-2, 1000, 
	                    w, &result, &error); 

	gsl_integration_workspace_free (w);

	return result;

}

double dndeeq(double E, double r ){

	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();

	///////before algorithm

	double dndeeq = (1 / c.bloss(E,r))*gslInt_diffusion(E, r);	
	////////after algorithm
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	if ( isnan(dndeeq) == 1)
	std::cout << "dndeeq(E = "<< E  << " , r = " << r/kpc2cm << " ) = "<< dndeeq <<" --> " << duration << std::endl;
	

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
	double nu_em = ( 1 + c.z )* nu/1000; // (observing freq)*(1+z)/100 convert from Ghz to MHz 


	//std::cout << E << " , " << r/mpc2cm*1000 << std::endl;

	double x = x0 *nu_em / ( c.bfield_model( r ) * pow( E , 2) );
	//std::cout << "x = " << x  << " B = " << c.bfield_model( r ) << std::endl;
	double dpsyn = psyn0 * c.bfield_model(r) * 0.5 * pow(  sin(theta) , 2) * fff( x  /sin(theta) ); 
	
	
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

	gsl_integration_qags (&F, 1e-16, pi, 0, 1e-3, 1000, //x?
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




void create_vLUT(){
		// iteration timer start
	std::clock_t vstart;
	double vduration;
	vstart = std::clock();
	///////before algorithm

	std::cout << "creating LUT..." <<std::endl; 
	
	for (int j = 0 ; j < c.vsize ; ++j ){

		c.vscale = p.mx/c.vsize;
		
		c.vlookup[j] = v(j*c.vscale);
		
	}

	////////after algorithm
	vduration = (std::clock()  -  vstart)/(double) CLOCKS_PER_SEC;
	std::cout << "vlookup time = " << vduration <<std::endl;
}

void createGLUT(){
	// iteration timer start
	std::clock_t Gstart;
	double Greensduration;
	Gstart = std::clock();
	///////before algorithm

	std::cout << "creating GLUT..." << std::endl; 


	double rh = c.rh*kpc2cm;
	double r_scale = rh/c.n_r;

	double rdv = sqrt(v(me));
	c.rootdv_max = (int)(rdv/kpc2cm + 1);
	std::cout << c.rootdv_max<<std::endl;
	//c.n_rootdv = (c.rootdv_max*20 + 1);
	double rootdv_scale = c.rootdv_max*kpc2cm/c.n_rootdv;

	c.GLUT[0][0] = 1;
	
	for (int i = 1 ; i < c.n_r ; ++i ){
				//		std::cout << i*r_scale/kpc2cm <<std::endl;																																																																																									
		for(int j = 1; j < c.n_rootdv ; ++j){


			c.GLUT[i][j] = Greens(i*r_scale, j*rootdv_scale);

			c.GLUT[0][j] = 0;
			c.GLUT[i][0] = 1;
			c.GLUT[c.n_r - 1][j] = 0;

		}

		int a = (int)(0.1*kpc2cm/rootdv_scale);
		//std::cout << i << "/" << c.n_r - 1 << " "<< c.GLUT[i][a] << std::endl;
	};
	////////after algorithm
	Greensduration = (std::clock()  -  Gstart)/(double) CLOCKS_PER_SEC;
	std::cout << "GLUT time = " << Greensduration <<std::endl;

}

void runFlux(double mx){ 
	

	p.mx = mx;
	int n = 50;
	if(c.SD == 1) 
		create_vLUT();


	double r = c.rh*mpc2cm;
	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();
	int a ; 
	///////before algorithm

	std::ostringstream makefilename;
	makefilename <<c.name << "_SynchFlux_" << p.mx << "Gev_coma.txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());

	for (int i = 0 ; i < n + 1; ++i  ){

		double nu_min = 10;
		double nu_max = 10e8;
	
		double nu = nu_min * ( exp(    (log(nu_max) - log(nu_min))/ n * i));
		double data = gslInt_ssyn( nu , r);
		file << nu << "\t" <<  data <<std::endl;
		std::cout << "\n" <<c.name << ", "<< p.mx <<" : "<< nu << "\t" << data << std::endl;

	};
	
	////////after algorithm
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "Total time:  " << duration <<std::endl;
}




main(){
	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();

	///////before algorithm
	dsinit_(); //initialixe DarkSUSY






		c.name = "A1";
		c.rh = 1000;
		c.B0 = 2;
		c.Bmu = 1;
		c.gamma = 0.5;
		createGLUT();	

		p.ch = 25;
		p.sv = 4.7e-25;								
		runFlux( 40);
		p.ch = 13;	
		p.sv = 8.8e-26;	
		runFlux( 81 ) ;



		c.name = "A2";
		c.rh = 1000;
		c.B0 = 40;
		c.Bmu = 25;
		c.gamma = 0.5;
		createGLUT();	
		p.ch = 25;
		p.sv = 4.7e-25;								
		runFlux( 40);
		p.ch = 13;	
		p.sv = 8.8e-26;	
		runFlux( 81 ) ;



		c.name = "d1";
		c.z = 1e-6;
		c.rh = 1.4;
		c.rcore = 0.1;

		c.B0 = 1.0;
		c.Bmu = 0.5;

		c.DM = 0;
		c.rhos_NFW = 7.1	;	
		c.rs_NFW = 0.28;	

		c.bsynch = 0.02549;
		c.bIC = 0.2499;
		c.bbrem = 0;
		c.bcoul = 0;

		c.gamma = 0.7;
		c.db = 1;
		c.D0 = 	pow(c.Bmu, 1/3)*3e26;		
		createGLUT();	
		p.ch = 25;
		p.sv = 4.7e-25;								
		runFlux( 40);
		p.ch = 13;	
		p.sv = 8.8e-26;	
		runFlux( 81 ) ;

		c.name = "d2";
		c.z = 1e-6;
		c.rh = 1.4;
		c.rcore = 0.1;

		c.B0 = 0.6;
		c.Bmu = 0.3;

		c.DM = 0;
		c.rhos_NFW = 7.1	;	
		c.rs_NFW = 0.28;	

		c.bsynch = 0.02549;
		c.bIC = 0.2499;
		c.bbrem = 0;
		c.bcoul = 0;

		c.gamma = 0.7;
		c.db = 1;
		c.D0 = 	pow(c.Bmu, 1/3)*3e26;		
		createGLUT();	
		p.ch = 25;
		p.sv = 4.7e-25;								
		runFlux( 40);
		p.ch = 13;	
		p.sv = 8.8e-26;	
		runFlux( 81 ) ;

	////////after algorithm
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "Total time:  " << duration <<std::endl;
}





/*
	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();
	int a ; 
	///////before algorithm




	////////after algorithm
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "Total time:  " << duration <<std::endl;




	*/