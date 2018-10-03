/*
 R1rho simulation program with explicit exchange between two exchanging spin systems.
 This simulation applies a 1H pi pulse in the center of the spinlock to suppress CSA/D CCR.
 1H pulse is delta pulse.
*/

/* Dominique Marion, Diego F. Gauto, Isabel Ayala, Audrey Hessel, Paul Schanda 
Univ. Grenoble Alpes, CEA CNRS, Institute for Structural Biology (IBS) 71 avenue des martyrs,
38044 Grenoble (France); E-mail: dominique.marion@ibs.fr, paul.schanda@ibs.fr
*/

#include "gamma.h"
#include <sys/time.h>
#include <sys/resource.h>

#define MAXSPINS 8

using namespace std;

complex powl(const complex& z, int m)  
{ if(norm(z) == 0)
    return complex(0);
  else
    return complex(exp(m * log(z)));
}

super_op powl(const super_op& LOp, int power) 
{ super_op LOp1;
  complex z;
  LOp.set_EBR();              // First put LOp in its eigenbasis
  LOp1=LOp;
  for(int i=0; i<LOp.LS(); i++)
  { z = LOp.get(i,i);
    LOp1.put(i,i, powl(z, power));
  }
  return LOp1;
}

int main(int argc, char *argv[])

{
  sys_dynamic axc[2];
  multi_sys ax;
  gen_op temp, zz,Up[4],Upp[4],Ham, H[2][5], sigmax, sigma1, sigma; 
  
  gen_op detect, detect_both, detectSxIa, detectSxIb, Hrf; // added Hrf
  gen_op detectNz, detectHz, detectHy, detectHx;
  super_op Iop,K,L,U,UT,Un,Un1,L0;
  space_T Adip[2][MAXSPINS][MAXSPINS], Adip_R[2][MAXSPINS][MAXSPINS];
  space_T Acsa[2][MAXSPINS], Acsa_R[2][MAXSPINS];
  double J[2][MAXSPINS][MAXSPINS],D[2][MAXSPINS][MAXSPINS];
  double iso_CSA[2][MAXSPINS],eta_CSA[2][MAXSPINS], delta_CSA[2][MAXSPINS];
  int i,j,k,m,n,Fnp,count,qu,steps,nrotor;
  int coherence, decoupling, orientation, Hstate, val1, val2;
  string name, names, Spin_Sys_Name;
  const double thetam=54.73561032;
  double mas_freq;
  double ltime, time, scale;
  double delta_omega, delta_dipolar;
  int nstart, nspins[2];
  double alpha,beta,gamma;
  double Population[2];
  double alpha_CSA[2][MAXSPINS],beta_CSA[2][MAXSPINS];
  double gamma_CSA[2][MAXSPINS];
  double alpha_D[2][MAXSPINS][MAXSPINS],beta_D[2][MAXSPINS][MAXSPINS];
  double gamma_D[2][MAXSPINS][MAXSPINS];
  double tauc, kex,k_ex,kd;
  double Ttotal;   /* Total relax delay */

  double J_coupling, D_coupling, D_alpha, D_beta, D_gamma;
  double N_shift, H_shift;
  double N_CSA_delta, N_CSA_eta, N_CSA_alpha, N_CSA_beta, N_CSA_gamma; 
  double H_CSA_delta, H_CSA_eta, H_CSA_alpha, H_CSA_beta, H_CSA_gamma;


  struct rusage me;
  struct timeval tp;
  int sampling_taur;  

  double gamB1;	


  int value1[] = {1, 50, 100, 144, 200, 300, 538, 1154};
  int value2[] = {1,  7,  27,  11,  29,  37,  55,  107};
  int value3[] = {1, 11,  41,  53,  79,  61, 229,  271};

  cout << "\nR1rho simulation under MAS with stochastic motional processes for NH        \n";
  cout << "==========================================================\n\n";
  cout << "Program version: " << __FILE__ << " compiled at " << __DATE__ ", "
       << __TIME__ << "\n\n";

  
  cout << " Environnement variables \n\n";
  J_coupling = atof(std::getenv("J_coupling"));
  D_coupling = atof(std::getenv("D_coupling"));
  D_alpha = atof(std::getenv("D_alpha"));
  D_beta = atof(std::getenv("D_beta"));
  D_gamma = atof(std::getenv("D_gamma"));
  N_shift = atof(std::getenv("N_shift"));
  H_shift = atof(std::getenv("H_shift"));
  N_CSA_delta = atof(std::getenv("N_CSA_delta"));
  N_CSA_eta = atof(std::getenv("N_CSA_eta"));
  N_CSA_alpha = atof(std::getenv("N_CSA_alpha"));
  N_CSA_beta = atof(std::getenv("N_CSA_beta"));
  N_CSA_gamma = atof(std::getenv("N_CSA_gamma"));
  H_CSA_delta = atof(std::getenv("H_CSA_delta"));
  H_CSA_eta = atof(std::getenv("H_CSA_eta"));
  H_CSA_alpha = atof(std::getenv("H_CSA_alpha"));
  H_CSA_beta = atof(std::getenv("H_CSA_beta"));
  H_CSA_gamma = atof(std::getenv("H_CSA_gamma"));
  Spin_Sys_Name = std::getenv("Spin_Sys_Name");
 
  cout << "Parameters:\n";
  cout << "Magic Angle Rotation          : " << thetam << " Degree\n";
   

  std::cout << "J_coupling     " << J_coupling << " Hz \n";
  std::cout << "D_coupling     " << D_coupling << " Hz \n";
  std::cout << "D Euler        (" << D_alpha << " , " << D_beta << " , " << D_gamma << ") \n";
  std::cout << "Iso shifts     [N] " << N_shift << " Hz  [H] " << H_shift << " Hz \n";
  std::cout << "CSA       [N]  (" << N_CSA_delta << " , " << N_CSA_eta << ") Hz \n";
  std::cout << "CSA Euler [N]  (" << N_CSA_alpha << " , " << N_CSA_beta << " , " << N_CSA_gamma << ") \n";
  std::cout << "CSA       [H]  (" << H_CSA_delta << " , " << H_CSA_eta << ") Hz \n";
  std::cout << "CSA Euler [H]  (" << H_CSA_alpha << " , " << H_CSA_beta << " , " << H_CSA_gamma << ") \n" ;

  count = 1;
  /*
  query_parameter(argc,argv,count++,"Spin System Name            ? ", names);
//setup for the spin system

  names = "2spin_2perc.sys";
  */
  ax.read(Spin_Sys_Name);
  axc[0] = ax.Comp(0);
  axc[1] = ax.Comp(1);

  nspins[0]=axc[0].spins();
  nspins[1]=axc[1].spins();
  
  cout << "size of spin system            : " << nspins[0]+nspins[1] << " spins\n";

  
  for(int conformers=0;conformers<2;++conformers)
  { for(i=0;i<nspins[k]-1;++i)
    { for(j=i+1;j<nspins[k];++j)
      { 
      
       J[conformers][i][j] = J_coupling;
       D[conformers][i][j] = D_coupling;
       alpha_D[conformers][i][j] = D_alpha;
       beta_D[conformers][i][j] = D_beta;
       gamma_D[conformers][i][j] = D_gamma;
      }
    }
  }
  for(int conformers=0;conformers<2;++conformers)
  { 
    iso_CSA [conformers][0] = H_shift;
    delta_CSA [conformers][0]= H_CSA_delta;
    eta_CSA [conformers][0] = H_CSA_eta;
    alpha_CSA [conformers][0]= H_CSA_alpha;
    beta_CSA [conformers][0]= H_CSA_beta;
    gamma_CSA [conformers][0]= H_CSA_gamma;
    
    iso_CSA [conformers][1] = N_shift;
    delta_CSA [conformers][1]= N_CSA_delta;
    eta_CSA [conformers][1] = N_CSA_eta;
    alpha_CSA [conformers][1]= N_CSA_alpha;
    beta_CSA [conformers][1]= N_CSA_beta;
    gamma_CSA [conformers][1]= N_CSA_gamma;
  
  }
  
  qu = atof(std::getenv("Powder_quality"));
  steps = atof(std::getenv("Time_steps"));
  mas_freq = atof(std::getenv("Spin_speed"));
  sampling_taur = atof(std::getenv("Sampling_rate"));
  Ttotal = atof(std::getenv("Relax_delay"));
  kd = atof(std::getenv("Random_field"));
  
  coherence = atoi(std::getenv("Coherence"));
  Hstate = atoi(std::getenv("Hstate"));
  decoupling = atoi(std::getenv("Decoupling"));
  orientation = atoi(std::getenv("Orientation"));

 
/* 
  query_parameter(argc,argv,count++,"Powder Quality (cheng)      ? ", qu);
  query_parameter(argc,argv,count++,"Number of time steps        ? ", steps);
  query_parameter(argc,argv,count++,"spinning speed              ? ", mas_freq);
  query_parameter(argc,argv,count++,"Rate of sampling, in units of tau_r ?", sampling_taur); // To avoid sampling too many points every few microsec.
  query_parameter(argc,argv,count++,"Total relaxation delay      ? ", Ttotal);
  query_parameter(argc,argv,count++,"Random field constant       ? ", kd);

  query_parameter(argc,argv,count++,"Population                  ? ", Population[0]);
  Population[0] = Population[0]/100.0;
  Population[1] = 1.0 - Population[0];
*/
  query_parameter(argc,argv,count++,"Jump rate constant          ? ", k_ex);
  query_parameter(argc,argv,count++,"RF field strength           ? ", gamB1);
  query_parameter(argc,argv,count++,"Chemical shift diff.        ? ", delta_omega);
  
  iso_CSA [1][1] = N_shift + delta_omega;
  
  query_parameter(argc,argv,count++,"Dipolar Tensor Angle        ? ", delta_dipolar);
 
  gamma_D[1][0][1] = D_gamma + delta_dipolar;
  gamma_CSA [1][1]= N_CSA_gamma + delta_dipolar;
 
  query_parameter(argc,argv,count++,"Output Filename             ? ", name);
/*
  query_parameter(argc,argv,count++,"Det. coherence (Na Nb Nx)   ? ", coherence);
  query_parameter(argc,argv,count++,"NoDec Dec MultiDec          ? ", decoupling);
  query_parameter(argc,argv,count++,"Orientation                 ? ", orientation);
*/
  time       = (1.0/mas_freq)/steps;
  kex        = 1.0/2.0*k_ex;  // There is a factor of two offset in what GAMMA thinks is kex
                                // and what I think is kex (namely kforward+kbackward)
//    cout << "Ttotal:         " << Ttotal << "\n";
  Fnp = Ttotal * mas_freq / (sampling_taur * 2);	// data points to be sampled. ; factor /2 comes from splitting the spinlock in two halves.
/*
                1/mas_freq
            |------------------|  x sampling_taur * 2
            										^
                 step                              echo
*/
    

  
/*    
  for(conformers=0;conformers<2;++conformers)
  { for(i=0;i<nspins[conformers]-1;++i)
    { for(j=i+1;j<nspins[conformers];++j)
      { cout << "J       coupling constant (" << i << "," << j << ") : " << 
                J[conformers][i][j] << " Hz\n";
        cout << "Dipolar coupling constant (" << i << "," << j << ") : " << 
                D[conformers][i][j] << " Hz\n";
        cout << "Relative orientation of D tensor: (" << alpha_D[conformers][i][j] << "," <<
                 beta_D[conformers][i][j] << "," << gamma_D[conformers][i][j] << ")\n";
      }
    }
  }
  for(conformers=0;conformers<2;++conformers)
  { for(i=0;i<nspins[conformers];++i)
    { cout << "Isotropic chem shift(" << i << ") : " << 
              iso_CSA[conformers][i] << " Hz\n";
      cout << "Delta of CSA tensor (" << i << ") : " << 
              delta_CSA[conformers][i] << " Hz\n";
      cout << "Eta of CSA tensor   (" << i << ") : " << 
              eta_CSA[conformers][i] << "\n";
      cout << "Relative orientation of CSA tensor: (" << alpha_CSA[conformers][i] << "," <<
               beta_CSA[conformers][i] << "," << gamma_CSA[conformers][i] << ")\n";
    }
  }

*/
  cout << "Powder Quality Number:         " << qu << "  (" << value1[qu] << " orientations)\n";
  cout << "Ttotal:                        " << Ttotal << "\n";
  cout << "Maximum relaxation delay:      " << Ttotal << " s\n";
  cout << "Sampling every                 " << sampling_taur << "th rotor period\n";
  cout << "Number of data points (Fnp):   " << Fnp << " points\n";
  cout << "MAS frequency:                 " << mas_freq << " Hz\n";
  cout << "Time increments:               " << time << "s\n";
  cout << "Exchange rate constant         " << 2.0*kex << " s-1\n";
  cout << "Random field constant          " << kd << " s-1\n";
  cout << "Output filename:               " << name << "\n";
  cout << "RF fields (gamma*B1, in Hz)    " << gamB1 << "\n";

  cout << "Populations                    [" << ax.pop(0) << "]  [" << ax.pop(1) << "] (from file)\n";
  cout << "15N chem diff                  " << delta_omega << " Hz\n";
  cout << "Dip. tensor angular change     " << delta_dipolar << " deg\n";


  cout << "\n";
  cout << "Nb of delays = " << Fnp <<"\n"; 
  cout << "Delays = [0:"<< 2.0*sampling_taur/mas_freq << ":" << 2.0*(Fnp-1)*sampling_taur/mas_freq <<"];\n";
  cout << "\n";

    
  for(int conformers=0;conformers<2;++conformers)
  { for(i=0;i<nspins[conformers]-1;++i)
    { for(j=i+1;j<nspins[conformers];++j)
      { cout << "J       coupling constant (" << i << "," << j << ") : " << 
                J[conformers][i][j] << " Hz\n";
        cout << "Dipolar coupling constant (" << i << "," << j << ") : " << 
                D[conformers][i][j] << " Hz\n";
        cout << "Relative orientation of D tensor: (" << alpha_D[conformers][i][j] << "," <<
                 beta_D[conformers][i][j] << "," << gamma_D[conformers][i][j] << ")\n";
      }
    }
  }
  for(int conformers=0;conformers<2;++conformers)
  { for(i=0;i<nspins[conformers];++i)
    { cout << "Isotropic chem shift(" << i << ") : " << 
              iso_CSA[conformers][i] << " Hz\n";
      cout << "Delta of CSA tensor (" << i << ") : " << 
              delta_CSA[conformers][i] << " Hz\n";
      cout << "Eta of CSA tensor   (" << i << ") : " << 
              eta_CSA[conformers][i] << "\n";
      cout << "Relative orientation of CSA tensor: (" << alpha_CSA[conformers][i] << "," <<
               beta_CSA[conformers][i] << "," << gamma_CSA[conformers][i] << ")\n";
    }
  }


  cout.flush();

  ax.Kex(kex,0);
  zz = Fx(axc[0]);
    
  Iop = commutator(multize(zz,ax,0));
    
  matrix mx(ax.LS(),ax.LS(),i_matrix_type);
  Iop.put_mx(mx);
    
  block_2D data(2,Fnp);
  for(m=0;m<Fnp;++m)
    {
      data(0,m) = 0;
      data(1,m) = float (m) * 2.0*sampling_taur/mas_freq;
    }

    
//setup for the space tensor
  matrix help(3,3,0);
  for(k=0;k<2;++k)
  { for(i=0;i<nspins[k]-1;++i)
    { for(j=i+1;j<nspins[k];++j)
      { help.put_h(-1.0/2.0,0,0);
        help.put_h(-1.0/2.0,1,1);
        help.put_h( 1.0,2,2);
        help   = - (complex) D[k][i][j] * help;
        Adip[k][i][j] = A2(help);
        Adip[k][i][j] = Adip[k][i][j].rotate(alpha_D[k][i][j],beta_D[k][i][j],gamma_D[k][i][j]);
      }
    }
  }
  for(k=0;k<2;++k)
  { for(i=0;i<nspins[k];++i)
    { help.put_h(-1.0/2.0*(1.0+eta_CSA[k][i]),0,0);
      help.put_h(-1.0/2.0*(1.0-eta_CSA[k][i]),1,1);
      help.put_h( 1.0,2,2);
      help = (complex) delta_CSA[k][i] * help;
      Acsa[k][i] = A2(help);
      Acsa[k][i] = Acsa[k][i].rotate(alpha_CSA[k][i],beta_CSA[k][i],gamma_CSA[k][i]);
    }
  }

  string name1 = name+".mat";
  string name2 = name;


	/* 1H pi pulse */
	zz=Ixypuls_U(axc[0],"1H",0.0,180.0);
	Upp[0] = multize(zz,ax,0);
	zz=Ixypuls_U(axc[1],"1H",0.0,180.0);
	Upp[0] += multize(zz,ax,1);
	
    zz = gamB1 * Fx(axc[0],"15N");
    Hrf = multize(zz,ax,0);
    zz = gamB1 * Fx(axc[1],"15N");
    Hrf += multize(zz,ax,1);
   
// Random field contribution
   
    zz = Fz(axc[0],"1H");
    temp = multize(zz,ax,0);
    zz = Fz(axc[1],"1H");
    temp += multize(zz,ax,1);
    L0  = kd*d_commutator(temp);
    
    zz = Fx(axc[0],"1H");
    temp = multize(zz,ax,0);
    zz = Fx(axc[1],"1H");
    temp += multize(zz,ax,1);
    L0  += kd*d_commutator(temp);
    
    zz = Fy(axc[0],"1H");
    temp = multize(zz,ax,0);
    zz = Fy(axc[1],"1H");
    temp += multize(zz,ax,1);
    L0  += kd*d_commutator(temp);



//here starts the powder loop
//reference JCP 59 (8), 3992 (1973).

if (orientation == 0)
	{ val1 = 1; 
	  val2 = value1[qu]; 
	}
	else
	{ if (orientation > value1[qu])
	  {  orientation = value1[qu];
	     cout << "Orientation set to " << value1[qu]<< "\n";
	  }
	  val1 = orientation;
	  val2 = orientation;
	}
	

//  for(count=1; count<=value1[qu]; ++count)
   for(count=val1; count<=val2; ++count)
  { beta  = 180.0 * count/value1[qu];
    alpha = 360.0 * ((value2[qu]*count) % value1[qu])/value1[qu];
    gamma = 360.0 * ((value3[qu]*count) % value1[qu])/value1[qu];
    if ((count % 10 == 1) || (orientation >0))
    { getrusage(0, & me);
      cout << count << "\tbeta = " << beta << "\talpha = "
           << alpha << "\tgamma = " << gamma
           << ",\ttime used: " << me.ru_utime.tv_sec << " seconds\n";
      cout.flush();
    }
    scale = sin(beta/180.0*PI);
    sigmax = Fx(axc[0],"15N");
    sigma1 = multize(sigmax,ax,0)*ax.pop(0); // in welcher Basis ist der Dichteoperator?
    sigmax = Fx(axc[0],"15N"); //hier könnte 1 stehen
    sigma1 += multize(sigmax,ax,1)*ax.pop(1);
    
// Adding 1H initial state Mz 
    
    if (Hstate == 1)
    {
    sigmax = Fz(axc[0], "1H");
    sigma1 += multize(sigmax,ax,0)*ax.pop(0);
    sigmax = Fz(axc[0], "1H");
    sigma1 += multize(sigmax,ax,1)*ax.pop(1);
	}
    if (Hstate == -1)
    {
    sigmax = -Fz(axc[0], "1H");
    sigma1 += multize(sigmax,ax,0)*ax.pop(0);
    sigmax = -Fz(axc[0], "1H");
    sigma1 += multize(sigmax,ax,1)*ax.pop(1);
	}

 // Detection: Set to either detect both states, or only one state. Comment out what is not needed:
    // Case 1: detect both states
    detect_both = multize(Fm, "15N", ax);

      
    // Case 2: detect only one state
    zz     = Fm(axc[0],"15N");
    detect = multize(zz, ax, 0);

    zz     = Ix(axc[0],1)*Ia(axc[0],0);
    detectSxIa = multize(zz, ax, 0);
    zz     = Ix(axc[0],1)*Ib(axc[0],0);
    detectSxIb = multize(zz, ax ,0);

    zz     = Fz(axc[0],"15N");
    detectNz = multize(zz, ax, 0);

    zz     = Fz(axc[0],"1H");
    detectHz = multize(zz, ax, 0);

    zz     = Fy(axc[0],"1H");
    detectHy = multize(zz, ax, 0);

    zz     = Fx(axc[0],"1H");
    detectHx = multize(zz, ax, 0);
    

      
//////////////////////////////////
      
//now we rotate the space tensor
    for(k=0;k<2;++k)
    { for(i=0;i<nspins[k]-1;++i)
      { for(j=i+1;j<nspins[k];++j)
        { Adip_R[k][i][j] = Adip[k][i][j].rotate(alpha,beta,gamma);
        }
      }
    }
    for(k=0;k<2;++k)
    { for(i=0;i<nspins[k];++i)
      { Acsa_R[k][i] = Acsa[k][i].rotate(alpha,beta,gamma);
      }
    }

//zero all components
    for(k=0;k<2;++k)
      for(i=0;i<5;++i)
        H[k][i] = gen_op();

    for(k=0;k<2;++k)
    {  for(j=-2;j<=2;++j)
      { for(m=0;m<nspins[k]-1;++m)
        { for(n=m+1;n<nspins[k];++n)
          { if(axc[k].isotope(m) != axc[k].isotope(n))
            { H[k][j+2] += Adip_R[k][m][n].component(2,j) * d2(j,0,thetam) * 1/sqrt(6.0)*2*Iz(axc[k],m)*Iz(axc[k],n);
              if(j==0)
              { H[k][j+2] += J[k][m][n]*Iz(axc[k],m)*Iz(axc[k],n);
              }
            }
            else
            { H[k][j+2] += Adip_R[k][m][n].component(2,j) * d2(j,0,thetam) * 1/sqrt(6.0)* 
                        (2*Iz(axc[k],m)*Iz(axc[k],n)-(Ix(axc[k],m)*Ix(axc[k],n)+Iy(axc[k],m)*Iy(axc[k],n)));
              if(j==0)
              { H[k][j+2] += J[k][m][n]*(Iz(axc[k],m)*Iz(axc[k],n)+Ix(axc[k],m)*Ix(axc[k],n)+Iy(axc[k],m)*Iy(axc[k],n));
              }
            }
          }
        }
        for(m=0;m<nspins[k];++m)
        { H[k][j+2] += Acsa_R[k][m].component(2,j) * d2(j,0,thetam) * 1/sqrt(6.0)*2*Iz(axc[k],m);
          if(j==0)
          { H[k][j+2] += iso_CSA[k][m]*Iz(axc[k],m);
          }
        }
      }
    }

    U = Iop;  // Identity Operator
//now we calculate the propagator over one cycle of the MAS
    for(ltime=0.5*time;ltime<=steps*time;ltime += time)
    { //Ham=gen_op();
        Ham = Hrf;
      for(k=0;k<2;++k)
      { for(i=-2;i<=2;++i)
          
          Ham += exp(complex(0,i*2.0*PI*ltime*mas_freq)) * multize(H[k][i+2],ax,k);
      }
      L = complex(0,2.0*PI)*commutator(Ham);
      if ( kd != 0.0 ) 
      { L = L + L0;
      }
        
      K = Xnm(ax); // macht die Austauschmatrix.
      L += K; // hier haben wir vollen 32x32 Liouville Operator
      UT = exp(L,-time);  
      UT.set_HBR();
      U = UT*U;  
    }  
	  
    U.set_HBR();
    sigma1.set_DBR(); 
    sigma=sigma1;
      
      Un=powl(U,sampling_taur);   // make a longer propagator, in order to sample less often than every taur
	  Un.set_HBR();	    // set to Hilbert space
     Un1=Un;   // set the propagator
	  
	// Spin-lock with 1H pi pulse in the center:
	  
    for(i=0;i<Fnp;++i) // incremented by one, but the minimum increment is in fact 2*sampling_taur, because inside we loop twice over k. 
    { 
      if (coherence == -1)
      data(0,i) += proj(sigma,detect_both)*scale;
      if (coherence == 0)
      data(0,i) += proj(sigma,detectSxIa)*scale;
      if (coherence == 1)
      data(0,i) += proj(sigma,detectSxIb)*scale;
      if (coherence == 2)
      data(0,i) += proj(sigma,detect)*scale;
      if (coherence == 3)
      data(0,i) += proj(sigma,detectNz)*scale;
      if (coherence == 4)
      data(0,i) += proj(sigma,detectHz)*scale;
      if (coherence == 5)
      data(0,i) += proj(sigma,detectHy)*scale;
      if (coherence == 6)
      data(0,i) += proj(sigma,detectHx)*scale;
      
       
      
      if (decoupling == 0)  
      // No Dec
      {
         sigma=Un*sigma;
         sigma=Un*sigma;
      }
      if (decoupling == 1)  
      // Dec  
      {
         sigma=sigma1;	// set sigma back to initial value, before applying the whole spinlock-pulse-spinlock sequence.
         sigma=Un1*sigma;
	     sigma.sim_trans_ip(Upp[0]);	// 1H pi pulse
	     sigma=Un1*sigma;
	     Un1=Un*Un1;
      }
      if (decoupling == 2)  
	  // Multi Dec		
	  {
	     sigma=Un1*sigma;
         sigma.sim_trans_ip(Upp[0]);    // 1H pi pulse
         sigma=Un1*sigma;
      }    
		
    }	  

  } // end of powder loop

  
  
  MATLAB(name1,name2,data,1);
  cout << "1st point  " << data(0,0) << data(0,1) << data(0,2) << "\n";
  cout << "Delays     " << data(1,0)*1e6 << data(1,1)*1e6 << data(1,2)*1e6 << "\n";
  cout.flush();
  exit(0);
}
