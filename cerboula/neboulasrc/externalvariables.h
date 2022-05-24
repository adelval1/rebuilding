//newtonloop
extern double macheps;
extern int newton_it;

//geompar
extern int neq,nspec,ngrid,nwrite,nreact,nrwall;
extern double deta,dcsi;
extern double sutres;
extern int lew,ram,stef;
extern int sutit,diff_type;  /* diff_type selects Stefan Maxwell or Ramshaw */
extern int axisymm,two_d,cone,flat_plate; /* geometry of the body */
extern int only_stag,only_stag_mod,finite_thick; /* modifications to follow 
                                                      Kolesnikov in the
                                                      stagnation point  */
extern int choice, ga_choice; /*wall catalycity modells introduced by pietro
				rini*/
extern int total_therm_cond; /*if it is 1 the diffusion fluxes are taken into
                               account trough the total thermal conductivity
			       only (oper should be 1 in this case) - NANNI*/

//edgepar
extern double udelta,*cdelta,hdelta,pdelta,rhodelta,mudelta,duedelta,dhdelta;
extern double dcsidelta,rdelta,beta,tdelta,dtedelta,twalldelta;
extern double hwall,*hspwall,cpwall,kwall,muwall,*gamdelta;
extern double v_e,duedy,a_e;
extern double delta_bl,chi,K_bl;
extern int n_elements,**lambda_matrx, eq_flg;

//const_def
extern double R,kb,pi,Na;

//charge
extern double *charge;
extern int chflag;

//catalysis
/*definition of variables for the wall reactions*/
extern double **nu,***mu; /*stoichiometric matrices*/
/* Some paremeters to choose the model for the wall b.c. */
extern int lew;  /* Local equilibrium wall */
extern int stef; /* b.c. with Stefan-Maxwell equations */
extern int ram;  /* b.c. with Ramshaw approximation */

