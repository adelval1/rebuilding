/* Variables from edgepar.h */
double udelta,*cdelta,hdelta,pdelta,rhodelta,mudelta,duedelta,dhdelta,twalldelta,*gamdelta;
double dcsidelta,rdelta,beta,tdelta,dtedelta;
double hwall,*hspwall,cpwall,kwall,muwall;
double v_e,duedy,a_e;
double delta_bl,chi,K_bl;
int n_elements,**lambda_matrx, eq_flg;

/* Variables from geompar.h */
int neq,nspec,ngrid,nwrite,nreact,nrwall;
double deta,dcsi;
double sutres;
int sutit,diff_type;
int axisymm,two_d,cone,flat_plate,only_stag,only_stag_mod,finite_thick,choice,ga_choice;

/* Variables from const_def.h*/
double R,kb,pi,Na;

/* Variables from peg_def.h  */
int flg_anha,mode,flg_neq,flg_stop,flg_termo,flg_traco,oper,total_therm_cond/*NANNI*/;
int oper_traco,sonine; 

/* Variables from newtloop.h */
double macheps;
int newton_it;

/* Variables from catalysis.h */
double **nu,***mu;
int lew,stef,ram,lfw;
/* Variables from charge.h */
double *charge;
int chflag;
/*Added by NANNI*/
int n_elements;
