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
