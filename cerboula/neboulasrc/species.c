# include <math.h>
# define min(x,y)  (((x) <= (y)) ? (x) : (y))
# define max(x,y)  (((x) >= (y)) ? (x) : (y))

extern void ableit(double *,double *,double,int);
extern double *vector(int);
extern double **matrix(int,int);
extern double ***trimatrix(int,int,int);
extern void free_vector(double *);
extern void free_matrix(double **,int);
extern void free_trimatrix(double ***,int,int);
void block_tridiag(double **,double **,double **,double **,double **,double ***,double *,
                   double *,double *,double **,double **,double **,int,int);


/********************************************************     
     subroutine species
     Integrating concentration equations
     Returns species concentrations
*********************************************************/

void species(double **c,double **cs,double **css,double **cp1,double *u,double *V,
             double **J,double **Js,double ***Dij,double *l0,double *l0s,double **mapr,
             double ***dmapr,double *M,double *Ms,double *rho,double *Mw,double **result,
             double xit2,double duedx,int jstart,double lam0,double *Mw_scaled,
             double **csbis, double *Jwall,int iwr/*NANNI*/)
{

   extern int neq,nspec;
   extern int lew,stef,ram; /* Boundary conditions */
   extern int chflag;
   extern double deta,udelta,rdelta,dcsidelta,rhodelta,mudelta,K_bl;
   extern int axisymm,two_d,cone,flat_plate; /* Geometric informations */
   void wall_ci(double *);
   void wall_ci2(double *);
   void wall_stef(double *,double *,double *,double **,double,double *,
                  double,double,double *,double **,double,double *);
   void wall_ram(double *,double *,double *,double **,double,double *,
                 double *,double,double,double *,double **,double);

   int i,j,k,nheavy;
   double **aa,**bb,**cc,**dd,**F,***dF,**res;
   double *f1,*g1,*h1,**dh1dc,*rhoiw;
   double *an,*bn,*cn,*dn,**Jb,**Jbs,*pt;
   double fact,fact1,fact2,fact3,deta22,K_bl2;
   const double Pr=0.72,Le=1.2;

   aa=matrix(neq,nspec),bb=matrix(neq,nspec); cc=matrix(neq,nspec);
   dd=matrix(neq,nspec); F=matrix(neq,nspec); dF=trimatrix(neq,nspec,nspec);
   res=matrix(neq,nspec);
   f1=vector(nspec); g1=vector(nspec); h1=vector(nspec);
   an=vector(nspec); bn=vector(nspec); cn=vector(nspec); dn=vector(nspec);
   Jb=matrix(nspec,neq); Jbs=matrix(nspec,neq); pt=vector(nspec);
   dh1dc=matrix(nspec,nspec); rhoiw=vector(nspec);

   deta22 = deta*deta;
   K_bl2=K_bl*K_bl;


   for(i=0;i<neq;i++)
      for(j=0;j<nspec;j++)
         Jb[j][i] = l0[i]*Le/Pr*csbis[j][i];

   for(i=0;i<nspec;i++) ableit(Jb[i],Jbs[i],deta,neq);

   if(jstart != 0)
   {
      fact = xit2/(udelta*dcsidelta);
      fact1 = sqrt(xit2)*rdelta/dcsidelta;
      fact2 = sqrt(xit2)/(udelta*rdelta);
   }
   else if(jstart == 0)
   {
      if(axisymm)
         if(cone)
         {
            fact = 0.0;
            fact1 = 1.0; /* It is included in the definition of Ji */
            fact2 = 0.0;
         }
         else
         {
            fact = 1.0/(2.0*duedx);
            fact1 = 1.0/sqrt(2.0*rhodelta*mudelta*duedx);
            fact2 = 1.0/sqrt(2.0*duedx/rhodelta/mudelta);
         }

      if(two_d)
         if(flat_plate)
         {
            fact = 0.0;
            fact1 = 1.0; /* It is included in the definition of Ji */
            fact2 = 0.0;
         }
         else
         {
            fact = 1.0/(duedx);
            fact1 = 1.0/sqrt(rhodelta*mudelta*duedx);
            fact2 = 1.0/sqrt(duedx/rhodelta/mudelta);
         }
   }

   for(i=0;i<nspec;i++)
      for(j=0;j<nspec;j++) dh1dc[i][j] = 0.0;
   
   for(i=0;i<neq;i++)
   {
      for(j=0;j<nspec;j++)
      {
/*    Here the equation is written in the form
      an*dw/dcsi+bn*dw/deta=cn*d^2w/deta^2+dn+F   */
         an[j] = (xit2*lam0)*u[i];
         bn[j] = V[i]-K_bl2*Le/Pr*l0s[i];
         cn[j] = K_bl2*l0[i]*Le/Pr;
         dn[j] = -fact1*K_bl*Js[j][i]-K_bl2*Jbs[j][i];
      }
      for(j=0;j<nspec;j++)
      {
/*    Here the equation is written in the form
      ai*d^2w/deta^2+bi*dw/deta+ci*w=d   */
         aa[i][j] = -cn[j]/deta22; 
         bb[i][j] = bn[j]/deta;
         cc[i][j] = an[j];
         dd[i][j] = dn[j]-u[i]*cp1[j][i];
         F[i][j] = mapr[i][j]*fact/rho[i];
         for(k=0;k<nspec;k++)
         {
            dF[i][j][k] = dmapr[i][j][k]*fact/rho[i];
         }
      }
   }

/******************  Boundary conditions begin ****************************/

  /* Paolo. This is the new modified boundary condition, that replaces the
     old stef boundary condition.
     USE THIS ONE AT THE PLACE OF ram WHEN DEALING WITH CHEMICAL
     NONEQUILIBRIUM */
   if(stef)
   {
      double *cw,*Jbw,*Jiw;
      cw=vector(nspec);
      Jbw=vector(nspec);
      Jiw=vector(nspec);
      for(i=0;i<nspec;i++)
      {
         cw[i]=c[i][0];
         Jbw[i]=Jb[i][0]*K_bl;
         Jiw[i]=J[i][0];
      }
     /* wall_stef(f1,g1,h1,dh1dc,rho[0],cw,M[0],Ms[0],Mw,Dij[0],(fact2/rho[0]),Jwall); */

      /* Paolo. New boundary condition */
      wall_bc(f1,g1,h1,dh1dc,rho[0],cw,M[0],Ms[0],Mw,Jiw,Jbw,l0[0],fact2);
      free_vector(cw);
      free_vector(Jbw);
      free_vector(Jiw);

   }
   else if(lew)
   {
      for(i=0;i<nspec;i++)
      {
         f1[i]=0.0;
         g1[i]=1.0;
      }
      wall_ci(h1);
   }
   else if(ram)
   {
      double *cw,*csw;
      cw=vector(nspec); csw=vector(nspec);
      for(i=0;i<nspec;i++)
      {
         cw[i]=c[i][0];
         csw[i]=cs[i][0];
      }
      wall_ram(f1,g1,h1,dh1dc,rho[0],cw,
               csw,M[0],Ms[0],Mw,Dij[0],(fact2/rho[0]));
      free_vector(cw); free_vector(csw);
   }
  
/******************  Boundary conditions end ****************************/

   for(i=0;i<neq;i++)
   {
      for(j=0;j<nspec;j++)
      {
         aa[i][j]=aa[i][j]/Mw_scaled[j];
         bb[i][j]=bb[i][j]/Mw_scaled[j];
         cc[i][j]=cc[i][j]/Mw_scaled[j];
         dd[i][j]=dd[i][j]/Mw_scaled[j];
         F[i][j]=F[i][j]/Mw_scaled[j];
         for(k=0;k<nspec;k++)
            dF[i][j][k]=dF[i][j][k]/Mw_scaled[j];
      }
   }

   for(i=0;i<nspec;i++)
   {
      f1[i]=f1[i]/Mw_scaled[i];
      g1[i]=g1[i]/Mw_scaled[i];
      h1[i]=h1[i]/Mw_scaled[i];
      for(j=0;j<nspec;j++)
         dh1dc[i][j]=dh1dc[i][j]/Mw_scaled[i];
   }

   /* Paolo */
   nheavy=nspec;
   if(chflag) nheavy=nspec-1;
   block_tridiag(aa,bb,cc,dd,F,dF,f1,g1,h1,dh1dc,c,res,neq,nheavy);

/* Paolo */
/*   for(i=0;i<neq-1;i++)
   {
      for(j=0;j<nheavy;j++)
      {
         if(res[i][j]<0.0)
         {
	    res[i][j]=1.e-20;
         }
         result[j][i]=res[i][j];
	}
   } */

   for(i=0;i<neq-1;i++)
   {
      for(j=0;j<nheavy;j++)  /*Paolo*/
      {
         result[j][i]=res[i][j];
      }
   }

   free_matrix(aa,neq); free_matrix(bb,neq); free_matrix(cc,neq);
   free_matrix(dd,neq); free_matrix(F,neq); free_trimatrix(dF,neq,nspec);
   free_matrix(res,neq);
   free_vector(f1); free_vector(g1); free_vector(h1);
   free_vector(an); free_vector(bn); free_vector(cn); free_vector(dn);
   free_matrix(Jb,nspec); free_matrix(Jbs,nspec); free_vector(pt);
   free_matrix(dh1dc,nspec); free_vector(rhoiw);

   /*printf("\n species is running");*/

}



void neutrality(double **c)
{
   int i,j;
   double tmp;

   extern int neq,nspec;
   extern double *charge;


   for(i=0;i<neq;i++)
   {
      tmp = 0.0;
      for(j=0;j<nspec-1;j++) tmp = tmp+charge[j]*c[j][i];
      c[nspec-1][i] = -tmp/charge[nspec-1];
   }

}
