! //////////////////////////////////////////////////////////////////////
! //                                                                  //
! //             P E G A S E   4.  TERMO   L I B R A R Y              //
! //                                                                  //
! //  SUBROUTINE eq_constants                                         //
! //                                                                  //
! //  Input:  mode        type of computation p,T or rho,T            //
! //          T           temperature                       K         //
! //          p           pressure                          Pa        //
! //          rho         density                           kg/m3     //
! //          flg_anha    anharmonicity corrections flag needed to    //
! //                      compute the anharmonicity corrections       //
! //          flg_stop    if 1, stops on subroutine errors            //
! //          flg_termo   needed to know the mode of computation of   //
! //                      the equilibrium constants                   //
! //                                                                  //
! //  Output: k(maxr,4)   array of equilibrium constants:             //
! //                      k(..,1) = Kp                                //
! //                      k(..,2) = Kc                                //
! //                      k(..,3) = Kx                                //
! //                      k(..,4) = Ky                                //
! //                                                                  //
! //  This subroutine is called to compute the natural logarithms of  //
! //  all equilibrium constants at the defined conditions p,T,rho.    //
! //  Notes: nusum is the sum of the stoechiometric coefficients of   //
! //         the equation that is being processed                     //
! //         sum is used to build the sum of products NU*muo          //
! //         muo is the chemical potential divided by RUNIV T         //
! //                                                                  //
! //  Benoit Bottin, 09/10/96. Modified and cleaned F90, 28/7/97.     //
! //                                                                  //
! //////////////////////////////////////////////////////////////////////
      SUBROUTINE EQ_CONSTANTS (T,TNEQ,p,rho,k,&
     &flg_anha,flg_neq,flg_stop,flg_termo)
!
      USE global_thermo
      IMPLICIT NONE
      INTEGER :: r,s,i,nusum,flg_anha,flg_neq,flg_stop,flg_termo
      REAL(kind=8) :: T,TNEQ(4),p,rho,k(:,:),muo(nsp),sum,msum
      REAL(kind=8) :: po,pot(8)
!      
      po=1.0d0
!     
      DO s=1,nsp
          CALL SPECIES_GIBBS(s,po,T,TNEQ,flg_anha,1,flg_neq,flg_stop,&
     &    flg_termo,pot)
          muo(s)=(pot(1)+pot(8))/RUNIV/T
      END DO
      DO r=1,nr
          sum=0
          nusum=0
          DO s=1,nsp
              nusum=nusum+NU(r,s)
              sum=sum+NU(r,s)*muo(s)
          END DO
          k(r,1)=-sum
          k(r,2)=-sum-nusum*LOG(RUNIV*T)
          k(r,3)=-sum-nusum*LOG(p)
          msum=0
          DO i=1,nsp
              msum=msum+NU(r,i)*LOG(MMOL(i))
          END DO
          k(r,4)=-sum-nusum*LOG(rho*RUNIV*T)+msum
      END DO
!
      RETURN
      END
! //////////////////////////////////////////////////////////////////////
