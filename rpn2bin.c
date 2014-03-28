#include <string.h>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* malloc, free */
//#include <outputParams.h>
/*
int initStep, finaStep, Step_rsti, outrep
*/
#include <rpnmacros.h>
extern convip();  /* rmnlib.a */

int nlev = 102 ;
int initStep = 0, finaStep = 6000, Step_rsti = 2000;
int tmstep = 178 ;
char outrep[] = "/home/deji/" ;

#include <ctype.h>
char *trim(char *str)
{
  if(str) { /* Don't forget to check for NULL! */
    /* First remove trailing spaces */
    size_t n;
    n = strlen(str);
    while (n > 0 && isspace((unsigned char)str[n - 1])) {
      n--;
    }
    str[n] = '\0';

    /* Then remove leading spaces */
    n = 0;
    while (str[n] != '\0' && isspace((unsigned char)str[n])) {
      n++;
    }
    memmove(str, str + n, strlen(str) - n + 1);
  }

  return str;
}

void main() {

   int ct, ct0, ct1  ;
   char cnum0[10], cnum[10] ;
   char rpn_prefix[] = "2007030300-00-00_", basenom[80] ;
   char dir_prefix[100], dir_prefix2[100], staging[100], flnom[100], flnom3[100] ;

   int Schm_phys = 1, Schm_dyn = 1 ;
   int do_3d = 0, do_2d = 1 ;

   int iun, iun3, handle, key, ier, ni, nj, nk, ip1, ip2, ip3, ig1, ig2, ig3, ig4, dateo, datev, datyp, deet,nbits, npak, npas, swa, lng, dltf, ubc, extra1, extra2, extra3, rewrit;
   char etiket[13], nomvar[5], typvar[3], grtyp[2] ;

   int nliste, liste[nlev] ;
   int cod, kin, k ;
   long int lgc ;
   int Level_ip1[nlev] ;
   float sig2[nlev], sig[nlev] ;

   /*
   float * ptemp, * puwind, *pvwind ;
   float * ppres, * psurf ;
   */
   float ppres[91*45], psurf[91*45] ;

   ct    = 0 ;
   ct1   = initStep ;
    /* Start the main loop */
   strcpy(dir_prefix, strcat(trim(outrep), "OUT_laststep_")) ;
   for (ct0 = 2000; ct0 <= finaStep; ct0 += Step_rsti) {
      snprintf(cnum0, 5, "%04d", ct0) ;
       /* The output sub-directory (GEMv4 style) */
      strcpy(staging, dir_prefix) ;
      strcat(strcat(dir_prefix, cnum0), "/input/000-000/") ;

       /* Step through the rpn files in the output directory */
      while (ct1 <= (ct0 + Step_rsti)) { 
                        /* Concatenate the rpn_prefix with the */
                        /* integration number to form the rpn file basename */
        snprintf(cnum, 6, "%05d", ct1) ;
        strcpy(basenom, rpn_prefix) ;
        strcat(strcat(basenom, cnum), "p") ;

        strcpy(dir_prefix2, dir_prefix) ;
         /* Physics output */
        if (Schm_phys) {
                        /* Get the full pathname and store in flnom */
           strcpy(flnom, strcat(strcat(dir_prefix, "pm"), basenom)) ;
                        /* Check to ensure the file is present, then */
                        /* open it using functions from rmnlib */
           FILE *fp = fopen(flnom, "r") ;
           if (fp) {
             ier = c_fnom(&iun, flnom, "STD+RND", 0) ;
             ier = c_fstouv(iun, "") ;
             fclose(fp) ;
           } else {
             printf ("FILE DOES NOT EXIST %s\n", flnom) ;
             break ;
           }
             printf ("was here in schm_phy %s %d\n", flnom, iun) ;
        }
        strcpy(dir_prefix, dir_prefix2) ;

         /* Dynamical core output */
        if (Schm_dyn) {
                        /* Get the full pathname and store in flnom */
           strcpy(flnom3, strcat(strcat(trim(dir_prefix2), "dm"), basenom)) ;
                        /* Check to ensure the file is present, then */
                        /* open it using functions from rmnlib */
           FILE *fp = fopen(flnom3, "r") ;
           if (fp) {
             ier = c_fnom(&iun3, flnom3, "STD+RND", 0) ;
             ier = c_fstouv(iun3, "") ;
             fclose(fp) ;
           } else {
             printf ("FILE DOES NOT EXIST %s\n", flnom3) ;
             break ;
           }
             printf ("was here in schm_dyn %s %d\n", flnom3, iun3) ;
        }

        printf("Stdout: %d %d %s\n", ct0, ct1, flnom) ;
        if (ct == 0) {

          /* Get the record information parameters using rmnlib functions*/
          /* First, get the levels for extracting 3-D variables */
          if (do_3d) {
             ier = c_fstinl(iun3, &ni, &nj, &nk, -1, " ", -1, -1, -1, "P", 
                                                 "TT", liste, &nliste, nlev) ;
             cod = -1 ;
             lgc = 0 ;
             for (k = nlev-1; k >= 0; k--) {
               ier = c_fstprm(liste[k], &dateo, &deet, &npas, &ni, &nj, &nk, 
                                &nbits, &datyp, &ip1, &ip2, &ip3, typvar,
                                nomvar, etiket,grtyp, &ig1, &ig2, &ig3, &ig4,
                          &swa, &lng, &dltf, &ubc, &extra1, &extra2, &extra3) ;
               Level_ip1[k] = ip1 ;
                      /* The f77name macro is defined in rpnmacros.h */
               //printf(" ip1s = %d, %d %s %d\n",ip1, cod, etiket, lgc) ;
               f77name(convip) (&ip1, &sig[k], &kin, &cod, &etiket, &lgc) ;
               sig2[k] = 1000.0 * sig[k] ;
               //printf(" after ip1s = %d, sig=%f, %d %s %d\n",ip1, sig[k], cod, etiket, kin) ;
               //printf(" Sigma level %d %f %F \n", k, sig2[k], sig[k]) ;
             }
          } 
           /* Secondly, get the information parameter for 2-D variables */
          else if (do_2d) {
             key = c_fstinf(iun, &ni, &nj, &nk, -1, "", -1, -1, -1, "P", "TG");
             ier = c_fstprm(key, &dateo, &deet, &npas, &ni, &nj, &nk, &nbits,
                            &datyp, &ip1, &ip2, &ip3, typvar, nomvar, etiket,
                            grtyp, &ig1, &ig2, &ig3, &ig4, &swa, &lng, &dltf,
                                           &ubc, &extra1, &extra2, &extra3) ;
              printf("Stdout surfsace key: %d %d\n", ni, nj) ;
          }
        } /* Finish record information parameters @ ct==0 */

          /* Allocate spaces for needed variables*/
        /*
        psurf =  (float *) malloc(ni * nj) ;
        ppres =  (float *) malloc(ni * nj) ;
        ptemp  = (float *) malloc(ni * nj) ;
        puwind = (float *) malloc(ni * nj) ;
        pvwind = (float *) malloc(ni * nj) ;

        if ((ptemp == NULL) || (puwind == NULL) || (pvwind == NULL))
            exit(1) ;
        if ((psurf == NULL) || (ppres == NULL))
            exit(1) ;

        if (do_2d) {
          ier = c_fstlir(ppres, iun3, &ni, &nj, &nk, -1, etiket, -1, -1,
                                                      -1, typvar, "P0") ;
          ier = c_fstlir(psurf, iun, &ni, &nj, &nk, -1, etiket, -1, -1,
                                                      -1, typvar, "TG") ;
        }
        printf("psurf is %f, psurf is %f", *(ppres+20), *(psurf+20)) ;

        if (do_3d) {
          for (k = nlev-1; k >= 100; k--) { 
             ip1 = Level_ip1[k] ;
             ier = c_fstlir(ptemp, iun3, &ni, &nj, &nk, -1, etiket, ip1, -1,
                                                      -1, typvar, "TT") ;
             ier = c_fstlir(puwind, iun3, &ni, &nj, &nk, -1, etiket, ip1, -1,
                                                      -1, typvar, "UU") ;
             ier = c_fstlir(pvwind, iun3, &ni, &nj, &nk, -1, etiket, ip1, -1,
                                                      -1, typvar, "VV") ;
          }
        }
        */

        if (do_2d) {
          ier = c_fstlir(ppres, iun3, &ni, &nj, &nk, -1, etiket, -1, -1,
                                                      -1, typvar, "P0") ;
          ier = c_fstlir(psurf, iun, &ni, &nj, &nk, -1, etiket, -1, -1,
                                                      -1, typvar, "TG") ;
        }
        printf("psurf is %f, psurf is %f", ppres[20], psurf[20]) ;


        /* Close the opened file(s) with fst function  */
        if (Schm_phys) {
          ier = c_fstfrm(iun) ;
          ier = c_fclos(iun) ;
        }
        if (Schm_dyn) {
          ier = c_fstfrm(iun3) ;
          ier = c_fclos(iun3) ;
        }

        ct++ ; /* Output counter */
        ct1 += tmstep ;    /* Next timestep in output */

      } /* End ct1 loop */
      strcpy(dir_prefix, staging) ;

      /*
      free (psurf) ;
      free (ppres) ;
      free (ptemp) ;
      free (puwind) ;
      free (pvwind) ;
      */

   } /* End the main loop */

}
