/***************************************************************************** 
 *                            RPN2NC                                         *
 *                                                                           * 
 *Object                                                                     * 
 *   Standalone program to read standard GEM model output RPN files from a   *
 *   typical GEMv4 simulation output folder, into a netCDF file              * 
 *                                                                           *
 * IN/OUT                                                                    * 
 *   In current form, the model output directory, netCDF output file,        * 
 *   and the various variables to read from model output and write into the  * 
 *   netCDF files are hardcoded.                                             *
 *   MODEL VARIABLES (and their model identifier) WRITTEN TO NC DATA;        *
 *      Surface pressure          "P0"                                       *
 *      Surface temperature       "TG"                                       *
 *      Atmospheric temperature   "TT"                                       *
 *      Atmospheric water vapor   "HU"                                       *
 *      Zonal Wind                "UU"                                       *
 *      Meridional wind           "VV"                                       *
 *   DERIVED PRODUCT WRITTEN TO NC DATA;                                     *
 *      Atmospheric total number density                                     *
 *                                                                           * 
 * NOTES:                                                                    * 
 *   The netCDF library and headers as well as working ARMNLIB environment   *
 *   (rmnlib library and headers) are needed to compile and run this program *
 *   The rmnlib library is provided by Environment Canada's CMC              *
 *                                                                           *
 * AUTHOR                                                                    *
 *   Deji Akingunola (dakingun@gmail.com)                                    *
 *                                                                           *
 *****************************************************************************/
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* calloc, free */
#include <string.h>     /* strcpy, strcat */

#include <rpnmacros.h>
extern convip();        /* rmnlib.a */

#include <netcdf.h>
/* Define th netCDF file to create. */
#define NC_FILE_NAME "rpnout1.nc"
/* Handle netCDF errors and exiting with a non-zero status. */
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;}

#define bz 1.3807E-23
#define nlev 102        /* Number of vertical levels in the model */

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

int main() {

  /* Work paramters, determined by the model information */
  /* (normally read from output configuration file) */
  int initStep = 0, finaStep = 6000, Step_rsti = 2000;
  int tmstep = 178 ;

  /* Number of vertical levels (normally read from output configuration file) */

  /* Work paramters; loop counters */
  int ct, ct0, ct1  ;
  char cnum0[10], cnum[10] ;

  /* Model output directory */
  char outrep[] = "/home/deji/" ;
  /* Parameters defining the output filename and location */
  char rpn_prefix[] = "2007030300-00-00_", basenom[40] ;
  char dir_prefix[100], dir_prefix2[100], staging[100] ;
  char flnom[100], flnom3[100] ;

  /* Conditionals for retrieving physics and dynamic outputs */
  int Schm_phys = 1, Schm_dyn = 1 ;

  /* Conditionals for retrieving 3-D and 2-D variables */
  int do_3d = 1, do_2d = 1 ;

  /* Work paramters*/
  int ier, k ;

  /* Work paramters; file unit numbers */
  int iun = 16, iun3 = 18;

  /* RPN file record parameters */
  int key, ni, nj, nk, ip1, ip2, ip3, ig1, ig2, ig3, ig4, dateo, datyp, deet, nbits, npas, swa, lng, dltf, ubc, extra1, extra2, extra3;
   char etiket[13], nomvar[5], typvar[3], grtyp[2] ;
   int nliste, liste[nlev] ;

  /* Work paramters; sigma levels */
  int Level_ip1[nlev] ;
  float sig2[nlev], sig[nlev] ;

  char styp[8] = "STD+RND", styp2[2] ;
  int lgc, cod, kin = 0 ;

  /* Atmospheric variables to extract from the RPN files */
   int surf_dim, ver_dim ;       /* surface and 3-D spatia dimensions */
   float * ptemp, pwvap ;        /* Atmospheric temperature and water vapour */ 
   float * puwind, * pvwind ;    /* Atmospheric zonal and meridional wind */
   float * ppres, * ptsurf;      /* Surface pressure and temperature */
   float * plat, *plon ;    /* Model latitude and longitude grid points */
   float * ptnd ;                   /* Total number density */   


  int ndims = 4 ;
  /** Netcdf definitions **/
  int ncid, ni_id, nj_id, nk_id, nt_id ;
      /* Coordinate and data variable ids */
  int lat_id, lon_id, lev_id ;
  int ps_id, st_id, te_id, hu_id, uu_id, vv_id, tnd_id ; 
  int dimids[ndims];
  int retval ;                  /* Error handling. */    
 
  size_t nc_start[ndims], nc_ct[ndims];     /* The start and count arrays */

  char varnom[25] ;
  char var_unit[10] ;

  /* TODO: Check ndims matches do_3d settings */
  ct    = 0 ;
  ct1   = initStep ;
  //strcat(strcpy(dir_prefix, outrep), "OUT_laststep_") ;
  strcat(strcpy(dir_prefix, trim(outrep)), "OUT_laststep_") ;

    /**** Start the main loop (looping through the output folders) ****/
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

       strcpy(dir_prefix2, dir_prefix) ;    /* Preserve dir_prefix */

         /** Physics output **/
       if (Schm_phys) {
                        /* Get the full pathname and store in flnom */
         strcpy(flnom, strcat(strcat(dir_prefix, "pm"), basenom)) ;
                        /* Check to ensure the file is present, and then */
                        /* open it using functions from rmnlib */
         FILE *fp = fopen(flnom, "r") ;
         if (fp) {
           ier = c_fnom(iun, flnom, styp, 0) ;
           fclose(fp) ;
         } else {
             printf ("FILE DOES NOT EXIST %s\n", flnom) ;
             break ;
         }
         ier = c_fstouv(iun, "") ;
       }
       strcpy(dir_prefix, dir_prefix2) ;

         /** Dynamical core output **/
       if (Schm_dyn) {
                        /* Get the full pathname and store in flnom */
         strcpy(flnom3, strcat(strcat(dir_prefix2, "dm"), basenom)) ;
                        /* Check to ensure the file is present, and then */
                        /* open it using functions from rmnlib */
         FILE *fp = fopen(flnom3, "r") ;
         if (fp) {
           ier = c_fnom(iun3, flnom3, styp, 0) ;
           fclose(fp) ;
         } else {
             printf ("FILE DOES NOT EXIST %s\n", flnom3) ;
             break ;
         }
         ier = c_fstouv(iun3, "") ;
       }

        /**** Get the file record information parameters using rmnlib functions
              on the first set of output files read ****/
       if (ct == 0) {
          /** First, get the levels for extracting 3-D variables **/
         if (do_3d) {
           strcpy(nomvar, "TT") ;
           ier = c_fstinl(iun3, &ni, &nj, &nk, -1, " ", -1, -1, -1, "", 
                                               nomvar, liste, &nliste, nlev) ;

           cod = -1 ;
           lgc = 0 ;
           for (k = nlev-1; k >= 0; k--) {
              ier = c_fstprm(liste[k], &dateo, &deet, &npas, &ni, &nj, &nk, 
                                &nbits, &datyp, &ip1, &ip2, &ip3, typvar,
                                nomvar, etiket,grtyp, &ig1, &ig2, &ig3, &ig4,
                          &swa, &lng, &dltf, &ubc, &extra1, &extra2, &extra3) ;
              Level_ip1[k] = ip1 ;

             /* Get the sigma levels */
              f77name(convip) (&ip1, &sig[k], &kin, &cod, styp, &lgc) ;
              sig2[k] = 1000.0 * sig[k] ;
              printf("ip1 is %d, sig is %f\n", ip1, sig2[k]) ;
           }
         }
           /** Obtain the information parameter for 2-D variables **/
         else {
            strcpy(nomvar, "P0") ;
            key = c_fstinf(iun3, &ni, &nj, &nk, -1, "", -1, -1, -1, "", nomvar);
            ier = c_fstprm(key, &dateo, &deet, &npas, &ni, &nj, &nk, &nbits,
                            &datyp, &ip1, &ip2, &ip3, typvar, nomvar, etiket,
                            grtyp, &ig1, &ig2, &ig3, &ig4, &swa, &lng, &dltf,
                                           &ubc, &extra1, &extra2, &extra3) ;
         }

          /** Allocate spaces for needed variables **/
         surf_dim = ni * nj ;                 /* Surface dimension */
         ver_dim = ni * nj * nlev ;           /* Vertical dimension */
         plon = (float *) calloc(ni, sizeof(plon)) ;
         plat = (float *) calloc(nj, sizeof(plat)) ;

         ppres =  (float *) calloc(surf_dim, sizeof(ppres)) ;
         if ((ppres == NULL) ) {
            printf("Memory allocation error for surface variables \n") ;
            exit(1) ;
         }

         if (do_2d) {
           ptsurf =  (float *) calloc(surf_dim, sizeof(ptsurf)) ;
           if ((ptsurf == NULL)) {
             printf("Memory allocation error for surface variables \n") ;
             exit(1) ;
           }
         }
         if (do_3d) {
           ptemp  = (float *) calloc(ver_dim, sizeof(ptemp)) ;
           pwvap  = (float *) calloc(ver_dim, sizeof(pwvap)) ;
           puwind = (float *) calloc(ver_dim, sizeof(puwind)) ;
           pvwind = (float *) calloc(ver_dim, sizeof(pvwind)) ;

           ptnd = (float *) calloc(ver_dim, sizeof(pvwind)) ;
           if ((ptemp == NULL) || (puwind == NULL) || (pvwind == NULL) || (ptnd == NULL) || (pwvap == NULL)) {
             printf("Memory allocation error for 3-D variables \n") ;
             exit(1) ;
           }
         }

                    /* Read model latitude and longitude grid points
                       from the RPN output file; reading once is sufficient */
         strcpy(nomvar, "^^") ;  strcpy(styp2, "X") ;
         ier = c_fstlir(plat, iun, &ni, &nj, &nk, -1, etiket, -1, -1,
                                                      -1, styp2, nomvar) ;
         strcpy(nomvar, ">>") ;
         ier = c_fstlir(plon, iun, &ni, &nj, &nk, -1, etiket, -1, -1,
                                                      -1, styp2, nomvar) ;

         /** netCDF file creation and variables/ coordinate definitions  **/
         /* Create the file. */
         if ((retval = nc_create(NC_FILE_NAME, NC_CLOBBER, &ncid)))
           ERR(retval);

         /* Define the dimensions and coordinate netCDF variables. */
         strcpy(varnom, "Longitude") ;
         if ((retval = nc_def_dim(ncid, varnom, ni, &ni_id)))
           ERR(retval);
         if ((retval = nc_def_var(ncid, varnom, NC_FLOAT, 1, &ni_id, &lon_id)))
           ERR(retval);
         strcpy(varnom, "Latitude") ;
         if ((retval = nc_def_dim(ncid, varnom, nj, &nj_id)))
           ERR(retval);
         if ((retval = nc_def_var(ncid, varnom, NC_FLOAT, 1, &nj_id, &lat_id)))
           ERR(retval);
         strcpy(varnom, "Sigma Levels") ;
         if ((retval = nc_def_dim(ncid, varnom, nlev, &nk_id)))
           ERR(retval);
         if ((retval = nc_def_var(ncid, varnom, NC_FLOAT, 1, &nk_id, &lev_id)))
           ERR(retval);
         strcpy(varnom, "Time") ;
         if ((retval = nc_def_dim(ncid, varnom, NC_UNLIMITED, &nt_id)))
           ERR(retval);

         /* Define units attributes for coordinate vars. */
         strcpy(varnom, "units") ; strcpy(var_unit, "degrees") ;
         if ((retval = nc_put_att_text(ncid, lat_id, varnom, 
				 strlen(var_unit), var_unit)))
           ERR(retval);
         if ((retval = nc_put_att_text(ncid, lon_id, varnom, 
				 strlen(var_unit), var_unit)))
           ERR(retval);
 
        /* Define the netCDF variables. */
         dimids[ndims-1] = ni_id;
         dimids[ndims-2] = nj_id;
         if (do_3d) 
           dimids[1] = nk_id;
         dimids[0] = nt_id;
         strcpy(varnom, "Surface pressure") ;
         if ((retval = nc_def_var(ncid, varnom, NC_FLOAT, 3, dimids, &ps_id)))
           ERR(retval);
         if (do_2d) {
           strcpy(varnom, "Surface temperature") ;
           if ((retval = nc_def_var(ncid, varnom, NC_FLOAT, 3, dimids, &st_id)))
             ERR(retval);
         }
         if (do_3d) {
           strcpy(varnom, "Atmospheric temperature") ;
           if ((retval = nc_def_var(ncid, varnom, NC_FLOAT, 4, dimids, &te_id)))
             ERR(retval);
           strcpy(varnom, "Atmospheric Water vapor") ;
           if ((retval = nc_def_var(ncid, varnom, NC_FLOAT, 4, dimids, &hu_id)))
             ERR(retval);
           strcpy(varnom, "Zonal wind") ;
           if ((retval = nc_def_var(ncid, varnom, NC_FLOAT, 4, dimids, &uu_id)))
             ERR(retval);
           strcpy(varnom, "Meridional wind") ;
           if ((retval = nc_def_var(ncid, varnom, NC_FLOAT, 4, dimids, &vv_id)))
             ERR(retval);
           strcpy(varnom, "Total Number Density") ;
           if ((retval = nc_def_var(ncid, varnom, NC_FLOAT, 4, dimids, &tnd_id)))
             ERR(retval);
         }

         /* Define units attributes for vars. */
         strcpy(varnom, "units") ; strcpy(var_unit, "mb") ; /* Surface pres.*/
         if ((retval = nc_put_att_text(ncid, ps_id, varnom,
                                       strlen(var_unit), var_unit)))
           ERR(retval);
         if (do_2d) {
           strcpy(var_unit, "celsius") ;                  /* Surface temp.*/
           if ((retval = nc_put_att_text(ncid, st_id, varnom,
                                       strlen(var_unit), var_unit)))
             ERR(retval);
         }
         if (do_3d) {
           strcpy(var_unit, "celsius") ;                  /* Atmos. temp.*/
           if ((retval = nc_put_att_text(ncid, te_id, varnom,
                                       strlen(var_unit), var_unit)))
             ERR(retval);
           strcpy(var_unit, "mole^-1") ;                  /* Atmos. temp.*/
           if ((retval = nc_put_att_text(ncid, hu_id, varnom,
                                       strlen(var_unit), var_unit)))
             ERR(retval);
           strcpy(var_unit, "knots") ;                  /* Zonal wind */
           if ((retval = nc_put_att_text(ncid, uu_id, varnom,
                                       strlen(var_unit), var_unit)))
             ERR(retval);
           if ((retval = nc_put_att_text(ncid, vv_id, varnom, /* V wind */
                                       strlen(var_unit), var_unit)))
             ERR(retval);
           strcpy(var_unit, "m^-3") ;               /* Number density */
           if ((retval = nc_put_att_text(ncid, tnd_id, varnom,
                                       strlen(var_unit), var_unit)))
             ERR(retval);
         }

         /* End define mode. */
         if ((retval = nc_enddef(ncid)))
           ERR(retval);

         /* Write the coordinate variable data. */
         if ((retval = nc_put_var_float(ncid, lat_id, plat)))
           ERR(retval);
         if ((retval = nc_put_var_float(ncid, lon_id, plon)))
           ERR(retval);

         nc_ct[0] = 1 ;
         if (do_3d)
            nc_ct[1] = nlev ;
         nc_ct[ndims-2] = nj ;
         nc_ct[ndims-1] = ni ;
         nc_start[1] = 0;
         nc_start[2] = 0;
         nc_start[3] = 0;
         /** End netCDF file creation and initialization **/

       } /*** End of record information parameters' retrievals @ ct==0 ***/

       /**** Read the data from RPN file into memory ****/
       /* Always read the surface pressure data */
       strcpy(nomvar, "P0") ;
       ier = c_fstlir(ppres, iun3, &ni, &nj, &nk, -1, etiket, -1, -1,
                                                      -1, typvar, nomvar) ;
       /* Read 2-D vars */
       if (do_2d) {    
         strcpy(nomvar, "TG") ;
         ier = c_fstlir(ptsurf, iun, &ni, &nj, &nk, -1, etiket, -1, -1,
                                                      -1, typvar, nomvar) ;
       }

       /* Read 3-D vars, reading a 2-D slab at each vertical level */
       if (do_3d) {
         for (k = nlev-1; k >= 0; --k) {  
            ip1 = Level_ip1[k] ;
            strcpy(nomvar, "TT") ;
            ier = c_fstlir(ptemp+(k*surf_dim), iun3, &ni, &nj, &nk, -1, 
                                    etiket, ip1, -1, -1, typvar, nomvar) ;
            strcpy(nomvar, "HU") ;
            ier = c_fstlir(pwvap+(k*surf_dim), iun3, &ni, &nj, &nk, -1, 
                                    etiket, ip1, -1, -1, typvar, nomvar) ;
            strcpy(nomvar, "UU") ;
            ier = c_fstlir(puwind+(k*surf_dim), iun3, &ni, &nj, &nk, -1,
                                    etiket, ip1, -1, -1, typvar, nomvar) ;
            strcpy(nomvar, "VV") ;
            ier = c_fstlir(pvwind+(k*surf_dim), iun3, &ni, &nj, &nk, -1,
                                    etiket, ip1, -1, -1, typvar, nomvar) ;
         }

         /* Calculate the number density */
         int ij ;
         for (ij = 0; ij < surf_dim ; ++ij) {
           *(ptnd+(ij*(nlev-1))) = *(ppres+ij) * sig2[nlev-1] / ((*(ptemp+(ij*(nlev-1))) + 273.16) * bz) ;          
         }
         
       }  /* End reading 3-D rpn data*/

       /** Write each timestep to netCDF **/
       nc_start[0] = ct;

       if ((retval = nc_put_vara_float(ncid, ps_id, nc_start, nc_ct, ppres)))
	 ERR(retval);
       if (do_2d) {
         if ((retval = nc_put_vara_float(ncid, st_id, nc_start, nc_ct, ptsurf)))
	   ERR(retval);
       }
       if (do_3d) {
         if ((retval = nc_put_vara_float(ncid, te_id, nc_start, nc_ct, ptemp)))
	   ERR(retval);
         if ((retval = nc_put_vara_float(ncid, hu_id, nc_start, nc_ct, pwvap)))
	   ERR(retval);
         if ((retval = nc_put_vara_float(ncid, uu_id, nc_start, nc_ct, puwind)))
	   ERR(retval);
         if ((retval = nc_put_vara_float(ncid, vv_id, nc_start, nc_ct, pvwind)))
	   ERR(retval);

         if ((retval = nc_put_vara_float(ncid, tnd_id, nc_start, nc_ct, ptnd)))
	   ERR(retval);
       }


        /**** Close the opened file(s) with fst function  ****/
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

  } /* End the main loop */

  /* Close the netCDF file. */
  if ((retval = nc_close(ncid)))
     ERR(retval) ; 

  /* Free the allocated memory spaces */
  free (plat) ;
  free (plon) ;
  free (ppres) ;
  if (do_2d) {
    free (ptsurf) ;
  }
  if (do_3d) {
    free (ptemp) ;
    free (pwvap) ;
    free (puwind) ;
    free (pvwind) ;

    free (ptnd) ;
  } 

  return 0 ;
}
