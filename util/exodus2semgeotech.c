/* This program converts the Binary (provided that ncdump command exists)
or ASCII exodus file exported from the CUBIT/Trelis to several mesh files
required by the fem3d and sem3d_solid packages. Basic steps starting from
the CUBIT:
in CUBIT:
- blocks defines only the material regions
  -> actual material properties should not be defined within the CUBIT.
  actual material properties can be listed later corresponding to each block 
  (i.e., material region).
- nodal boundary conditions must be defined using node set
  -> each node set name must contain the corresponding BC names as defined in 
  char *ns_bcname[] below
  e.g., node set name can be front_nsbcux or front_nsbcux_nsbcuy etc.
- surface boundary conditions must be defined using side set
  -> each side set name must contain the corresponding BC names as defined in 
  char *ss_bcname[] below
  e.g., side set name can be front_ssbcux or front_ssbcux_ssbcuy etc.

step1: export mesh file as exodus file say "mesh.e"
step2: convert "mesh.e" to ASCII file using
  >>ncdump mesh.e > mesh.txt
step3: produce mesh and BC files
  >>exodus2semgeotech mesh.txt
  OR
  >>exodus2semgeotech mesh.txt 1000.0

There will be several output files:
*_coord_? : total number of nodes followed by nodal coordinate ? (? -> x, y, z)
*_connectivity : total number of elements followed by connectivity list
*_material_id : total number of elements followed by material IDs
*_??bcu? : node IDs which have u? = 0 as the boundary conditions (?? -> ns or ss, ? -> x, y, z) 
------------------------------------------------------
DEVELOPER:
  Hom Nath Gharti
  NORSAR
  homnath_AT_norsar_DOT_no
DEPENDENCY:
  stringmanip.c: string manipulation routines
COMPILE:
  gcc exodus2semgeotech.c -o exodus2semgeotech
USAGE: 
  exodus2semgeotech <inputfile> <OPTIONS>
  Example: exodus2semgeotech sloep3d_mest.txt
  or
  exodus2semgeotech slope3d_mesh.e -fac=0.001 -bin=1
OPTIONS:
  -fac: use this option to multiply coordinates. this is importantn for unit 
        conversion, e.g., to convert m to km use -fac=0.001
  -bin: use this option if you want to convert exodus binary directly, provided
        that the command ncdump is in the path. ncdump is a part of netCDF library
        that can be downloaded freely from 
        http://www.unidata.ucar.edu/downloads/netcdf/index.jsp. use -bin=1 for binary 
        or -bin=0 for ascii file.
HISTORY: 
  HNG,Apr 23,2010;HNG,Apr 17,2010;HNG,Feb 08,2009
TODO:
  - add support to side sets,i.e., surface boundary conditions
-------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "stringmanip.c"

#define OFF 0
#define ON 1

/* auxiliary routines */
void removeExtension(char *, char *);
int get_int(int *, char *, char *);
int look_int(int *, char *, char *);
int getfirstquote(char *, char *);

/* main routine */
int main(int argc,char **argv){
int i,itmp,j,k;
int ndim; /* geometry dimension */ 
int nnode,nelmt; /* number of nodes, number of elements */
int nblk,nns,nss; /* number of blocks, number of node sets */
int elmt_count,node_count; /* element, node count */
int node_countx,node_county,node_countz; /* element, node count */
int blk_count,ns_count,ss_count; /* block, node set count */
int ns_nbc; /* number of bc types in each node set */
int ss_nbc; /* number of bc types in each side set */
int dim_stat,ns_stat,ss_stat,con_stat,coord_stat,mat_stat; /* status */
int coordx_stat,coordy_stat,coordz_stat; /* status */
int *blk_nelmt,*blk_nenod; /* number of elements, number of nodes per element in each node */
int *ns_nnode,*ss_nside; /* number of nodes in each node set */

double fac,dtmp; /* multiplication factor for coordinates, temporary float */
char *bulk,line[100],token[62],dumc[250],stag[62];
double **coord;
char **coord_name; /* coordinates name */
char fonly[62],infname[62],outfname[62];
char **ns_name; /* node set names */
char **ss_name; /* side set names */
int ns_maxnbc; /* maximum number of BC types */
int ss_maxnbc; /* maximum number of BC types */

/* change this line to look for other BCs in node sets */
/* e.g., char *ns_bcname[4]={"nsbcux","nsbcuy","nsbcuz","nsbcfx"};*/
char *ns_bcname[]={"nsbcux","nsbcuy","nsbcuz"};
int *ns_bcfilestat;
int *ns_bc_nnode; /* number of nodes in each nodal bc */

/* change this line to look for other BCs in side sets */
/* e.g., char *ss_bcname[4]={"ssbcux","ssbcuy","ssbcuz","ssbcfx"};*/
char *ss_bcname[]={"ssbcux","ssbcuy","ssbcuz","ssbcphi","ssbcphix","ssbcphiy"};
int *ss_bcfilestat;
int *ss_elmt,*ss_side;
int *ss_bc_nside; /* number of sides in each side bc */

int isbin; /* test if binary */

FILE *inf,*outf_dum,*outf_mat,*outf_con,*outf_coord[3],**outf_nsbc,**outf_ssbc;

/* default factor and binary switch*/    
fac=1.0; isbin=OFF;

if(argc<2){
  fprintf(stderr,"ERROR: no input file!\n");
  exit(-1);
}

/* scan command */
if(argc>2){
  for(i=2;i<argc;i++){
    if(look_double(&dtmp,"-fac=",argv[i])==0){
      fac=dtmp;     
      continue;
    }
    else if(look_int(&itmp,"-bin=",argv[i])==0){
      isbin=itmp;
      continue;
    }else{
      printf("ERROR: unrecognized option \"%s\"\n",argv[i]);
      exit(-1);
    }
  }
}

printf("input file: %s\n",argv[1]);
printf("fac: %f\n",fac);
printf("isbin: %d\n",isbin);
printf("--------------------------------\n");
/* default input file name is argv[1]*/
strcpy(infname,argv[1]);
removeExtension(argv[1],fonly);

if (isbin){
  printf("converting binary to ascii...");
  /* set input file name */
  sprintf(infname,"%s.txt",fonly);

  /* convert binary netCDF file to ascii file */
  sprintf(dumc,"ncdump %s > %s.txt",argv[1],fonly);  
  if (system(dumc)!=0){
    printf("ERROR: command \"%s\" cannot be executed! use -bin=0 or no option for ascii input file! \n",dumc);
    exit(-1);
  }
  printf("complete!\n");
}

/* open input file */
inf=fopen(infname,"r");
if(inf==NULL){
  fprintf(stderr,"ERROR: file \"%s\" not found!\n",argv[1]);
  exit(-1);
}
/*printf("--------------------------------\n");*/

bulk=malloc(100000); /* bulk string */
      
/* initialize some variables to 0 */
ndim=0; nns=0; nblk=0; nnode=0; nelmt=0; nss=0;

/* intialize count to 0 */
blk_count=0; ns_count=0; ss_count=0; node_count=0; elmt_count=0;  
node_countx=0; node_county=0; node_countz=0;

/* set default status to OFF */
dim_stat=OFF; ns_stat=OFF; ss_stat=OFF; con_stat=OFF; coord_stat=OFF;
coordx_stat=OFF; coordy_stat=OFF; coordz_stat=OFF;

/* number of BC types to check for nodeset */
ns_maxnbc = sizeof(ns_bcname)/sizeof(char *);
outf_nsbc=malloc(ns_maxnbc*sizeof(FILE *));
ns_bcfilestat=malloc(ns_maxnbc*sizeof(int));
for(i=0;i<ns_maxnbc;i++)ns_bcfilestat[i]=0;
ns_bc_nnode=malloc(ns_maxnbc*sizeof(int));
for(j=0;j<ns_maxnbc;j++){     
  ns_bc_nnode[j]=0;
}

/* number of BC types to check for sideset */
ss_maxnbc = sizeof(ss_bcname)/sizeof(char *);
outf_ssbc=malloc(ss_maxnbc*sizeof(FILE *));
ss_bcfilestat=malloc(ss_maxnbc*sizeof(int));
for(i=0;i<ss_maxnbc;i++)ss_bcfilestat[i]=0;
ss_bc_nside=malloc(ss_maxnbc*sizeof(int));
for(j=0;j<ss_maxnbc;j++){     
  ss_bc_nside[j]=0;
}
fscanf(inf,"%s",token);
if(strcmp(token,"netcdf")!=0){
  printf("ERROR: invalid exodus file or wrong -bin value!\n");
  printf("HINT: try correct value for -bin option or use valid exodus file!\n");
  exit(-1);
}

while(!feof(inf)){
  fscanf(inf,"%s",token);
  /* read dimensions */
  if(dim_stat!=ON && strcmp(token,"dimensions:")==0){
    printf("reading dimensions...");
    while(strstr(fgets(line,100,inf),"variables:") == NULL){    
      strncat(bulk,line,strcspn(line,";")+1);              
    }    
    get_int(&ndim,"num_dim =",bulk);
    if(ndim>0){   
      /* allocate memory */
      coord_name=malloc(ndim*sizeof(char *));
      for(i=0;i<ndim;i++){      
      coord_name[i]=malloc(62*sizeof(char)); /* each name has maximum of 62 characters */   
      }   
    }else{
      printf("ERROR: illegal value of dimension!\n");
      exit(-1);
    }
    get_int(&nnode,"num_nodes =",bulk);
    get_int(&nelmt,"num_elem =",bulk);
    get_int(&nblk,"num_el_blk =",bulk);
    
    /* allocate memory */
    blk_nelmt=malloc(nblk*sizeof(int));
    blk_nenod=malloc(nblk*sizeof(int));
    /* nodeset information */
    if (look_int(&nns,"num_node_sets =",bulk)!=0){
      nns=0;
    }else{
      /* allocate memory */
      ns_name=malloc(nns*sizeof(char *));
      for(i=0;i<nns;i++){
        ns_name[i]=malloc(62*sizeof(char)); /* each name has maximum of 62 characters */
      }
      ns_nnode=malloc(nns*sizeof(int));
    }
    
    if(nns>0){ /* This segment has a significance only if nns has legitimate value */
      for(i=0;i<nns;i++){
        sprintf(stag,"num_nod_ns%d =",i+1);
        get_int(&ns_nnode[i],stag,bulk);          
      }
    }
    coord=malloc(3*sizeof(double *));
    for(i=0;i<3;i++){
        coord[i]=malloc(nnode*sizeof(double));
        if(coord[i]==NULL){
            printf("ERROR: not enough memory!");
            exit(-1);
        }
    }
  
    /* sideset information */
    if (look_int(&nss,"num_side_sets =",bulk)!=0){
      nss=0;
    }else{
      /* allocate memory */
      ss_name=malloc(nss*sizeof(char *));
      for(i=0;i<nss;i++){
        ss_name[i]=malloc(62*sizeof(char)); /* each name has maximum of 62 characters */
      }
      ss_nside=malloc(nss*sizeof(int));
    }
    
    if(nss>0){ /* This segment has a significance only if nss has legitimate value */
      for(i=0;i<nss;i++){
        sprintf(stag,"num_side_ss%d =",i+1);
        get_int(&ss_nside[i],stag,bulk);          
      }
    }

    /* block information */        
    for(i=0;i<nblk;i++){
      sprintf(stag,"num_el_in_blk%d =",i+1);
      get_int(&blk_nelmt[i],stag,bulk); 
      
      sprintf(stag,"num_nod_per_el%d =",i+1);
      get_int(&blk_nenod[i],stag,bulk);
    }    
      
    dim_stat=ON;
    free(bulk);
    printf("complete!\n");
    printf(" geometry dimension: %d\n",ndim);    
    printf(" number of nodes: %d\n",nnode);
    printf(" number of elements: %d\n",nelmt);
    printf(" number of blocks: %d\n",nblk);
    continue;
  }  
  
  /* read coordinate names */
  if(strcmp(token,"coor_names")==0){    
    fscanf(inf,"%s",dumc); /* = */      
    for (i=0; i<ndim; i++){
      fscanf(inf,"%s",dumc);    
      getfirstquote(dumc,coord_name[i]);    
    } 
    continue;
  }
  
  /* read and write nodal boundary conditions */
  if(strcmp(token,"ns_names")==0){
    printf("saving nodal BCs...");
    fscanf(inf,"%s",dumc); /* = */      
    for (i=0; i<nns; i++){
      fscanf(inf,"%s",dumc);    
      getfirstquote(dumc,ns_name[i]);

      /* count sides in each side BC */
      for(j=0;j<ns_maxnbc;j++){
        if (strstr(ns_name[i],ns_bcname[j])!=NULL){
          ns_bc_nnode[j]+=ns_nnode[i];
        }
      }
    }
   
    /* open bc nodal files */
    /*sprintf(outfname,"%s_bcux",fonly);
    outf_nsbc[0]=fopen(outfname,"w");
    sprintf(outfname,"%s_bcuy",fonly);
    outf_nsbc[1]=fopen(outfname,"w");
    sprintf(outfname,"%s_bcuz",fonly);
    outf_nsbc[2]=fopen(outfname,"w");*/      
    continue;
  }
  
  if(ns_stat!=ON){
    for(i=0;i<nns;i++){   
      sprintf(stag,"node_ns%d",i+1);
      if(strcmp(token,stag)==0){
        ns_nbc=0;
        for(j=0;j<ns_maxnbc;j++){   
          if (strstr(ns_name[i],ns_bcname[j])!=NULL){
            sprintf(outfname,"%s_%s",fonly,ns_bcname[j]);       
            if(ns_bcfilestat[j]==1){
              /* already opened */
              outf_nsbc[ns_nbc]=fopen(outfname,"a");
            }else{
              /* create new */
              outf_nsbc[ns_nbc]=fopen(outfname,"w");
              fprintf(outf_nsbc[ns_nbc],"%d\n",ns_bc_nnode[j]);
            }

            ns_bcfilestat[j]=1; /* this file is now opened */
            ns_nbc+=1;
          }
        } /* for(j=0 ..) */
    
        if(ns_nbc==0){
          printf("WARNING: no BC name found in node side \"%s\"!\n",ns_name[i]);
        }

        fscanf(inf,"%s",dumc); /* = */
        for(j=0;j<ns_nnode[i]; j++){      
          fscanf(inf,"%d,",&itmp); /* read comma separated data */
          for(k=0;k<ns_nbc;k++){        
            fprintf(outf_nsbc[k],"%d\n",itmp);
          }
        }

        for(j=0;j<ns_nbc;j++){
          fclose(outf_nsbc[j]);
        }

        /*if (strstr(ns_name[i],"bcux")!=NULL){
          ind_outf[ns_nbc]=0;
          ns_nbc+=1;
        }
        if (strstr(ns_name[i],"bcuy")!=NULL){
          ind_outf[ns_nbc]=1;
          ns_nbc+=1;
        }
        if (strstr(ns_name[i],"bcuz")!=NULL){
          ind_outf[ns_nbc]=2;
          ns_nbc+=1;
        }
        fscanf(inf,"%s",dumc);
        for(j=0; j<ns_nnode[i]; j++){       
          fscanf(inf,"%d,",&itmp);
          for(k=0;k<ns_nbc;k++){
            fprintf(outf_nsbc[ind_outf[k]],"%d\n",itmp);
          }
        }*/

        ns_count++;
        /* close bc files if writing has finished */
        if(ns_count==nns){
          /* free memory */
          for(i=0;i<nns;i++){
            free(ns_name[i]);
          }
          free(ns_bcfilestat);
          free(ns_nnode);
          ns_stat=ON;
            
          /*fclose(outf_nsbc[0]);
          fclose(outf_nsbc[1]);
          fclose(outf_nsbc[2]);*/
          printf("complete!\n");          
        } 
        continue;
      }          
    }
  }

  /* read and write side boundary conditions */
  if(strcmp(token,"ss_names")==0){
    printf("saving side BCs...");
    fscanf(inf,"%s",dumc); /* = */      
    for (i=0; i<nss; i++){  
      fscanf(inf,"%s",dumc);    
      getfirstquote(dumc,ss_name[i]);
      
      /* count sides in each side BC */
      for(j=0;j<ss_maxnbc;j++){     
        if (strstr(ss_name[i],ss_bcname[j])!=NULL){
          ss_bc_nside[j]+=ss_nside[i];
        }
      }
    }
 
    /* open bc nodal files */
    /*sprintf(outfname,"%s_bcux",fonly);
    outf_nsbc[0]=fopen(outfname,"w");
    sprintf(outfname,"%s_bcuy",fonly);
    outf_nsbc[1]=fopen(outfname,"w");
    sprintf(outfname,"%s_bcuz",fonly);
    outf_nsbc[2]=fopen(outfname,"w");*/      
    continue;
  }
 
  if(ss_stat!=ON){
    for(i=0;i<nss;i++){   
      sprintf(stag,"elem_ss%d",i+1);    
      if(strcmp(token,stag)==0){      
        ss_nbc=0;
        for(j=0;j<ss_maxnbc;j++){     
          if(strstr(ss_name[i],ss_bcname[j])!=NULL){ 
            sprintf(outfname,"%s_%s",fonly,ss_bcname[j]);       
            if(ss_bcfilestat[j]==1){
              /* already opened */
              outf_ssbc[ss_nbc]=fopen(outfname,"a");
            }else{
              /* create new */
              outf_ssbc[ss_nbc]=fopen(outfname,"w");
              fprintf(outf_ssbc[ss_nbc],"%d\n",ss_bc_nside[j]);
            }

            ss_bcfilestat[j]=1; /* this file is now opened */ 
            ss_nbc+=1;
          }
        }
    
        /*if(ss_nbc==0){
          printf("WARNING: no BC name found in sideset name \"%s\"!\n",ss_name[i]);
          // Open a filename with a name of sideset name 
          outf_dum=fopen(ss_name[i],"w");
        }*/

        ss_elmt=malloc(ss_nside[i]*sizeof(int));
        fscanf(inf,"%s",dumc); /* = */
        for(j=0;j<ss_nside[i]; j++){        
          fscanf(inf,"%d,",&ss_elmt[j]); /* read comma separated data */
          /*for(k=0;k<ns_nbc;k++){        
            fprintf(outf_nsbc[k],"%d\n",itmp);
          }*/
        }
        
        fscanf(inf,"%s",dumc); /* ; */    
        fscanf(inf,"%s",token);
        /*printf("%s\n",token);
        exit(-1)*/
        sprintf(stag,"side_ss%d",i+1);    
        if(strcmp(token,stag)==0){
          ss_side=malloc(ss_nside[i]*sizeof(int));
          fscanf(inf,"%s",dumc); /* = */
          for(j=0;j<ss_nside[i]; j++){        
            fscanf(inf,"%d,",&ss_side[j]); /* read comma separated data */
            /*for(k=0;k<ns_nbc;k++){        
            fprintf(outf_nsbc[k],"%d\n",itmp);
            }*/
          }
        }
        if(ss_nbc>0){
          /* write to appropriate BC files */
          for(k=0;k<ss_nbc;k++){
            for(j=0;j<ss_nside[i];j++){
              fprintf(outf_ssbc[k],"%d %d\n",ss_elmt[j],ss_side[j]);
            } 
          }
        }else{
          /* write to a dummy file with a name of ss_name */
          printf("WARNING: no BC name found in sideset name \"%s\"!\n",ss_name[i]);
          /* Open a filename with a name of sideset name */
          sprintf(outfname,"%s_%s",fonly,ss_name[i]);
          outf_dum=fopen(outfname,"w");
          /*outf_dum=fopen(ss_name[i],"w");*/
          fprintf(outf_dum,"%d\n",ss_nside[i]);
          for(j=0;j<ss_nside[i];j++){
            fprintf(outf_dum,"%d %d\n",ss_elmt[j],ss_side[j]);
          }
          fclose(outf_dum);
        }

        free(ss_elmt);
        free(ss_side);

        for(j=0;j<ss_nbc;j++){
          fclose(outf_ssbc[j]);
        }

        /*if (strstr(ns_name[i],"bcux")!=NULL){
          ind_outf[ns_nbc]=0;
          ns_nbc+=1;
        }
        if (strstr(ns_name[i],"bcuy")!=NULL){
          ind_outf[ns_nbc]=1;
          ns_nbc+=1;
        }
        if (strstr(ns_name[i],"bcuz")!=NULL){
          ind_outf[ns_nbc]=2;
          ns_nbc+=1;
        }
        fscanf(inf,"%s",dumc);
        for(j=0; j<ns_nnode[i]; j++){       
          fscanf(inf,"%d,",&itmp);
          for(k=0;k<ns_nbc;k++){
            fprintf(outf_nsbc[ind_outf[k]],"%d\n",itmp);
          }
        }*/

        ss_count++;
        /* close bc files if writing has finished */
        if(ss_count==nss){
          /* free memory */
          for(i=0;i<nss;i++){
            free(ss_name[i]);
          }
          free(ss_bcfilestat);      
          free(ss_nside);
          ss_stat=ON;
            
          /*fclose(outf_nsbc[0]);
          fclose(outf_nsbc[1]);
          fclose(outf_nsbc[2]);*/
          printf("complete!\n");          
        } 
        continue;
      }          
    }
  }
 
  /* Connectivity */
  if(nblk>0 && con_stat!=ON){   
   
    /* write connectivity and material id */  
    for(i=0;i<nblk;i++){
      sprintf(stag,"connect%d",i+1);
      if(strcmp(token,stag)==0){
        blk_count++;
        
        /* open connectivity and material files */
        if(blk_count==1){          
          printf("saving connectivity and materials..."); 
          sprintf(outfname,"%s_connectivity",fonly);
          outf_con=fopen(outfname,"w");          
          fprintf(outf_con,"%d\n",nelmt);
          
          sprintf(outfname,"%s_material_id",fonly);
          outf_mat=fopen(outfname,"w");       
          fprintf(outf_mat,"%d\n",nelmt);
        }
        
        fscanf(inf,"%s",dumc); /* = */
        for(j=0;j<blk_nelmt[i];j++){
          for(k=0;k<blk_nenod[i];k++){
            fscanf(inf,"%d,",&itmp);
            fprintf(outf_con,"%d ",itmp);       
          }          
          elmt_count++;
          fprintf(outf_con,"\n");/* new line */
          fprintf(outf_mat,"%d\n",i+1);
        }        
        
        if(blk_count==nblk){
          con_stat=ON;
          mat_stat=ON;
          free(blk_nelmt);
          free(blk_nenod);
          printf("complete!\n");
          fclose(outf_con);
          fclose(outf_mat);
        }        
        continue;
      }        
    }
  }   
 
 /* Coordinates */
   if(strcmp(token,"coordx")==0){
        printf("reading x coordinates...");
        fscanf(inf,"%s",dumc);
        for(j=0;j<nnode;j++){
            fscanf(inf,"%lf,",&dtmp); /* read comma separated data */ 
            coord[0][j]=fac*dtmp;
            node_countx++;
        }   
        coordx_stat=ON;
        printf("complete!\n");
        continue;
    }   
    if(strcmp(token,"coordy")==0){
        printf("reading y coordinates...");
        fscanf(inf,"%s",dumc);
        for(j=0;j<nnode;j++){
            fscanf(inf,"%lf,",&dtmp); /* read comma separated data */ 
            coord[1][j]=fac*dtmp;
            node_county++;
        }   
        coordy_stat=ON;
        printf("complete!\n");
        continue;
    }   
    if(strcmp(token,"coordz")==0){
        printf("reading z coordinates...");
        fscanf(inf,"%s",dumc);
        for(j=0;j<nnode;j++){
            fscanf(inf,"%lf,",&dtmp); /* read comma separated data */ 
            coord[2][j]=fac*dtmp;
            node_countz++;
        }   
        coordz_stat=ON;
        printf("complete!\n");
        continue;
    }
}
  /* write coordinates file */
    printf("writing coordinates...");
    for(i=0;i<ndim;i++){
        sprintf(outfname,"%s_coord_%s",fonly,coord_name[i]);
        outf_coord[i]=fopen(outfname,"w");
        fprintf(outf_coord[i],"%d\n",nnode);
        for(j=0;j<nnode;j++)fprintf(outf_coord[i],"%.15lf\n",coord[i][j]);
        fclose(outf_coord[i]);
    }
    printf("complete!\n");
    for(i=0;i<ndim;i++){
        free(coord_name[i]);
    }
    free(coord_name);
    for(i=0;i<3;i++)free(coord[i]);
    free(coord);
  /* Coordinates */
  /*if(strcmp(token,"coord")==0){
    printf("saving coordinates...");
    fscanf(inf,"%s",dumc);      
    for(i=0;i<ndim;i++){
      sprintf(outfname,"%s_coord_%s",fonly,coord_name[i]);
      outf_coord[i]=fopen(outfname,"w");
      fprintf(outf_coord[i],"%d\n",nnode);
      for(j=0;j<nnode;j++){
        fscanf(inf,"%lf,",&ftmp);*/ /* read comma separated data */
    /*    fprintf(outf_coord[i],"%.15lf\n",fac*ftmp);
        if(i==0)node_count++;
      }
      fclose(outf_coord[i]);
    }
    f(r(i=0;i<ndim;i++){
      free(coord_name[i]);
    }
  
    coord_stat=ON;
    printf("complete!\n");
    continue;
  } 
 */ 

/* check status */
if(coordx_stat!=ON && coordx_stat!=ON && coordx_stat!=ON){
  printf("ERROR: coordinates cannot be read!\n");
  printf("HINT: it seems that the EXODUS file format is \"set large exodus file off\"!\n");
  printf("      Please try using exodusold2semgeotech instead!\n");
  exit(-1);
}
if(nnode!=node_countx || nnode!=node_county || nnode!=node_countz){
  printf("ERROR: number of nodes inconsistent!\n");
  printf("%d %d %d %d\n",nnode,node_countx,node_county,node_countz);
  exit(-1);
}
if(nelmt!=elmt_count){
  printf("ERROR: number of elements inconsistent!\n");
  exit(-1);
}
if(dim_stat!=ON){
  printf("ERROR: dimensions cannot be read!\n");
  exit(-1);
}
if(ns_stat!=ON){
  printf("WARNING: nodal boundary conditions cannot be read!\n");
}
if(ss_stat!=ON){
  printf("WARNING: side boundary conditions cannot be read!\n");
}
if(con_stat!=ON){
  printf("ERROR: connectivity cannot be read!\n");
  exit(-1);
}
if(mat_stat!=ON){
  printf("ERROR: material IDs cannot be read!\n");
  exit(-1);
}
if(coordx_stat!=ON){
  printf("ERROR: x coordinates cannot be read!\n");
  exit(-1);
}
if(coordy_stat!=ON){
  printf("ERROR: y coordinates cannot be read!\n");
  exit(-1);
}
if(coordz_stat!=ON){
  printf("ERROR: z coordinates cannot be read!\n");
  exit(-1);
}
/*if(coord_stat!=ON){
  printf("ERROR: coordinates cannot be read!\n");
  exit(-1);
}
*/
printf("--------------------------------\n");
return(0);
}
/*============================================================================*/
