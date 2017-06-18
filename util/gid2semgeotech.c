/** @file gid2semgeotech.c
*  @brief Converts ASCII Gid mesh file to SPECFEM3D_GEOTECH files.
*
*  This program converts the ASCII GiD mesh file to several mesh files required
*  by the SPECFEM3D_GEOTECH package. GiD (www.gidhome.com) is a commercial pre
*  and post processor for numerical simulations.
*
*  @author Hom Nath Gharti (hgharti_AT_princeton_DOT_edu), Zhenzhen Yan

* ## Dependencies:
*  stringmanip.c: string manipulation routines
*
* ## Compile:
*  gcc gid2semgeotech.c -o gid2semgeotech
*  
* ## Usage:
*  gid2semgeotech <inputfile> <OPTIONS>
*  Example: gid2semgeotech gid2semgeotech_example.dat
*  or
*  gid2semgeotech gid2semgeotech_example.dat -fac=0.001
*  
* ## Options:
* - -fac: use this option to multiply coordinates. this is importantn for unit
*        conversion, e.g., to convert m to km use -fac=0.001
* # Basic steps starting from GID:
*
* ### step1: export mesh file in ASCII format "mesh.dat"
* 
* ### step2: produce mesh and BC files
*  >>gid2semgeotech mesh.dat
*  OR
*  >>gid2semgeotech mesh.dat 1000.0
*
*There will be several output files:
* - coord_? : total number of nodes followed by nodal coordinate ? (? -> x, y, z)
*
* - _connectivity : total number of elements followed by connectivity list
*
* - _material_id : total number of elements followed by material IDs
*  
* - ??bcu? : node IDs which have u? = 0 as the boundary conditions (?? -> ns or ss, ? -> x, y, z)
*/
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
int ielmt,iface,imat,inode,ix,iy,iz,n1,n2,n3,n4,n5,n6,n7,n8;
int nelmt_bc,nface,nbcx,nbcy,nbcz;
int ndim; /* geometry dimension */
int nnode,nelmt; /* number of nodes, number of elements */
int nblk,nns,nss; /* number of blocks, number of node sets */
int elmt_count,node_count; /* element, node count */
int blk_count,ns_count,ss_count; /* block, node set count */
int ns_nbc; /* number of bc types in each node set */
int ss_nbc; /* number of bc types in each side set */
int dim_stat,ns_stat,ss_stat,con_stat,coord_stat,mat_stat; /* status */
int *blk_nelmt,*blk_nenod; /* number of elements, number of nodes per element in each node */
int *ns_nnode,*ss_nside; /* number of nodes in each node set */

double fac,ftmp,x,y,z; /* multiplication factor for coordinates, temporary float */
char *bulk,line[100],token[62],dumc[62],stag[62];
char **coord_name; /* coordinates name */
char fonly[62],infname[62],outfname[62],outfnamex[62],outfnamey[62],outfnamez[62];
char **ns_name; /* node set names */
char **ss_name; /* node set names */
int ns_maxnbc; /* maximum number of BC types */
int ss_maxnbc; /* maximum number of BC types */

/* change this line to look for other BCs in node sets */
/* e.g., char *ns_bcname[4]={"nsbcux","nsbcuy","nsbcuz","nsbcfx"};*/
char *ns_bcname[]={"nsbcux","nsbcuy","nsbcuz"};
int *ns_bcfilestat;
int *ns_bc_nnode; /* number of nodes in each nodal bc */

/* change this line to look for other BCs in side sets */
/* e.g., char *ss_bcname[4]={"ssbcux","ssbcuy","ssbcuz","ssbcfx"};*/
char *ss_bcname[]={"ssbcux","ssbcuy","ssbcuz"};
int *ss_bcfilestat;
int *ss_elmt,*ss_side;
int *ss_bc_nside; /* number of sides in each side bc */

int gid2exodus[]={5,4,1,2,3,6};
int imatnum,idomain;
double gamma,ym,nu,phi,coh,psi;
int isbin; /* test if binary */

FILE *inf,*tempf,*outf_dum,*outf_mat,*outf_con,*outf_coord[3],*outf_bc[3],**outf_nsbc,**outf_ssbc;

/* default factor and binary switch*/
fac=1.0; isbin=OFF;

if(argc<2){
  fprintf(stderr,"ERROR: input file not entered!\n");
  exit(-1);
}

/* scan command */
if(argc>2){
  for(i=2;i<argc;i++){
    if(look_double(&ftmp,"-fac=",argv[i])==0){
    fac=ftmp;
    continue;
  }else{
    printf("ERROR: unrecognized option \"%s\"",argv[i]);
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

/* open input file */
inf=fopen(infname,"r");
if(inf==NULL){
  fprintf(stderr,"ERROR: file \"%s\" not found!",argv[1]);
  exit(-1);
}

bulk=malloc(1000); /* bulk string */

/* initialize some variables to 0 */
ndim=0; nns=0; nblk=0; nnode=0; nelmt=0; nss=0;

/* intialize count to 0 */
blk_count=0; ns_count=0; ss_count=0; node_count=0; elmt_count=0;

/* set default status to OFF */
dim_stat=OFF; ns_stat=OFF; ss_stat=OFF; con_stat=OFF; coord_stat=OFF;

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

while(!feof(inf)){
  fgets(line,100,inf);
  sscanf(line,"%s",token);
  /* read problem size */
  if(dim_stat!=ON && strcmp(token,"Number")==0){
    fscanf(inf,"%d %d %d",&nelmt,&nnode,&nblk);     
  
    printf(" number of elements: %d\n",nelmt);
    printf(" number of nodes: %d\n",nnode);
    printf(" number of blocks: %d\n",nblk);
    dim_stat=ON;     
    continue;
  }
  /* read and save coordinates */
  if(strcmp(token,"$")==0 && strstr(line,"Coordinates")!=NULL){
    printf("saving coordinates...");
    fgets(line,100,inf);/* discard this line */
    sprintf(outfnamex,"%s_coord_x",fonly);
    sprintf(outfnamey,"%s_coord_y",fonly);
    sprintf(outfnamez,"%s_coord_z",fonly);
    outf_coord[0]=fopen(outfnamex,"w");
    outf_coord[1]=fopen(outfnamey,"w");
    outf_coord[2]=fopen(outfnamez,"w");
    fprintf(outf_coord[0],"%d\n",nnode);
    fprintf(outf_coord[1],"%d\n",nnode);
    fprintf(outf_coord[2],"%d\n",nnode);
    for(i=0;i<nnode;i++){
      fscanf(inf,"%d %lf %lf %lf",&inode,&x,&y,&z); /* read comma separated data */
      fprintf(outf_coord[0],"%.15lf\n",fac*x);
      fprintf(outf_coord[1],"%.15lf\n",fac*y);
      fprintf(outf_coord[2],"%.15lf\n",fac*z);
      node_count++;
    }
    fclose(outf_coord[0]);
    fclose(outf_coord[1]);
    fclose(outf_coord[2]);
    
    coord_stat=ON;
    printf("complete!\n");
    continue;
  }

  /* Connectivity */
  if(strcmp(token,"$")==0 && strstr(line,"$ Con")!=NULL){
    printf("saving connectivity & material IDs...");
    fgets(line,100,inf);/* discard this line */
    sprintf(outfname,"%s_connectivity",fonly);
    outf_con=fopen(outfname,"w");
    fprintf(outf_con,"%d\n",nelmt);
    sprintf(outfname,"%s_material_id",fonly);
    outf_mat=fopen(outfname,"w");
    fprintf(outf_mat,"%d\n",nelmt);
    for(i=0;i<nelmt;i++){
      fscanf(inf,"%d %d %d %d %d %d %d %d %d %d",&ielmt,
      &n1,&n2,&n3,&n4,&n5,&n6,&n7,&n8,&imat);
      fprintf(outf_con,"%d %d %d %d %d %d %d %d\n",n1,n2,n3,n4,n5,n6,n7,n8);
      fprintf(outf_mat,"%d\n",imat);
      elmt_count++;
    }
    con_stat=ON;
    mat_stat=ON;
    printf("complete!\n");
    fclose(outf_con);
    fclose(outf_mat);
    continue;
  }

  /* Material list */
  if(strcmp(token,"$")==0 && strstr(line,"$ Materials")!=NULL){
    printf("saving material list...");
    fgets(line,100,inf);/* discard this line */
    sprintf(outfname,"%s_material_list",fonly);
    outf_mat=fopen(outfname,"w");
    fprintf(outf_mat,"# material properties (id,domain,gamma,ym,nu,phi,coh,psi)\n");
    fprintf(outf_mat,"%d\n",nblk);
    fscanf(inf,"%d %d %lf %lf %lf %lf %lf %lf",&imatnum,&idomain,&gamma,&ym,&nu,&phi,&coh,&psi);
    fprintf(outf_mat,"%d %d %.6f %.6e %.6f %.6f %.6f %.6f\n",imatnum,idomain,gamma,ym,nu,phi,coh,psi);
    printf("complete!\n");
    fclose(outf_mat);
    continue;
  }

  /* Boundary conditions */
  if(strcmp(token,"$")==0 && strstr(line,"$ Boundary")!=NULL){
    printf("saving boundary conditions...");
    fgets(line,100,inf);/* discard this line */
    
    tempf=fopen("temp_bc","w");
    nelmt_bc=0;
    nbcx=0;
    nbcy=0;
    nbcz=0;
    while(fscanf(inf,"%d %d %d %d %d",&ielmt,&iface,&ix,&iy,&iz)>0){
      nelmt_bc++;
      if(ix==1)nbcx++;
      if(iy==1)nbcy++;
      if(iz==1)nbcz++;
      fprintf(tempf,"%d %d %d %d %d\n",ielmt,gid2exodus[iface-1],ix,iy,iz);
    }
    fclose(tempf);
    sprintf(outfnamex,"%s_ssbcux",fonly);
    sprintf(outfnamey,"%s_ssbcuy",fonly);
    sprintf(outfnamez,"%s_ssbcuz",fonly);
    outf_bc[0]=fopen(outfnamex,"w");
    outf_bc[1]=fopen(outfnamey,"w");
    outf_bc[2]=fopen(outfnamez,"w");
    fprintf(outf_bc[0],"%d\n",nbcx);
    fprintf(outf_bc[1],"%d\n",nbcy);
    fprintf(outf_bc[2],"%d\n",nbcz);
    tempf=fopen("temp_bc","r");
    for(i=0;i<nelmt_bc;i++){
      fscanf(tempf,"%d %d %d %d %d",&ielmt,&iface,&ix,&iy,&iz);
      if(ix==1)fprintf(outf_bc[0],"%d %d\n",ielmt,iface);
      if(iy==1)fprintf(outf_bc[1],"%d %d\n",ielmt,iface);
      if(iz==1)fprintf(outf_bc[2],"%d %d\n",ielmt,iface);
    }
    fclose(tempf);
    unlink("temp_bc");
    fclose(outf_bc[0]);    
    fclose(outf_bc[1]);    
    fclose(outf_bc[2]);    
    printf("complete!\n");
    continue;
exit(0);
  }
}
/* check status */
if(nnode!=node_count){
  printf("ERROR: number of nodes inconsistent!\n");
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
  exit(-1);
}
if(ss_stat!=ON){
  printf("WARNING: side boundary conditions cannot be read!\n");
  exit(-1);
}
if(con_stat!=ON){
  printf("ERROR: connectivity cannot be read!\n");
  exit(-1);
}
if(mat_stat!=ON){
  printf("ERROR: material IDs cannot be read!\n");
  exit(-1);
}
if(coord_stat!=ON){
  printf("ERROR: coordinates cannot be read!\n");
  exit(-1);
}

printf("--------------------------------\n");
return(0);
}
/*======================================*/
