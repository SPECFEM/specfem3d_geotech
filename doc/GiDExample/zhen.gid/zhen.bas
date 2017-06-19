==================================================================
                        General Data File
==================================================================

%%%%%%%%%%%%%%%%%% Problem Size  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Number of Elements & Nodes , and nmats:
*nelem *npoin *nmats

%%%%%%%%%%%%%%%%%%% Mesh Database  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
$ Coordinates
$  NODO  COOR.-X   COOR.-Y   COOR.-Z
*loop nodes
*format "%6i %12.6e %12.6e %12.6e"
  *NodesNum *NodesCoord
*end nodes
.................................................................
$ Conectivities
$ ELEM. SECUENCIA DE CONECTIVIDADES MATER.
*loop elems
*elemsnum *elemsConec *elemsmat
*end elems
.................................................................
$ Materials
$ MatNum domain gamma ym nu phi coh psi
*loop materials
   *matnum() *MatProp(1) *MatProp(2) *MatProp(3) *MatProp(4) *MatProp(5) *MatProp(6) *MatProp(7)
*end
.................................................................
$ Boundary Condition
*Set Cond Line-Constraints *nodes *or(1,int) *or(3,int)
*Add Cond Point-Constraints *nodes *or(1,int) *or(3,int)
*loop nodes *OnlyInCond
*NodesNum *cond(1) *cond(2) *cond(3)
*end nodes
$ ELEMNUM  FaceID  X  Y  Z
*Set elems(tetrahedra)
*Set Cond Surface-constraints *elems *or(1,int) *or(2,int) *or(3,int) *CanRepeat
*loop elems *OnlyInCond
*if(CondNumEntities(int)>0)
*format "%i %i %i %i %i"
*ElemsNum *CondElemFace *cond(1) *cond(2) *cond(3)
*endif
*end elems
*Set elems(hexahedra)
*Set Cond Surface-constraints *elems *or(1,int) *or(2,int) *or(3,int) *CanRepeat
*loop elems *OnlyInCond
*if(CondNumEntities(int)>0)
*format "%i %i %i %i %i"
*ElemsNum *CondElemFace *cond(1) *cond(2) *cond(3)
*endif
*end elems