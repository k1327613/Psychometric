 # Psychometric model, code by Michel Nivard

rm(list = ls())

# Load Library
require(OpenMx)
require(psych)

# ------------------------------------------------------------------------------
# PREPARE DATA

allVars <- c('ID_fam','zyg','sexzyg','sex1','sex2','cchaos91', 'cchaos92', 'cchaos121', 'cchaos122', 'cchaos141','cchaos142', 'cchaos161', 'cchaos162', 'pchaos9', 'pchaos12', 'pchaos14', 'cpfee91', 'cpfee92', 'cpfee121', 'cpfee122', 'cpfee141', 'cpfee142', 'cpfeeneg91', 'cpfeeneg92', 'cpfeeneg121', 'cpfeeneg122', 'cpfeeneg141', 'cpfeeneg142', 'cpfeepos91', 'cpfeepos92', 'cpfeepos121', 'cpfeepos122', 'cpfeepos141', 'cpfeepos142', 'ppfee91', 'ppfee92', 'ppfee121', 'ppfee122', 'ppfee141', 'ppfee142', 'ppfeeneg91', 'ppfeeneg92', 'ppfeeneg121', 'ppfeeneg122', 'ppfeeneg141', 'ppfeeneg142', 'ppfeepos91', 'ppfeepos92', 'ppfeepos121', 'ppfeepos122', 'ppfeepos141', 'ppfeepos142',  'cpdis91', 'cpdis92', 'cpdis121', 'cpdis122', 'cpdis141', 'cpdis142',  'cpdis161', 'cpdis162', 'cpdisneg91', 'cpdisneg92', 'cpdisneg121', 'cpdisneg122', 'cpdisneg141', 'cpdisneg142', 'cpdisneg161', 'cpdisneg162', 'cpdispos91', 'cpdispos92', 'cpdispos121', 'cpdispos122', 'cpdispos141', 'cpdispos142', 'cpdispos161', 'cpdispos162','ppdis91', 'ppdis92', 'ppdis121', 'ppdis122', 'ppdis141', 'ppdis142', 'ppdisneg91', 'ppdisneg92', 'ppdisneg121', 'ppdisneg122', 'ppdisneg141', 'ppdisneg142', 'ppdispos91', 'ppdispos92', 'ppdispos121', 'ppdispos122', 'ppdispos141', 'ppdispos142', 'cmfq121', 'cmfq122', 'cbhmfq161', 'cbhmfq162', 'pmfq121', 'pmfq122', 'pbhmfq161', 'pbhmfq162' )				#Name your variables
data <- read.table ('/Users/lauriejhannigan/Dropbox/Documents (Dropbox)/PhD/My Documents (KCL)/Projects/TEDS/Data/ALLDATA_all_ages_by reporter_for_MX_v2.raw', header=FALSE, na.strings='.', col.names=allVars)	#THIS DIRECTS OPENMX TO YOUR DATAFILE IN THE WORKING DIRECTORY

selVars  <- c('cpdis91','ppdis91','cpdis121','ppdis121','cpdis141','ppdis141','cpdis92', 'ppdis92','cpdis122', 'ppdis122','cpdis142', 'ppdis142')  #SELECT VARIABLES FOR ANALYSIS


nv <- 6
nlv <- nv/2
ntv <- 2*nv

mzmData <- subset(data, sexzyg==1, selVars)	#IDENTIFIES SUBGROUPS BASED ON ZYGOSITY - in this dataset use sexzyg
dzmData <- subset(data, sexzyg==2, selVars)	#1MZM 2DZM 3MZF 4DZF 5DZO 7unknown
mzfData <- subset(data, sexzyg==3, selVars)
dzfData <- subset(data, sexzyg==4, selVars)
dzmfData <- subset(data, sexzyg==5, selVars)

MxdataMZm   <-mxData( observed= mzmData, type="raw" )
MxdataDZm   <-mxData( observed= dzmData, type="raw" )
MxdataMZf   <-mxData( observed= mzfData, type="raw" )
MxdataDZf   <-mxData( observed= dzfData, type="raw" )
MxdataDZmf   <-mxData( observed= dzmfData, type="raw" )

# Algebra for expected Mean and Variance/Covariance Matrices in MZ & DZ twins
MeansMZm		<-mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values= .2, label=c("mm1","mm1","mm2", "mm2","mm3","mm3"), name="Mm" )
MeansMZmm		<-mxAlgebra( expression= cbind(Mm,Mm), name="MMZm")
MeansMZf		<-mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values= .2, label=c("mf1","mf1","mf2", "mf2","mf3","mf3"), name="Mf" )
MeansMZff		<-mxAlgebra( expression= cbind(Mf,Mf), name="MMZf")
MeansDZm		<-mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values= .2, label=c("mm1","mm1","mm2", "mm2","mm3","mm3"), name="Mm" )
MeansDZmm		<-mxAlgebra( expression= cbind(Mm,Mm), name="MDZm")
MeansDZf		<-mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values= .2, label=c("mf1","mf1","mf2", "mf2","mf3","mf3"), name="Mf" )
MeansDZff		<-mxAlgebra( expression= cbind(Mf,Mf), name="MDZf")
MeansDZmf 		<-mxAlgebra(expression= cbind(Mm,Mf), name="MDZmf")

### Here we specify the a,c and e path for the rater agreement ###
a <- mxMatrix("Lower", nrow=nlv, ncol=nlv, free= T, values= .1 , 	  labels= c(   "a11"
																																	,"a21","a22"
																																	,"a31","a32","a33"), byrow=T,  name="X")

c <- mxMatrix("Lower", nrow=nlv, ncol=nlv, free= T, values= .1,    labels= c(     "c11"
																																	,"c21","c22"
																																	,"c31","c32","c33"), byrow=T,  name="Y")  
																							 								
																							 								
e <- mxMatrix("Lower", nrow=nlv, ncol=nlv, free= T, values= .1,  labels= c(  "e11"
																																,"e21","e22"
																																,"e31","e32","e33"), byrow=T,  name="Z") 

### Here we square the a,c and e path for the rater agreement to obtain the varience of A,C and E  ###
VA <- mxAlgebra(X %*% t(X), name="A")
VC <- mxAlgebra(Y %*% t(Y), name="C")
VE <- mxAlgebra(Z %*% t(Z), name="E")  

# squared paths
a2 <- mxAlgebra( expression=X * X, name='a2')
c2 <- mxAlgebra( expression=Y * Y, name='c2')
e2 <- mxAlgebra( expression=Z * Z, name='e2')

### Here we specify the loadings of the latent "agreement" varaible on the observed variables ###
Lambda <- mxMatrix(type="Full", nrow=ntv, ncol=nv, free=F,values=			c( 1,0,0,0,0,0,
																																		1,0,0,0,0,0,
																															 			0,1,0,0,0,0,
																																 		0,1,0,0,0,0,
			  																													 		0,0,1,0,0,0,
			  																															0,0,1,0,0,0,
			  																															0,0,0,1,0,0,
			  																															0,0,0,1,0,0,
				  																												  		0,0,0,0,1,0,
				  																												  		0,0,0,0,1,0,
				  																												  		0,0,0,0,0,1,
																											  					  		0,0,0,0,0,1), byrow=T, name="Ly")

#Matrices for scalar - creating a scalar for males (Agreement & child-/parent-report)
scalarA    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=1, label=c("k","l","m","k","l","m"), name="scalarA" ) 
scalarAmf  <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=c(T,T,T,F,F,F), values=1, label=c("k","l","m","one","one","one"), name="scalarAmf" )  #create a scalar males if DZ pairs are ordered male, female

scalarC    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=1, label=c("kc","lc","mc","kc","lc","mc"), name="scalarC" ) 
scalarCmf  <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=c(T,T,T,F,F,F), values=1, label=c("kc","lc","mc","one","one","one"), name="scalarCmf" )  #create a scalar males if DZ pairs are ordered male, female

scalarP    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=1, label=c("kp","lp","mp","kp","lp","mp"), name="scalarP" ) 
scalarPmf  <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=c(T,T,T,F,F,F), values=1, label=c("kp","lp","mp","one","one","one"), name="scalarPmf" )  #create a scalar males if DZ pairs are ordered male, female


### Here we create the latent covariance matrix for MZ twins ####
PsyMZf <-  mxAlgebra(rbind(cbind(A+C+E   , A+C),
                          cbind(A+C     , A+C+E)), name="PsyMZf")

PsyMZm <-  mxAlgebra( scalarA %*% (rbind( cbind(A+C+E , A+C),
                                cbind(A+C    , A+C+E))), name="PsyMZm")

### Here we create the latent covariance matrix for DZ twins ####
PsyDZf <-  mxAlgebra(rbind(cbind( A+C+E , (.5%x%A)+C),
                          cbind((.5%x%A)+C        , A+C+E)), name="PsyDZf")
                          
PsyDZm <-  mxAlgebra( scalarA %*% (rbind( cbind(A+C+E ,(.5%x%A)+C),
                                cbind((.5%x%A)+C   , A+C+E))), name="PsyDZm")
                          
PsyDZmf <-  mxAlgebra( scalarAmf %*% (rbind( cbind(A+C+E , (.5%x%A)+C),
                                cbind((.5%x%A)+C   , A+C+E))), name="PsyDZmf")                                                    

### Here we specify the a,c and e paths for the rater DISagreement for the child ratings  ###
ap <- mxMatrix("Lower", nrow=nlv, ncol=nlv, free= T, values= .2 , 	  labels= c("ap11"
																																	,"ap21","ap22"
																																	,"ap31","ap32","ap33"), byrow=T,  name="Xp")

cp <- mxMatrix("Lower", nrow=nlv, ncol=nlv, free= T, values= .2,    labels= c( "cp11"
																																	,"cp21","cp22"
																																	,"cp31","cp32","cp33"), byrow=T,  name="Yp")  
																					 																														 								
ep <- mxMatrix("Lower", nrow=nlv, ncol=nlv, free= T, values= .1,  labels= c("ep11"
																																,"ep21","ep22"
																																,"ep31","ep32","ep33"), byrow=T,  name="Zp") 

### Here we square the a,c and e path for the rater DISagreement to obtain the varience of A,C and E for the child ratings  ###
VAp <- mxAlgebra(Xp %*% t(Xp), name="Ap")
VCp <- mxAlgebra(Yp %*% t(Yp), name="Cp")
VEp <- mxAlgebra(Zp %*% t(Zp), name="Ep")  

ap2 <- mxAlgebra( expression=Xp * Xp, name='ap2')
cp2 <- mxAlgebra( expression=Yp * Yp, name='cp2')
ep2 <- mxAlgebra( expression=Zp * Zp, name='ep2')

### Here we specify the a,c and e paths for the rater DISagreement for the child ratings  ###
ac <- mxMatrix("Lower", nrow=nlv, ncol=nlv, free= T, values= .2 , 	  labels= c("ac11"
																																	,"ac21","ac22"
																																	,"ac31","ac32","ac33"), byrow=T,  name="Xc")

cc <- mxMatrix("Lower", nrow=nlv, ncol=nlv, free= T, values= .2,    labels= c( "cc11"
																																	,"cc21","cc22"
																																	,"cc31","cc32","cc33"), byrow=T,  name="Yc")  
																			 																															 								
ec <- mxMatrix("Lower", nrow=nlv, ncol=nlv, free= T, values= .1,  labels= c("ec11"
																																,"ec21","ec22"
																																,"ec31","ec32","ec33"), byrow=T,  name="Zc") 

### Here we square the a,c and e path for the rater DISagreement to obtain the variance of A,C and E for the child ratings  ###
VAc <- mxAlgebra(Xc %*% t(Xc), name="Ac")
VCc <- mxAlgebra(Yc %*% t(Yc), name="Cc")
VEc <- mxAlgebra(Zc %*% t(Zc), name="Ec")  

ac2 <- mxAlgebra( expression=Xc * Xc, name='ac2')
cc2 <- mxAlgebra( expression=Yc * Yc, name='cc2')
ec2 <- mxAlgebra( expression=Zc * Zc, name='ec2')

																	 

### Here we create and expand (using Ly) the child-report part of the residual matrix for MZ twins  ####
cMZf <-  mxAlgebra( Ly%&%(rbind(cbind(Ac+Cc+Ec   , Ac+Cc),
                          cbind(Ac+Cc     , Ac+Cc+Ec))), name="cMZf")

cMZm <-  mxAlgebra( Ly%&%(scalarC %*% (rbind( cbind(Ac+Cc+Ec , Ac+Cc),
                                cbind(Ac+Cc    , Ac+Cc+Ec)))), name="cMZm")

### Here we create expand (using Ly) the child-report part of the residual matrix for DZ twins ####
cDZf <-  mxAlgebra(Ly%&%(rbind(cbind( Ac+Cc+Ec , (.5%x%Ac)+Cc),
                          cbind((.5%x%Ac)+Cc       , Ac+Cc+Ec))), name="cDZf")
                          
cDZm <-  mxAlgebra( Ly%&%(scalarC %*% (rbind( cbind(Ac+Cc+Ec , (.5%x%Ac)+Cc),
                                cbind((.5%x%Ac)+Cc   , Ac+Cc+Ec)))), name="cDZm")
                          
cDZmf <-  mxAlgebra(Ly%&% (scalarCmf %*% (rbind( cbind(Ac+Cc+Ec , (.5%x%Ac)+Cc),
                                cbind((.5%x%Ac)+Cc    , Ac+Cc+Ec)))), name="cDZmf")      
                                
### Here we create expand (using Ly) the parent-report part of the residual matrix for MZ twins ####
pMZf <-  mxAlgebra(Ly%&%(rbind(cbind(Ap+Cp+Ep   , Ap+Cp),
                          cbind(Ap+Cp     , Ap+Cp+Ep))), name="pMZf")

pMZm <-  mxAlgebra( Ly%&%(scalarP %*% (rbind( cbind(Ap+Cp+Ep , Ap+Cp),
                                cbind(Ap+Cp    , Ap+Cp+Ep)))), name="pMZm")

### Here we create expand (using Ly) the parent-report part of the residual matrix for DZ twins ####
pDZf <-  mxAlgebra(Ly%&%(rbind(cbind( Ap+Cp+Ep , (.5%x%Ap)+Cp),
                          cbind((.5%x%Ap)+Cp        , Ap+Cp+Ep))), name="pDZf")
                          
pDZm <-  mxAlgebra( Ly%&%(scalarP %*% (rbind( cbind(Ap+Cp+Ep , (.5%x%Ap)+Cp),
                                cbind((.5%x%Ap)+Cp   , Ap+Cp+Ep)))), name="pDZm")
                          
pDZmf <-  mxAlgebra( Ly%&%(scalarPmf %*% (rbind( cbind(Ap+Cp+Ep , (.5%x%Ap)+Cp),
                                cbind((.5%x%Ap)+Cp    , Ap+Cp+Ep)))), name="pDZmf")      

### Need to offset and combine parent- and child-rep to create residuals representing rater disagreement

offsetC <- mxMatrix("Full", nrow=ntv, ncol=ntv, free=F, values = 			  c(1,0,1,0,1,0,1,0,1,0,1,0,
																							 0,0,0,0,0,0,0,0,0,0,0,0,
																							 1,0,1,0,1,0,1,0,1,0,1,0,
																							 0,0,0,0,0,0,0,0,0,0,0,0,
																							 1,0,1,0,1,0,1,0,1,0,1,0,
																							 0,0,0,0,0,0,0,0,0,0,0,0,
																							 1,0,1,0,1,0,1,0,1,0,1,0,
																							 0,0,0,0,0,0,0,0,0,0,0,0,
																							 1,0,1,0,1,0,1,0,1,0,1,0,
																							 0,0,0,0,0,0,0,0,0,0,0,0,
																							 1,0,1,0,1,0,1,0,1,0,1,0,
																							 0,0,0,0,0,0,0,0,0,0,0,0), byrow=T, name='offC')

offsetP <- mxMatrix("Full", nrow=ntv, ncol=ntv,free=F, values = 			  c(0,0,0,0,0,0,0,0,0,0,0,0,
																							 0,1,0,1,0,1,0,1,0,1,0,1,
																							 0,0,0,0,0,0,0,0,0,0,0,0,
																							 0,1,0,1,0,1,0,1,0,1,0,1,
																							 0,0,0,0,0,0,0,0,0,0,0,0,
																							 0,1,0,1,0,1,0,1,0,1,0,1,
																							 0,0,0,0,0,0,0,0,0,0,0,0,
																							 0,1,0,1,0,1,0,1,0,1,0,1,
																							 0,0,0,0,0,0,0,0,0,0,0,0,
																							 0,1,0,1,0,1,0,1,0,1,0,1,
																							 0,0,0,0,0,0,0,0,0,0,0,0,
																							 0,1,0,1,0,1,0,1,0,1,0,1), byrow=T, name='offP')

                                                               
### Here we specify the residual matrix for MZ twins ###
ThetaMZf <- mxAlgebra((cMZf*offC)+(pMZf*offP), name="ThetaMZf")
                           
ThetaMZm <- mxAlgebra((cMZm*offC)+(pMZm*offP), name="ThetaMZm")                           
  
### Here we specify the residual matrix for DZ twins ###  
ThetaDZf <- mxAlgebra((cDZf*offC)+(pDZf*offP), name="ThetaDZf")  
                           
ThetaDZm <- mxAlgebra((cDZm*offC)+(pDZm*offP), name="ThetaDZm")             
                           
ThetaDZmf <- mxAlgebra((cDZmf*offC)+(pDZmf*offP), name="ThetaDZmf")             
                           
### Here we specify the total expected covariance, both agreement and disagreement for MZ's and then DZ's ###  
CovMZf <- mxAlgebra( expression=Ly%*%PsyMZf%*%t(Ly)+ThetaMZf,name="CovMZf")
CovMZm <- mxAlgebra( expression=Ly%*%PsyMZm%*%t(Ly)+ThetaMZm,name="CovMZm")
CovDZf <- mxAlgebra( expression=Ly%*%PsyDZf%*%t(Ly)+ThetaDZf,name="CovDZf")
CovDZm <- mxAlgebra( expression=Ly%*%PsyDZm%*%t(Ly)+ThetaDZm,name="CovDZm")                           
CovDZmf <- mxAlgebra( expression=Ly%*%PsyDZmf%*%t(Ly)+ThetaDZmf,name="CovDZmf")                              
                                          
### Here we specify the FIML Objectives for the MZ and DZ model ###
MZfObj <- mxFIMLObjective( covariance="CovMZf", means="MMZf",
                          dimnames=selVars) 
MZmObj <- mxFIMLObjective( covariance="CovMZm", means="MMZm",
                          dimnames=selVars) 
DZfObj <- mxFIMLObjective( covariance="CovDZf", means="MDZf",
                          dimnames=selVars)                            
DZmObj <- mxFIMLObjective( covariance="CovDZm", means="MDZm",
                          dimnames=selVars)        
DZmfObj <- mxFIMLObjective( covariance="CovDZmf", means="MDZmf",
                          dimnames=selVars)                                                       

pars <- list(a,c,e,a2,c2,e2,ap2,cp2,ep2,ac2,cc2,ec2,VA,VC,VE,ap,cp,ep,VAp,VCp,VEp,ac,cc,ec,VAc,VCc,VEc,scalarA,scalarAmf, scalarC,scalarCmf, scalarP, scalarPmf, offsetC, offsetP)                         
MZfModel <-  mxModel("MZfModel",pars,PsyMZf,MeansMZf, MeansMZff, MxdataMZf, CovMZf, MZfObj, Lambda,ThetaMZf, cMZf, pMZf)
MZmModel <-  mxModel("MZmModel",pars,PsyMZm,MeansMZm, MeansMZmm, MxdataMZm, CovMZm, MZmObj, Lambda,ThetaMZm, cMZm, pMZm)
DZfModel <-  mxModel("DZfModel",pars,PsyDZf,MeansDZf, MeansDZff, MxdataDZf, CovDZf, DZfObj, Lambda,ThetaDZf, cDZf, pDZf)
DZmModel <-  mxModel("DZmModel",pars,PsyDZm,MeansDZm, MeansDZmm, MxdataDZm, CovDZm, DZmObj, Lambda,ThetaDZm,cDZm, pDZm)
DZmfModel <-  mxModel("DZmfModel",pars,PsyDZmf,MeansMZf,MeansMZm,MeansDZmf, MxdataDZmf, CovDZmf, DZmfObj, Lambda,ThetaDZmf,cDZmf,pDZmf)
            
CombAlg    <- mxAlgebra(MZfModel.objective +MZmModel.objective +DZfModel.objective+DZmModel.objective+DZmfModel.objective, name="minus2loglikelihood")
CombAlgObj <- mxAlgebraObjective("minus2loglikelihood")      

### Build the final model ###
Psychometric <- mxModel("PsychometricModel",MZfModel,MZmModel,DZfModel,DZmModel,DZmfModel,CombAlg,CombAlgObj)
#Psychometric <- mxModel(Psychometric, mxCI(c("MZfModel.a2",'MZfModel.c2','MZfModel.e2'))) 
#Psychometric <- mxModel(Psychometric, mxCI(c("MZfModel.ap2",'MZfModel.cp2','MZfModel.ep2'))) 
#Psychometric <- mxModel(Psychometric, mxCI(c("MZfModel.ac2",'MZfModel.cc2','MZfModel.ec2'))) 

### Run the model ###
FittedPsych <- mxRun(mxRun(Psychometric))
#FittedPsych <- mxRun(Psychometric, intervals=F)
PsychSumm<- summary(FittedPsych)
source("GenEpiHelperFunctions.R")
parameterSpecifications(FittedPsych)
expectedMeansCovariances(FittedPsych)
tableFitStatistics(FittedPsych)                  

colnames <- c('pdis9_sta2', 'pdis12_sta2', 'pdis14_sta2','pdis9_stc2', 'pdis12_stc2', 'pdis14_stc2','pdis9_ste2', 'pdis12_ste2', 'pdis14_ste2')
rnames <- c( 'pdis9', 'pdis12', 'pdis14')

#Cbind squared a,c,e path estimates
latent.sq.paths <- cbind(FittedPsych$MZfModel.a2@result,FittedPsych$MZfModel.c2@result,FittedPsych$MZfModel.e2@result)
colnames(latent.sq.paths) <- colnames
rownames(latent.sq.paths) <- rnames
latent.sq.paths
parent.sq.paths <- cbind(FittedPsych$MZfModel.ap2@result,FittedPsych$MZfModel.cp2@result,FittedPsych$MZfModel.ep2@result)
colnames(parent.sq.paths) <- colnames
rownames(parent.sq.paths) <- rnames
parent.sq.paths
child.sq.paths <- cbind(FittedPsych$MZfModel.ac2@result,FittedPsych$MZfModel.cc2@result,FittedPsych$MZfModel.ec2@result)
colnames(child.sq.paths) <- colnames
rownames(child.sq.paths) <- rnames
child.sq.paths

#Calculate squared path estmates for child-specific and latent factor variance, standardised to variance of observed child-reports
sq.child.paths <- cbind(latent.sq.paths, child.sq.paths)
sq.child.paths
sq.child.path.sums <- matrix(rep(c((sum(sq.child.paths[1,])),(sum(sq.child.paths[2,])), (sum(sq.child.paths[3,]))), times = 2*(nlv*3)), nrow=3, ncol=(nv*3))
sq.child.path.sums
st.sq.child.paths <- sq.child.paths /sq.child.path.sums
st.sq.child.paths.only<- st.sq.child.paths[,10:18] 
st.sq.clatent.paths.only <- st.sq.child.paths[,1:9]
write.csv(st.sq.child.paths.only,'/Users/lauriejhannigan/Dropbox/Documents (Dropbox)/PhD/My Documents (KCL)/Projects/TEDS/Model fitting/Psychometric/CSVs for plots/pdis_child.csv',row.names=TRUE)
#Standardised path estimates for child-reports only
child <- st.sq.child.paths.only + st.sq.clatent.paths.only

#Contributions of total latent and common phenotypes to variance in child-reports
st.clatent.path.proportion.pheno <- matrix(rep(c((sum(st.sq.clatent.paths.only[1,])),(sum(st.sq.clatent.paths.only[2,])), (sum(st.sq.clatent.paths.only[3,]))), times = 1), nrow=3, ncol=1)
st.clatent.path.proportion.pheno
st.child.path.proportion.pheno <- matrix(rep(c((sum(st.sq.child.paths.only[1,])),(sum(st.sq.child.paths.only[2,])), (sum(st.sq.child.paths.only[3,]))), times = 1), nrow=3, ncol=1)
st.child.path.proportion.pheno

#Calculate squared path estmates for parent-specific and latent factor variance, standardised to variance of observed parent-reports
sq.parent.paths <- cbind(latent.sq.paths, parent.sq.paths)
sq.parent.paths
sq.parent.path.sums <- matrix(rep(c((sum(sq.parent.paths[1,])),(sum(sq.parent.paths[2,])), (sum(sq.parent.paths[3,]))), times = 2*(nlv*3)), nrow=3, ncol=(nv*3))
sq.parent.path.sums
st.sq.parent.paths <-sq.parent.paths /sq.parent.path.sums 
st.sq.parent.paths.only<- st.sq.parent.paths[,10:18] 
st.sq.platent.paths.only <- st.sq.parent.paths[,1:9]
write.csv(st.sq.parent.paths.only,'/Users/lauriejhannigan/Dropbox/Documents (Dropbox)/PhD/My Documents (KCL)/Projects/TEDS/Model fitting/Psychometric/CSVs for plots/pdis_parent.csv',row.names=TRUE)

#Standardised path estimates for parent-reports only
parent <- st.sq.parent.paths.only + st.sq.platent.paths.only
parent

#Contributions of total latent and common phenotypes to variance in parent-reports
st.platent.path.proportion.pheno <- matrix(rep(c((sum(st.sq.platent.paths.only[1,])),(sum(st.sq.platent.paths.only[2,])), (sum(st.sq.platent.paths.only[3,]))), times = 1), nrow=3, ncol=1)
st.platent.path.proportion.pheno
st.parent.path.proportion.pheno <- matrix(rep(c((sum(st.sq.parent.paths.only[1,])),(sum(st.sq.parent.paths.only[2,])), (sum(st.sq.parent.paths.only[3,]))), times = 1), nrow=3, ncol=1)
st.parent.path.proportion.pheno

#Squared path estimates for latent factor standardised to rater agreed variance
latent.sq.paths.sums <- matrix(rep(c((sum(latent.sq.paths[1,])),(sum(latent.sq.paths[2,])), (sum(latent.sq.paths[3,]))), times = 1*(nlv*3)), nrow=3, ncol=(nlv*3))
st.latent.sq.paths <- latent.sq.paths/latent.sq.paths.sums
st.latent.sq.paths

#Squared path estimates for latent factor standardised to variance in parent/child-reports (mean)
unique.sq.paths <- (child.sq.paths + parent.sq.paths)/2
unique.sq.paths
sq.unique.paths <- cbind(latent.sq.paths, unique.sq.paths)
sq.unique.paths
sq.unique.path.sums <- matrix(rep(c((sum(sq.unique.paths[1,])),(sum(sq.unique.paths[2,])), (sum(sq.unique.paths[3,]))), times = 2*(nlv*3)), nrow=3, ncol=(nv*3))
sq.unique.path.sums
st.sq.unique.paths <- (sq.unique.paths /sq.unique.path.sums) 
st.sq.unique.paths.only<- st.sq.unique.paths[,10:18] 
st.sq.ulatent.paths.only <- st.sq.unique.paths[,1:9]
st.sq.ulatent.paths.only
write.csv(st.sq.ulatent.paths.only,'/Users/lauriejhannigan/Dropbox/Documents (Dropbox)/PhD/My Documents (KCL)/Projects/TEDS/Model fitting/Psychometric/CSVs for plots/pdis_latent.csv',row.names=TRUE)


#sink("/Users/lauriejhannigan/Dropbox/Documents (Dropbox)/PhD/My Documents (KCL)/Projects/TEDS/Model fitting/Psychometric/psychometric_pdis.csv")
fits <- tableFitStatistics(FittedPsych)
output <- list(latent.sq.paths, parent.sq.paths, child.sq.paths, parent, child, fits)
output
confs <- PsychSumm$CI
confs
#sink()
#sink()

########################## PLOTS ##########################

library(devtools)
source_url('https://gist.github.com/menugget/7689145/raw/dac746aa322ca4160a5fe66c70fec16ebe26faf9/image.scale.2.r')
source_url('https://gist.github.com/menugget/7864454/raw/f698da873766347d837865eecfa726cdf52a6c40/plot.stream.4.R')
source_url('https://gist.github.com/menugget/7864471/raw/8127dfaae183233d203580bc247a73a564d5feab/plot.stacked.2.R')

########## STACKED#########

#plot.stacked makes a stacked plot where each y series is plotted on top
#of the each other using filled polygons
#
#Arguments include:
#'x' - a vector of values
#'y' - a matrix of data series (columns) corresponding to x
#'order.method' = c("as.is", "max", "first") 
#  "as.is" - plot in order of y column
#  "max" - plot in order of when each y series reaches maximum value
#  "first" - plot in order of when each y series first value > 0
#'col' - fill colors for polygons corresponding to y columns (will recycle)
#'border' - border colors for polygons corresponding to y columns (will recycle) (see ?polygon for details)
#'lwd' - border line width for polygons corresponding to y columns (will recycle)
#'...' - other plot arguments

source("plot_stacked.R") 

st.sq.ulatent.paths.only
 y <- st.sq.ulatent.paths.only
 z <- st.sq.unique.paths.only
 y1 <- y[1:3,1:3]
 y2 <- y[1:3,4:6]
 y3 <- y[1:3,7:9]
 y4 <- matrix(rep(c((sum(z[1,])),(sum(z[2,])), (sum(z[3,]))), times = 1, nrow=3, ncol=1))
 x <- c(9,12,14)
#Plot Ex. 1 - Color by max value
#pal <- colorRampPalette(c(rgb(0.85,0.85,1), rgb(0.2,0.2,0.7)))
#BREAKS <- pretty(apply(y,2,max),8)
#LEVS <- levels(cut(1, breaks=BREAKS))
#COLS <- pal(length(BREAKS )-1)

 atot <- matrix(c((sum(y1[1,1:3])),(sum(y1[2,1:3])),(sum(y1[3,1:3]))), nrow=3,ncol=1)
 ctot <- matrix(c((sum(y2[1,1:3])),(sum(y2[2,1:3])),(sum(y2[3,1:3]))), nrow=3,ncol=1)
 etot <- matrix(c((sum(y3[1,1:3])),(sum(y3[2,1:3])),(sum(y3[3,1:3]))), nrow=3,ncol=1)
 y11 <- cbind(y1,ctot,etot,y4)
 y22 <- cbind(atot, y2,etot,y4)
 y33 <- cbind(atot,ctot, y3,y4)
yall<-cbind(y1,y2,y3,y4)
yall
y11
y22
y33


png("/Users/lauriejhannigan/Dropbox/Documents (Dropbox)/PhD/My Documents (KCL)/Projects/TEDS/Plots/Psychometric/pdis_latent.png", width=15, height=4.5, units="in", res=400)
layout(matrix(1:6, nrow=2, ncol=3), widths=c(4.5,4.5,4.5,4.5,4.5,4.5), heights=c(4,.5,4,.5,4,.5), respect=T)
par(mar=c(0,2.5,2,2.5), cex=0.75, xpd=F)
plot.stacked(x,y11, xlim=c(9,14), ylim=c(0, 1), yaxs="i", col=c('navy','royalblue3','lightskyblue', 'gray87', 'gray87','white'), border='gray0', lwd=0.5, xaxt="n" ,yaxp=c(0,1,10), ylab="Proportion of phenotypic variance accounted for", las=2)
axis(1, at = c(9,12,14), labels= c('Age 9','Age 12','Age 14') )
mtext("Genetic", line=0.25, side=3)
mtext("Proportion of phenotypic variance accounted for", line=3, side=2)
box()
par(mar=c(1,1,1,2.5), cex=1)
plot(1, t="n", xlab="", ylab="", axes=FALSE)
par(mar=c(0,2.5,2,2.5), cex=0.75)
plot.stacked(x,y22, xlim=c(9,14), ylim=c(0, 1), yaxs="i", col=c('gray87','darkgreen','forestgreen','darkolivegreen1', 'gray87','white'), border="grey0", lwd=0.5, xaxt="n",yaxp=c(0,1,10), las=2)
axis(1, at = c(9,12,14), labels= c('Age 9','Age 12','Age 14') )
mtext("Shared environmental", line=0.25, side=3)
box()
par(mar=c(1,1,1,2.5), cex=1)
plot(1, t="n", xlab="", ylab="", axes=FALSE)
par(mar=c(0,2.5,2,2.5), cex=0.75)
plot.stacked(x,y33, xlim=c(9,14), ylim=c(0, 1), yaxs="i", col=c('gray87', 'gray87','darkred','firebrick2','darksalmon','white'), border="grey0", lwd=0.5, xaxt="n",yaxp=c(0,1,10), las=2)
axis(1, at = c(9,12,14), labels= c('Age 9','Age 12','Age 14') )
mtext("Unique environmental", line=0.25, side=3)
box()
par(mar=c(1,1,1,2.5), cex=1)
plot(1, t="n", xlab="", ylab="", axes=FALSE)
dev.off()

#Legend only A
png("/Users/lauriejhannigan/Dropbox/Documents (Dropbox)/PhD/My Documents (KCL)/Projects/TEDS/Plots/Psychometric/A_legend_common.png", width=15, height=12, units="cm", res=400)
layout(matrix(1:1, nrow=1, ncol=1), widths=lcm(10), heights=lcm(c(11)), respect=F)
par(mar=c(0,0,0,0))
plot(1, t="n", xlab="", ylab="", axes=FALSE)
legend("center", legend=c('Rater disagreement','Age 9 genetic factors', 'Age 12 genetic factors', 'Age 14 genetic factors' ), border="grey0", fill=c('white','navy','royalblue3','lightskyblue')) #pch=22, pt.bg=COL)
dev.off()

#Legend only C
png("/Users/lauriejhannigan/Dropbox/Documents (Dropbox)/PhD/My Documents (KCL)/Projects/TEDS/Plots/Psychometric/C_legend_common.png", width=15, height=12, units="cm", res=400)
layout(matrix(1:1, nrow=1, ncol=1), widths=lcm(10), heights=lcm(c(11)), respect=F)
par(mar=c(0,0,0,0))
plot(1, t="n", xlab="", ylab="", axes=FALSE)
legend("center", legend=c('Rater disagreement','Age 9 shared environmental factors', 'Age 12 shared environmental factors', 'Age 14 shared environmental factors'), border="grey0", fill=c('white','darkgreen','forestgreen','darkolivegreen1')) #pch=22, pt.bg=COL)
dev.off()

#Legend only E
png("/Users/lauriejhannigan/Dropbox/Documents (Dropbox)/PhD/My Documents (KCL)/Projects/TEDS/Plots/Psychometric/E_legend_common.png", width=15, height=12, units="cm", res=400)
layout(matrix(1:1, nrow=1, ncol=1), widths=lcm(10), heights=lcm(c(11)), respect=F)
par(mar=c(0,0,0,0))
plot(1, t="n", xlab="", ylab="", axes=FALSE)
legend("center",  legend=c('Rater disagreement','Age 9 unique environmental factors', 'Age 12 unique environmental factors', 'Age 14 unique environmental factors'), border="grey0", fill=c('white','darkred','firebrick2','darksalmon')) #pch=22, pt.bg=COL)
dev.off()
#############################


 y <- st.sq.parent.paths.only
 z <- st.sq.platent.paths.only
 y1 <- y[1:3,1:3]
 y2 <- y[1:3,4:6]
 y3 <- y[1:3,7:9]
 y4 <- matrix(rep(c((sum(z[1,])),(sum(z[2,])), (sum(z[3,]))), times = 1, nrow=3, ncol=1))
 x <- c(9,12,14)
#Plot Ex. 1 - Color by max value
#pal <- colorRampPalette(c(rgb(0.85,0.85,1), rgb(0.2,0.2,0.7)))
#BREAKS <- pretty(apply(y,2,max),8)
#LEVS <- levels(cut(1, breaks=BREAKS))
#COLS <- pal(length(BREAKS )-1)

 atot <- matrix(c((sum(y1[1,1:3])),(sum(y1[2,1:3])),(sum(y1[3,1:3]))), nrow=3,ncol=1)
 ctot <- matrix(c((sum(y2[1,1:3])),(sum(y2[2,1:3])),(sum(y2[3,1:3]))), nrow=3,ncol=1)
 etot <- matrix(c((sum(y3[1,1:3])),(sum(y3[2,1:3])),(sum(y3[3,1:3]))), nrow=3,ncol=1)
 y11 <- cbind(y4,y1,ctot,etot)
 y22 <- cbind(y4,atot, y2,etot)
 y33 <- cbind(y4,atot,ctot, y3)
yall<-cbind(y4,y1,y2,y3)
yall
y11
y22
y33


png("/Users/lauriejhannigan/Dropbox/Documents (Dropbox)/PhD/My Documents (KCL)/Projects/TEDS/Plots/Psychometric/pdis_parent.png", width=15, height=4.5, units="in", res=400)
layout(matrix(1:6, nrow=2, ncol=3), widths=c(4.5,4.5,4.5,4.5,4.5,4.5), heights=c(4,.5,4,.5,4,.5), respect=T)
par(mar=c(0,2.5,2,2.5), cex=0.75)
plot.stacked(x,y11, xlim=c(9,14), ylim=c(0, 1), yaxs="i", col=c('white','navy','royalblue3','lightskyblue', 'gray87', 'gray87'), border='gray0', lwd=0.5, xaxt="n" ,yaxp=c(0,1,10), ylab="Proportion of phenotypic variance accounted for", las=2)
axis(1, at = c(9,12,14), labels= c('Age 9','Age 12','Age 14') )
mtext("Genetic", line=0.25, side=3)
mtext("Proportion of phenotypic variance accounted for", line=3, side=2)
box()
par(mar=c(1,1,1,2.5), cex=1)
plot(1, t="n", xlab="", ylab="", axes=FALSE)
par(mar=c(0,2.5,2,2.5), cex=0.75)
plot.stacked(x,y22, xlim=c(9,14), ylim=c(0, 1), yaxs="i", col=c('white','gray87','darkgreen','forestgreen','darkolivegreen1', 'gray87'), border="grey0", lwd=0.5, xaxt="n",yaxp=c(0,1,10), las=2)
axis(1, at = c(9,12,14), labels= c('Age 9','Age 12','Age 14') )
mtext("Shared environmental", line=0.25, side=3)
box()
par(mar=c(1,1,1,2.5), cex=1)
plot(1, t="n", xlab="", ylab="", axes=FALSE)
par(mar=c(0,2.5,2,2.5), cex=0.75)
plot.stacked(x,y33, xlim=c(9,14), ylim=c(0, 1), yaxs="i", col=c('white','gray87', 'gray87','darkred','firebrick2','darksalmon'), border="grey0", lwd=0.5, xaxt="n",yaxp=c(0,1,10), las=2)
axis(1, at = c(9,12,14), labels= c('Age 9','Age 12','Age 14') )
mtext("Unique environmental", line=0.25, side=3)
box()
par(mar=c(1,1,1,2.5), cex=1)
plot(1, t="n", xlab="", ylab="", axes=FALSE)
dev.off()

#Legend only A
png("/Users/lauriejhannigan/Dropbox/Documents (Dropbox)/PhD/My Documents (KCL)/Projects/TEDS/Plots/Psychometric/A_legend_unique.png", width=15, height=12, units="cm", res=400)
layout(matrix(1:1, nrow=1, ncol=1), widths=lcm(10), heights=lcm(c(11)), respect=F)
par(mar=c(0,0,0,0))
plot(1, t="n", xlab="", ylab="", axes=FALSE)
legend("center", legend=c('Rater agreement','Age 9 genetic factors', 'Age 12 genetic factors', 'Age 14 genetic factors' ), border="grey0", fill=c('white','navy','royalblue3','lightskyblue')) #pch=22, pt.bg=COL)
dev.off()

#Legend only C
png("/Users/lauriejhannigan/Dropbox/Documents (Dropbox)/PhD/My Documents (KCL)/Projects/TEDS/Plots/Psychometric/C_legend_unique.png", width=15, height=12, units="cm", res=400)
layout(matrix(1:1, nrow=1, ncol=1), widths=lcm(10), heights=lcm(c(11)), respect=F)
par(mar=c(0,0,0,0))
plot(1, t="n", xlab="", ylab="", axes=FALSE)
legend("center", legend=c('Rater agreement','Age 9 shared environmental factors', 'Age 12 shared environmental factors', 'Age 14 shared environmental factors'), border="grey0", fill=c('white','darkgreen','forestgreen','darkolivegreen1')) #pch=22, pt.bg=COL)
dev.off()

#Legend only E
png("/Users/lauriejhannigan/Dropbox/Documents (Dropbox)/PhD/My Documents (KCL)/Projects/TEDS/Plots/Psychometric/E_legend_unique.png", width=15, height=12, units="cm", res=400)
layout(matrix(1:1, nrow=1, ncol=1), widths=lcm(10), heights=lcm(c(11)), respect=F)
par(mar=c(0,0,0,0))
plot(1, t="n", xlab="", ylab="", axes=FALSE)
legend("center",  legend=c('Rater agreement','Age 9 unique environmental factors', 'Age 12 unique environmental factors', 'Age 14 unique environmental factors'), border="grey0", fill=c('white','darkred','firebrick2','darksalmon')) #pch=22, pt.bg=COL)
dev.off()

##############


 y <- st.sq.child.paths.only
 z <- st.sq.clatent.paths.only
 y1 <- y[1:3,1:3]
 y2 <- y[1:3,4:6]
 y3 <- y[1:3,7:9]
 y4 <- matrix(rep(c((sum(z[1,])),(sum(z[2,])), (sum(z[3,]))), times = 1, nrow=3, ncol=1))
 x <- c(9,12,14)
#Plot Ex. 1 - Color by max value
#pal <- colorRampPalette(c(rgb(0.85,0.85,1), rgb(0.2,0.2,0.7)))
#BREAKS <- pretty(apply(y,2,max),8)
#LEVS <- levels(cut(1, breaks=BREAKS))
#COLS <- pal(length(BREAKS )-1)

 atot <- matrix(c((sum(y1[1,1:3])),(sum(y1[2,1:3])),(sum(y1[3,1:3]))), nrow=3,ncol=1)
 ctot <- matrix(c((sum(y2[1,1:3])),(sum(y2[2,1:3])),(sum(y2[3,1:3]))), nrow=3,ncol=1)
 etot <- matrix(c((sum(y3[1,1:3])),(sum(y3[2,1:3])),(sum(y3[3,1:3]))), nrow=3,ncol=1)
 y11 <- cbind(y4,y1,ctot,etot)
 y22 <- cbind(y4,atot, y2,etot)
 y33 <- cbind(y4,atot,ctot, y3)
yall<-cbind(y4,y1,y2,y3)
yall
y11
y22
y33


png("/Users/lauriejhannigan/Dropbox/Documents (Dropbox)/PhD/My Documents (KCL)/Projects/TEDS/Plots/Psychometric/pdis_child.png", width=15, height=4.5, units="in", res=400)
layout(matrix(1:6, nrow=2, ncol=3), widths=c(4.5,4.5,4.5,4.5,4.5,4.5), heights=c(4,.5,4,.5,4,.5), respect=T)
par(mar=c(0,2.5,2,2.5), cex=0.75)
plot.stacked(x,y11, xlim=c(9,14), ylim=c(0, 1), yaxs="i", col=c('white','navy','royalblue3','lightskyblue', 'gray87', 'gray87'), border='gray0', lwd=0.5, xaxt="n" ,yaxp=c(0,1,10), ylab="", las=2)
axis(1, at = c(9,12,14), labels= c('Age 9','Age 12','Age 14') )
mtext("Genetic", line=0.25, side=3)
mtext("Proportion of phenotypic variance accounted for", line=3, side=2)
box()
par(mar=c(1,1,1,2.5), cex=1)
plot(1, t="n", xlab="", ylab="", axes=FALSE)
par(mar=c(0,2.5,2,2.5), cex=0.75)
plot.stacked(x,y22, xlim=c(9,14), ylim=c(0, 1), yaxs="i", col=c('white','gray87','darkgreen','forestgreen','darkolivegreen1', 'gray87'), border="grey0", lwd=0.5, xaxt="n",yaxp=c(0,1,10), las=2)
axis(1, at = c(9,12,14), labels= c('Age 9','Age 12','Age 14') )
mtext("Shared environmental", line=0.25, side=3)
box()
par(mar=c(1,1,1,2.5), cex=1)
plot(1, t="n", xlab="", ylab="", axes=FALSE)
par(mar=c(0,2.5,2,2.5), cex=0.75)
plot.stacked(x,y33, xlim=c(9,14), ylim=c(0, 1), yaxs="i", col=c('white','gray87', 'gray87','darkred','firebrick2','darksalmon'), border="grey0", lwd=0.5, xaxt="n",yaxp=c(0,1,10), las=2)
axis(1, at = c(9,12,14), labels= c('Age 9','Age 12','Age 14') )
mtext("Unique environmental", line=0.25, side=3)
box()
par(mar=c(1,1,1,2.5), cex=1)
plot(1, t="n", xlab="", ylab="", axes=FALSE)
dev.off()